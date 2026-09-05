import { classifyContaminationAsync } from "./contamination";
import { runFastTreeIsolated } from "./biowasm";
import { runAlivibeMsa } from "./alivibe-msa-runtime";
import { prepareConsensusThresholdReanalysis, restoreStoredPipelineConfig, restoreUntouchedThresholdStatuses } from "./consensus-threshold-reanalysis";
import { collapsePostprocess, postprocess, type MsaRunner, type PanelFilterRunner, type PostprocessOutput } from "./postprocess";
import { OPTIONAL_STAGE_ORDER, markOptionalStageSkipped, stageCompleted, stageStatus, statusRecord } from "./optional-stages";
import { treeTipNames } from "./tree-names";
import type { OptionalStageName, PipelineTiming, ResultBundle, ThresholdSelection } from "./types";

export { buildStoredConsensusThresholdReview, prepareConsensusThresholdReanalysis } from "./consensus-threshold-reanalysis";

export interface DeferredAnalysisProgress { stage: OptionalStageName; fraction: number; detail: string }
export interface DeferredAnalysisOptions {
  signal?: AbortSignal;
  onProgress?(progress: DeferredAnalysisProgress): void;
  onCheckpoint?(bundle: ResultBundle): void;
  /** Internal sample-scoped replay controls. */
  sampleNames?: readonly string[];
  forceRecompute?: boolean;
  /** Injectable engines keep the replay path independently testable. */
  runMsa?: MsaRunner;
  panelMsa?: MsaRunner;
  panelFilter?: PanelFilterRunner;
  runTree?(fasta: string, signal?: AbortSignal): Promise<string>;
}

function explicitStatuses(bundle: ResultBundle): NonNullable<ResultBundle["optionalStages"]> {
  return Object.fromEntries(OPTIONAL_STAGE_ORDER.map((stage) => [stage, { ...stageStatus(bundle, stage) }])) as NonNullable<ResultBundle["optionalStages"]>;
}

function timing(timings: PipelineTiming[] | undefined, stage: PipelineTiming["stage"], seconds: number, workItems?: number) {
  const entry: PipelineTiming = { stage, seconds }; if (workItems != null) entry.workItems = workItems;
  return [...(timings ?? []).filter((row) => row.stage !== stage), entry];
}

function starTree(fasta: string) {
  const names = treeTipNames([...fasta.matchAll(/^>(.*)$/gm)].map((match) => match[1]));
  return `(${names.map((name) => `${name}:0.0`).join(",")});`;
}

function report(options: DeferredAnalysisOptions, stage: OptionalStageName, fraction: number, detail: string) {
  options.onProgress?.({ stage, fraction: Math.max(0, Math.min(1, fraction)), detail });
}

function checkpoint(options: DeferredAnalysisOptions, bundle: ResultBundle) { options.onCheckpoint?.(bundle); return bundle; }

function omitSampleKeys<T>(source: Record<string, T> | undefined, samples: ReadonlySet<string>): Record<string, T> {
  return Object.fromEntries(Object.entries(source ?? {}).filter(([key]) => !samples.has(key.split("/", 1)[0])));
}

function mergeSampleRows<T extends { sample: string }>(original: readonly T[], replacement: readonly T[], samples: ReadonlySet<string>, order: readonly string[]): T[] {
  const sampleOrder = new Map(order.map((name, index) => [name, index]));
  return [...original.filter((row) => !samples.has(row.sample)), ...replacement]
    .sort((left, right) => (sampleOrder.get(left.sample) ?? Number.MAX_SAFE_INTEGER) - (sampleOrder.get(right.sample) ?? Number.MAX_SAFE_INTEGER));
}

/** Apply new consensus thresholds and recompute through the prior completion boundary. */
export async function recomputeAfterConsensusThresholds(bundle: ResultBundle, selection: ThresholdSelection,
  options: DeferredAnalysisOptions = {}): Promise<ResultBundle> {
  const plan = prepareConsensusThresholdReanalysis(bundle, selection);
  if (!plan.target) return checkpoint(options, plan.bundle);
  if (plan.scope === "samples") {
    // Keep the public bundle valid and unchanged until every requested stage
    // has completed for the affected samples. This makes cancellation atomic.
    const computed = await computeThrough(plan.bundle, plan.target, { ...options,
      sampleNames: plan.affectedSamples, forceRecompute: true, onCheckpoint: undefined });
    return checkpoint(options, computed);
  }
  const wrapped: DeferredAnalysisOptions = { ...options,
    onCheckpoint: (value) => options.onCheckpoint?.(restoreUntouchedThresholdStatuses(value, plan)) };
  return restoreUntouchedThresholdStatuses(await computeThrough(plan.bundle, plan.target, wrapped), plan);
}

export { markOptionalStageSkipped };

/** Compute a requested optional output plus every missing prerequisite. */
export async function computeThrough(bundle: ResultBundle, target: OptionalStageName, options: DeferredAnalysisOptions = {}): Promise<ResultBundle> {
  const targetIndex = OPTIONAL_STAGE_ORDER.indexOf(target); if (targetIndex < 0) throw new Error(`Unknown optional stage: ${target}`);
  const selectedSamples = options.sampleNames ? new Set(options.sampleNames) : undefined;
  const sampleScoped = Boolean(selectedSamples), force = Boolean(options.forceRecompute);
  const contaminationReapply = targetIndex >= 1 && stageCompleted(bundle, "postprocessing")
    && bundle.postprocessingContaminationMode === "bypassed" && stageCompleted(bundle, "contamination")
    && bundle.config.parameters.contaminationFilter;
  if (targetIndex >= 1 && (!stageCompleted(bundle, "postprocessing") || contaminationReapply || force) && !bundle.downstreamResources)
    throw new Error("This project does not contain the panel/reference sequences required to resume downstream filtering. Re-run from the FASTQ with webPORPID 0.3.6 or later.");
  const config = restoreStoredPipelineConfig(bundle);
  const runMsa = options.runMsa ?? runAlivibeMsa, runTree = options.runTree ?? runFastTreeIsolated;
  if (selectedSamples) for (const sample of selectedSamples)
    if (!config.samples.some((candidate) => candidate.name === sample)) throw new Error(`Cannot recompute unknown sample ${sample}.`);
  let current: ResultBundle = { ...bundle, optionalStages: explicitStatuses(bundle) };

  if (target === "contamination" && !stageCompleted(current, "contamination")) {
    const started = performance.now(); report(options, "contamination", 0, "Preparing run-wide contamination checks");
    const contamination = await classifyContaminationAsync(current.consensuses, config, options.signal, (state) =>
      report(options, "contamination", state.fraction, state.detail));
    const summaries = current.summaries.map((summary) => ({ ...summary,
      contaminationPassed: summary.consensusSequences - contamination.filter((call) => call.sample === summary.sample && call.discarded).length }));
    const stamp = new Date().toISOString(), statuses = explicitStatuses(current);
    statuses.contamination = statusRecord("completed", `${contamination.filter((call) => call.discarded).length} consensus sequences identified for exclusion.`, stamp);
    if (stageCompleted(current, "postprocessing") && current.postprocessingContaminationMode === "bypassed")
      statuses.postprocessing = statusRecord("completed", `${statuses.postprocessing.detail} Contamination decisions are now available, but these stored downstream outputs remain explicitly unfiltered until reapplied.`, stamp);
    current = checkpoint(options, { ...current, contamination, summaries, optionalStages: statuses,
      timings: timing(current.timings, "contamination", (performance.now() - started) / 1000, current.consensuses.length),
      log: [...current.log, `${stamp} on-demand contamination: completed; ${contamination.filter((call) => call.discarded).length} discarded`] });
  }
  const shouldReapply = targetIndex >= 1 && stageCompleted(current, "postprocessing")
    && current.postprocessingContaminationMode === "bypassed" && stageCompleted(current, "contamination")
    && config.parameters.contaminationFilter;
  if (targetIndex >= 1 && (!stageCompleted(current, "postprocessing") || shouldReapply || force)) {
    const contaminationApplied = stageCompleted(current, "contamination") && config.parameters.contaminationFilter;
    const started = performance.now(); report(options, "postprocessing", 0, `Starting panel screening, retained-family alignments, and annotations${contaminationApplied ? " with contamination decisions applied" : " without contamination filtering"}`);
    const downstream = await postprocess(current.consensuses, contaminationApplied ? current.contamination : [], config, options.signal, runMsa,
      Math.max(1, current.provenance.workers), (state) => report(options, "postprocessing", state.fraction, state.detail),
      { collapse: false, sampleNames: selectedSamples, panelMsa: options.panelMsa, panelFilter: options.panelFilter });
    const baseBySample = new Map(current.summaries.map((summary) => [summary.sample, summary]));
    downstream.summaries.forEach((summary) => {
      const base = baseBySample.get(summary.sample); if (!base) throw new Error(`Stored summary is missing sample ${summary.sample}.`);
      summary.demultiplexedReads = base.demultiplexedReads;
      summary.selectedReads = base.selectedReads; summary.downsampledReads = base.downsampledReads;
      summary.observedUmis = base.observedUmis; summary.likelyRealUmis = base.likelyRealUmis;
      if (!contaminationApplied) delete summary.contaminationPassed;
    });
    const stamp = new Date().toISOString(), statuses = explicitStatuses(current);
    const seconds = (performance.now() - started) / 1000;
    statuses.postprocessing = statusRecord("completed", sampleScoped
      ? `Recomputed ${downstream.records.length} consensus-family records for ${selectedSamples!.size} changed sample(s); unchanged sample outputs were preserved.`
      : `${downstream.records.length} consensus-family records evaluated${contaminationApplied ? " with contamination decisions applied" : "; contamination was bypassed and excluded zero sequences"}.`, stamp);
    if (!sampleScoped) {
      statuses.collapse = statusRecord("deferred", "Waiting for the current downstream-filtering output.", stamp);
      statuses.tree = statusRecord("deferred", "Waiting for haplotype collapse.", stamp);
    }
    const records = sampleScoped ? mergeSampleRows(current.records, downstream.records, selectedSamples!, config.samples.map((sample) => sample.name)) : downstream.records;
    const summaries = sampleScoped ? mergeSampleRows(current.summaries, downstream.summaries, selectedSamples!, config.samples.map((sample) => sample.name)) : downstream.summaries;
    const alignments = sampleScoped ? { ...omitSampleKeys(current.alignments, selectedSamples!), ...downstream.alignments } : downstream.alignments;
    const referenceAlignments = sampleScoped ? { ...omitSampleKeys(current.referenceAlignments, selectedSamples!), ...downstream.referenceAlignments } : downstream.referenceAlignments;
    const collapseGroups = sampleScoped ? { ...Object.fromEntries(Object.entries(current.collapseGroups ?? {}).filter(([sample]) => !selectedSamples!.has(sample))) } : {};
    const functionalFilterErrors = sampleScoped
      ? Object.fromEntries(Object.entries(current.functionalFilterErrors ?? {}).filter(([sample]) => !selectedSamples!.has(sample))) : {};
    const trees = sampleScoped ? omitSampleKeys(current.trees, selectedSamples!) : {};
    const alignmentEdits = sampleScoped ? omitSampleKeys(current.alignmentEdits, selectedSamples!) : {};
    current = checkpoint(options, { ...current, records, summaries, alignments, referenceAlignments, collapseGroups,
      functionalFilterErrors: Object.keys(functionalFilterErrors).length ? functionalFilterErrors : undefined,
      trees, alignmentEdits, optionalStages: statuses, postprocessingContaminationMode: contaminationApplied ? "applied" : "bypassed",
      timings: sampleScoped ? current.timings : timing(current.timings, "postprocessing", seconds, downstream.records.length),
      log: [...current.log, `${stamp} on-demand postprocessing: completed for ${downstream.records.length} consensus families${sampleScoped ? ` in changed samples ${[...selectedSamples!].join(", ")}; ${seconds.toFixed(3)} s; unchanged sample products and edits preserved` : ""}; contamination=${contaminationApplied ? "applied" : "bypassed"}${shouldReapply ? "; prior alignment edits detached because their baseline was replaced" : ""}`] });
  }
  if (targetIndex >= 2 && (!stageCompleted(current, "collapse") || force)) {
    const started = performance.now(), input: PostprocessOutput = { records: current.records, summaries: current.summaries,
      alignments: current.alignments, referenceAlignments: current.referenceAlignments ?? {}, collapseGroups: current.collapseGroups ?? {},
      functionalFilterErrors: current.functionalFilterErrors ?? {}, collapseSeconds: 0 };
    const collapsed = await collapsePostprocess(input, config, options.signal,
      (state) => report(options, "collapse", state.fraction, state.detail), runMsa, selectedSamples);
    const count = Object.values(collapsed.collapseGroups).reduce((sum, groups) => sum + groups.length, 0), stamp = new Date().toISOString();
    const functionalCount = collapsed.summaries.reduce((sum, summary) => sum + (summary.functionalPassed ?? 0), 0);
    const selectedCount = selectedSamples ? [...selectedSamples].reduce((sum, sample) => sum + (collapsed.collapseGroups[sample]?.length ?? 0), 0) : count;
    const seconds = (performance.now() - started) / 1000, statuses = explicitStatuses(current);
    statuses.collapse = statusRecord("completed", sampleScoped
      ? `${selectedCount} haplotypes recomputed for ${selectedSamples!.size} changed sample(s); unchanged collapse/functional outputs were preserved.`
      : `${count} haplotypes; ${functionalCount} collapsed variants passed configured functional filters; multiplicities count UMI families.`, stamp);
    const errorLog = Object.entries(collapsed.functionalFilterErrors).filter(([sample]) => !selectedSamples || selectedSamples.has(sample))
      .map(([sample, message]) => `${stamp} functional filter error: ${sample}; ${message}; this sample's functional outputs were omitted and other samples continued`);
    current = checkpoint(options, { ...current, summaries: collapsed.summaries, alignments: collapsed.alignments,
      referenceAlignments: collapsed.referenceAlignments, collapseGroups: collapsed.collapseGroups,
      functionalFilterErrors: Object.keys(collapsed.functionalFilterErrors).length ? collapsed.functionalFilterErrors : undefined,
      optionalStages: statuses, timings: sampleScoped ? current.timings : timing(current.timings, "collapse", seconds, count),
      log: [...current.log, `${stamp} on-demand collapse: ${sampleScoped ? `${selectedCount} haplotypes in changed samples ${[...selectedSamples!].join(", ")}; ${seconds.toFixed(3)} s; unchanged samples preserved` : `${count} haplotypes; ${functionalCount} collapsed functional passes`}; multiplicities count UMI families`, ...errorLog] });
  }
  if (targetIndex >= 3 && (!stageCompleted(current, "tree") || force)) {
    const inputs = config.samples.filter((sample) => !selectedSamples || selectedSamples.has(sample.name)).flatMap((sample) => {
      const key = `${sample.name}/nucleotide`, fasta = current.alignments[key]; return fasta ? [[key, fasta] as const] : [];
    });
    const trees = { ...current.trees }, started = performance.now(); let cursor = 0, finished = 0;
    await Promise.all(Array.from({ length: Math.min(Math.max(1, current.provenance.workers), Math.max(1, inputs.length)) }, async () => {
      while (true) {
        if (options.signal?.aborted) throw new DOMException("Phylogeny inference skipped.", "AbortError");
        const index = cursor++; if (index >= inputs.length) return;
        const [key, fasta] = inputs[index]; report(options, "tree", finished / Math.max(1, inputs.length),
          `Inferring ${key.split("/")[0]} (${finished} of ${inputs.length} complete) with bundled double-precision FastTree`);
        try { trees[key] = await runTree(fasta, options.signal); }
        catch (cause) { if (options.signal?.aborted) throw cause; trees[key] = starTree(fasta); }
        finished++; report(options, "tree", finished / Math.max(1, inputs.length), `Finished ${finished} of ${inputs.length} collapsed phylogenies`);
      }
    }));
    const stamp = new Date().toISOString(), seconds = (performance.now() - started) / 1000, statuses = explicitStatuses(current);
    statuses.tree = statusRecord("completed", sampleScoped
      ? `${inputs.length} changed-sample collapsed phylogenies recomputed; unchanged trees were preserved.`
      : inputs.length ? `${inputs.length} collapsed phylogenies inferred.` : "No retained collapsed alignments required a tree.", stamp);
    current = checkpoint(options, { ...current, trees, optionalStages: statuses,
      timings: sampleScoped ? current.timings : timing(current.timings, "tree", seconds, inputs.length),
      log: [...current.log, `${stamp} on-demand phylogeny: ${inputs.length} collapsed trees inferred with bundled double-precision FastTree${sampleScoped ? ` for changed samples ${[...selectedSamples!].join(", ")}; ${seconds.toFixed(3)} s; unchanged trees preserved` : ""}`] });
  }
  return current;
}
