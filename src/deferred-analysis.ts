import { classifyContaminationAsync } from "./contamination";
import { runFastTreeIsolated } from "./biowasm";
import { runAlivibeMsa } from "./alivibe-msa-runtime";
import { collapsePostprocess, postprocess, type PostprocessOutput } from "./postprocess";
import { OPTIONAL_STAGE_ORDER, restorePipelineConfig, stageCompleted, stageStatus, statusRecord } from "./optional-stages";
import { treeTipNames } from "./tree-names";
import type { OptionalStageName, PipelineConfig, PipelineTiming, ResultBundle } from "./types";

export interface DeferredAnalysisProgress { stage: OptionalStageName; fraction: number; detail: string }
export interface DeferredAnalysisOptions {
  signal?: AbortSignal;
  onProgress?(progress: DeferredAnalysisProgress): void;
  onCheckpoint?(bundle: ResultBundle): void;
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

export function markOptionalStageSkipped(bundle: ResultBundle, stage: OptionalStageName): ResultBundle {
  const statuses = explicitStatuses(bundle), index = OPTIONAL_STAGE_ORDER.indexOf(stage), stamp = new Date().toISOString();
  statuses[stage] = statusRecord("skipped", "Skipped by user during on-demand computation; partial outputs from this stage were discarded.", stamp);
  for (const downstream of OPTIONAL_STAGE_ORDER.slice(index + 1))
    statuses[downstream] = statusRecord("deferred", `Waiting for ${stage} after it was skipped.`, stamp);
  return { ...bundle, optionalStages: statuses,
    log: [...bundle.log, `${stamp} on-demand ${stage}: skipped by user; partial stage outputs discarded`] };
}

/** Compute a requested optional output plus every missing prerequisite. */
export async function computeThrough(bundle: ResultBundle, target: OptionalStageName, options: DeferredAnalysisOptions = {}): Promise<ResultBundle> {
  const targetIndex = OPTIONAL_STAGE_ORDER.indexOf(target); if (targetIndex < 0) throw new Error(`Unknown optional stage: ${target}`);
  if (targetIndex >= 1 && !stageCompleted(bundle, "postprocessing") && !bundle.downstreamResources)
    throw new Error("This project does not contain the panel/reference sequences required to resume downstream filtering. Re-run from the FASTQ with webPORPID 0.3.6 or later.");
  const config: PipelineConfig = bundle.downstreamResources
    ? restorePipelineConfig(bundle.config, bundle.downstreamResources, bundle.contaminationReferences)
    : { dataset: bundle.config.dataset, contaminationPanel: bundle.config.contaminationPanel,
      contaminationPanelSequences: (bundle.contaminationReferences ?? []).map((record) => ({ ...record })),
      parameters: { ...bundle.config.parameters, deterministicSeed: BigInt(bundle.config.parameters.deterministicSeed) },
      samples: bundle.config.samples.map((sample) => ({ ...sample, panelSequences: [] })) };
  let current: ResultBundle = { ...bundle, optionalStages: explicitStatuses(bundle) };

  if (targetIndex >= 0 && !stageCompleted(current, "contamination")) {
    const started = performance.now(); report(options, "contamination", 0, "Preparing run-wide contamination checks");
    const contamination = await classifyContaminationAsync(current.consensuses, config, options.signal, (state) =>
      report(options, "contamination", state.fraction, state.detail));
    const summaries = current.summaries.map((summary) => ({ ...summary,
      contaminationPassed: summary.consensusSequences - contamination.filter((call) => call.sample === summary.sample && call.discarded).length }));
    const stamp = new Date().toISOString(), statuses = explicitStatuses(current);
    statuses.contamination = statusRecord("completed", `${contamination.filter((call) => call.discarded).length} consensus sequences excluded.`, stamp);
    current = checkpoint(options, { ...current, contamination, summaries, optionalStages: statuses,
      timings: timing(current.timings, "contamination", (performance.now() - started) / 1000, current.consensuses.length),
      log: [...current.log, `${stamp} on-demand contamination: completed; ${contamination.filter((call) => call.discarded).length} discarded`] });
  }
  if (targetIndex >= 1 && !stageCompleted(current, "postprocessing")) {
    const started = performance.now(); report(options, "postprocessing", 0, "Starting panel screening, alignments, functional checks, and annotations");
    const downstream = await postprocess(current.consensuses, current.contamination, config, options.signal, runAlivibeMsa,
      Math.max(1, current.provenance.workers), (state) => report(options, "postprocessing", state.fraction, state.detail), { collapse: false });
    downstream.summaries.forEach((summary, index) => {
      const base = current.summaries[index]; summary.demultiplexedReads = base.demultiplexedReads;
      summary.selectedReads = base.selectedReads; summary.downsampledReads = base.downsampledReads;
      summary.observedUmis = base.observedUmis; summary.likelyRealUmis = base.likelyRealUmis;
    });
    const stamp = new Date().toISOString(), statuses = explicitStatuses(current);
    statuses.postprocessing = statusRecord("completed", `${downstream.records.length} consensus-family records evaluated.`, stamp);
    current = checkpoint(options, { ...current, records: downstream.records, summaries: downstream.summaries,
      alignments: downstream.alignments, referenceAlignments: downstream.referenceAlignments, collapseGroups: {},
      trees: {}, optionalStages: statuses,
      timings: timing(current.timings, "postprocessing", (performance.now() - started) / 1000, downstream.records.length),
      log: [...current.log, `${stamp} on-demand postprocessing: completed for ${downstream.records.length} consensus families`] });
  }
  if (targetIndex >= 2 && !stageCompleted(current, "collapse")) {
    const started = performance.now(), input: PostprocessOutput = { records: current.records, summaries: current.summaries,
      alignments: current.alignments, referenceAlignments: current.referenceAlignments ?? {}, collapseGroups: current.collapseGroups ?? {}, collapseSeconds: 0 };
    const collapsed = await collapsePostprocess(input, config, options.signal,
      (state) => report(options, "collapse", state.fraction, state.detail));
    const count = Object.values(collapsed.collapseGroups).reduce((sum, groups) => sum + groups.length, 0), stamp = new Date().toISOString();
    const statuses = explicitStatuses(current); statuses.collapse = statusRecord("completed", `${count} haplotypes; multiplicities count UMI families.`, stamp);
    current = checkpoint(options, { ...current, summaries: collapsed.summaries, alignments: collapsed.alignments,
      referenceAlignments: collapsed.referenceAlignments, collapseGroups: collapsed.collapseGroups, optionalStages: statuses,
      timings: timing(current.timings, "collapse", (performance.now() - started) / 1000, count),
      log: [...current.log, `${stamp} on-demand collapse: ${count} haplotypes; multiplicities count UMI families`] });
  }
  if (targetIndex >= 3 && !stageCompleted(current, "tree")) {
    const inputs = config.samples.flatMap((sample) => {
      const key = `${sample.name}/nucleotide`, fasta = current.alignments[key]; return fasta ? [[key, fasta] as const] : [];
    });
    const trees = { ...current.trees }, started = performance.now(); let cursor = 0, finished = 0;
    await Promise.all(Array.from({ length: Math.min(Math.max(1, current.provenance.workers), Math.max(1, inputs.length)) }, async () => {
      while (true) {
        if (options.signal?.aborted) throw new DOMException("Phylogeny inference skipped.", "AbortError");
        const index = cursor++; if (index >= inputs.length) return;
        const [key, fasta] = inputs[index]; report(options, "tree", finished / Math.max(1, inputs.length),
          `Inferring ${key.split("/")[0]} (${finished} of ${inputs.length} complete) with bundled double-precision FastTree`);
        try { trees[key] = await runFastTreeIsolated(fasta, options.signal); }
        catch (cause) { if (options.signal?.aborted) throw cause; trees[key] = starTree(fasta); }
        finished++; report(options, "tree", finished / Math.max(1, inputs.length), `Finished ${finished} of ${inputs.length} collapsed phylogenies`);
      }
    }));
    const stamp = new Date().toISOString(), statuses = explicitStatuses(current);
    statuses.tree = statusRecord("completed", inputs.length ? `${inputs.length} collapsed phylogenies inferred.` : "No retained collapsed alignments required a tree.", stamp);
    current = checkpoint(options, { ...current, trees, optionalStages: statuses,
      timings: timing(current.timings, "tree", (performance.now() - started) / 1000, inputs.length),
      log: [...current.log, `${stamp} on-demand phylogeny: ${inputs.length} collapsed trees inferred with bundled double-precision FastTree`] });
  }
  return current;
}
