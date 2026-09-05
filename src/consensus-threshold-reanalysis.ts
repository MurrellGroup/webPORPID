import { OPTIONAL_STAGE_ORDER, restorePipelineConfig, stageCompleted, stageStatus, statusRecord } from "./optional-stages.ts";
import { applyThresholdSelection, buildConsensusThresholdReview } from "./threshold-review.ts";
import type { OptionalStageName, PipelineConfig, ResultBundle, ResultConfig, ThresholdReview, ThresholdSelection } from "./types.ts";

function explicitStatuses(bundle: ResultBundle): NonNullable<ResultBundle["optionalStages"]> {
  return Object.fromEntries(OPTIONAL_STAGE_ORDER.map((stage) => [stage, { ...stageStatus(bundle, stage) }])) as NonNullable<ResultBundle["optionalStages"]>;
}

export function restoreStoredPipelineConfig(bundle: ResultBundle): PipelineConfig {
  return bundle.downstreamResources
    ? restorePipelineConfig(bundle.config, bundle.downstreamResources, bundle.contaminationReferences)
    : { dataset: bundle.config.dataset, contaminationPanel: bundle.config.contaminationPanel,
      contaminationPanelSequences: (bundle.contaminationReferences ?? []).map((record) => ({ ...record })),
      parameters: { ...bundle.config.parameters, deterministicSeed: BigInt(bundle.config.parameters.deterministicSeed) },
      samples: bundle.config.samples.map((sample) => ({ ...sample, panelSequences: [] })) };
}

// Kept local so this replay planner stays directly testable without loading
// browser-only pipeline modules. It mirrors config.resultConfig exactly.
function compactResultConfig(config: PipelineConfig): ResultConfig {
  return {
    dataset: config.dataset,
    samples: config.samples.map(({ panelSequences: _panel, functionalReferenceSequence: _reference, ...sample }) => sample),
    contaminationPanel: config.contaminationPanel,
    parameters: { ...config.parameters, deterministicSeed: config.parameters.deterministicSeed.toString() },
  };
}

/** Reconstruct the post-consensus decision dialog entirely from the compact project. */
export function buildStoredConsensusThresholdReview(bundle: ResultBundle): ThresholdReview {
  const contaminationApplied = stageCompleted(bundle, "contamination") && bundle.config.parameters.contaminationFilter;
  const discardedIds = contaminationApplied
    ? new Set(bundle.contamination.filter((call) => call.discarded).map((call) => call.sequenceId))
    : new Set<string>();
  const review = buildConsensusThresholdReview(bundle.consensuses, discardedIds, restoreStoredPipelineConfig(bundle), bundle.umiFamilies);
  return { ...review, title: "Review consensus-family filters again",
    detail: "Stored consensus agreement and family abundance are loaded. Adjust agreement, outlier-quantile, and artefact-fraction thresholds. Accepting them will rerun only the previously completed downstream stages." };
}

export interface ConsensusThresholdReanalysisPlan {
  bundle: ResultBundle;
  /** Highest previously completed downstream stage; absent when nothing downstream has run yet. */
  target?: Exclude<OptionalStageName, "contamination">;
  /** Samples whose effective downstream thresholds changed. */
  affectedSamples: string[];
  /** Global changes require a run-wide replay; otherwise affectedSamples is exact. */
  scope: "none" | "samples" | "all";
  /** Statuses to restore for downstream stages that had deliberately remained deferred or skipped. */
  untouchedStatuses: Partial<NonNullable<ResultBundle["optionalStages"]>>;
}

/**
 * Apply an accepted consensus-filter decision and invalidate only products
 * derived after that checkpoint. UMI grouping, consensus calls, and computed
 * contamination calls remain byte-for-byte untouched.
 */
export function prepareConsensusThresholdReanalysis(bundle: ResultBundle, selection: ThresholdSelection): ConsensusThresholdReanalysisPlan {
  if (selection.phase !== "consensus-filters") throw new Error("Only consensus-filter thresholds can be replayed from a completed project.");
  const config = restoreStoredPipelineConfig(bundle), before = restoreStoredPipelineConfig(bundle);
  const families = bundle.umiFamilies.map((family) => ({ ...family }));
  const accepted = applyThresholdSelection(config, families, selection), stamp = accepted.acceptedUtc;
  const postprocessingDone = stageCompleted(bundle, "postprocessing"), collapseDone = stageCompleted(bundle, "collapse"), treeDone = stageCompleted(bundle, "tree");
  const priorTarget: ConsensusThresholdReanalysisPlan["target"] = treeDone ? "tree" : collapseDone ? "collapse" : postprocessingDone ? "postprocessing" : undefined;
  const globalChanged = before.parameters.artefactFraction !== config.parameters.artefactFraction
    || before.parameters.outlierQuantile !== config.parameters.outlierQuantile
    || before.parameters.agreementThreshold !== config.parameters.agreementThreshold;
  const effective = (candidate: PipelineConfig, index: number) => {
    const sample = candidate.samples[index];
    return [sample.artefactFractionOverride ?? candidate.parameters.artefactFraction,
      sample.outlierQuantileOverride ?? candidate.parameters.outlierQuantile,
      sample.agreementOverride ?? candidate.parameters.agreementThreshold] as const;
  };
  let affectedSamples = config.samples.flatMap((sample, index) => {
    const left = effective(before, index), right = effective(config, index);
    return left.some((value, field) => value !== right[field]) ? [sample.name] : [];
  });
  // Reapplying previously bypassed contamination to only a subset would make
  // the run-wide contamination-mode flag scientifically ambiguous.
  const contaminationReapply = priorTarget && bundle.postprocessingContaminationMode === "bypassed"
    && stageCompleted(bundle, "contamination") && config.parameters.contaminationFilter;
  const scope: ConsensusThresholdReanalysisPlan["scope"] = globalChanged || contaminationReapply ? "all"
    : affectedSamples.length ? "samples" : "none";
  if (scope === "all") affectedSamples = config.samples.map((sample) => sample.name);
  const target = scope === "none" ? undefined : priorTarget;
  const statuses = explicitStatuses(bundle), untouchedStatuses: ConsensusThresholdReanalysisPlan["untouchedStatuses"] = {};
  if (!collapseDone) untouchedStatuses.collapse = { ...statuses.collapse };
  if (!treeDone) untouchedStatuses.tree = { ...statuses.tree };
  if (target && !bundle.downstreamResources)
    throw new Error("This project does not contain the panel/reference sequences required to rerun downstream analysis. Re-run from the FASTQ with webPORPID 0.3.6 or later.");

  const detail = "Invalidated after an accepted run-wide consensus-threshold change; recomputation is starting from stored consensus calls.";
  if (target && scope === "all") {
    statuses.postprocessing = statusRecord("deferred", detail, stamp);
    statuses.collapse = collapseDone ? statusRecord("deferred", detail, stamp) : statuses.collapse;
    statuses.tree = treeDone ? statusRecord("deferred", detail, stamp) : statuses.tree;
  }
  const changed = accepted.changes.join("; ");
  const replayDetail = scope === "all" ? target ? `recomputing all samples through ${target}` : "global thresholds saved; no completed downstream stage required recomputation"
    : scope === "samples" ? target ? `recomputing only changed samples through ${target}: ${affectedSamples.join(", ")}` : `sample thresholds saved for ${affectedSamples.join(", ")}; no completed downstream stage required recomputation`
      : "no effective threshold changed; no downstream recomputation required";
  return { target, affectedSamples, scope, untouchedStatuses, bundle: {
    ...bundle,
    config: compactResultConfig(config),
    umiFamilies: families,
    thresholdSelections: [...(bundle.thresholdSelections ?? []), accepted],
    optionalStages: statuses,
    ...(target && scope === "all" ? { records: [], alignments: {}, referenceAlignments: {}, collapseGroups: {}, trees: {}, alignmentEdits: undefined,
      functionalFilterErrors: undefined,
      postprocessingContaminationMode: undefined,
      timings: bundle.timings?.filter((entry) => !["postprocessing", "collapse", "tree"].includes(entry.stage)) } : {}),
    log: [...bundle.log,
      `${stamp} interactive consensus filter reopened after completion: ${changed}`,
      `${stamp} consensus-threshold replay: UMI grouping, consensus calls, and contamination calls preserved; ${replayDetail}${target && scope === "all" ? "; all derived alignments, trees, and active edits detached" : target && scope === "samples" ? "; derived products and active edits will be replaced only for those samples" : ""}`],
  } };
}

export function restoreUntouchedThresholdStatuses(bundle: ResultBundle, plan: ConsensusThresholdReanalysisPlan): ResultBundle {
  const statuses = explicitStatuses(bundle);
  for (const [stage, status] of Object.entries(plan.untouchedStatuses))
    if (status) statuses[stage as OptionalStageName] = { ...status };
  return { ...bundle, optionalStages: statuses };
}
