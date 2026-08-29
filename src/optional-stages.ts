import type {
  DownstreamResources, OptionalStageName, OptionalStageStatus, PipelineConfig, ResultBundle, ResultConfig,
} from "./types.ts";

export const OPTIONAL_STAGE_ORDER: readonly OptionalStageName[] = ["contamination", "postprocessing", "collapse", "tree"];

export function stageCompleted(bundle: ResultBundle, stage: OptionalStageName): boolean {
  const explicit = bundle.optionalStages?.[stage];
  if (explicit) return explicit.state === "completed";
  // Results through 0.3.5 always ran the first three stages. Tree inference
  // alone could be deferred and was recorded in runOptions.
  return stage === "tree" ? !bundle.runOptions?.deferPhylogeny : true;
}

export function stageStatus(bundle: ResultBundle, stage: OptionalStageName): OptionalStageStatus {
  const explicit = bundle.optionalStages?.[stage];
  if (explicit) return explicit;
  const completed = stageCompleted(bundle, stage);
  return {
    state: completed ? "completed" : "deferred",
    detail: completed ? "Computed by the original analysis." : "Deferred by the original analysis.",
    updatedUtc: bundle.provenance.createdUtc,
  };
}

export function downstreamResources(config: PipelineConfig): DownstreamResources {
  return { samples: config.samples.map((sample) => ({
    name: sample.name,
    panelSequences: sample.panelSequences.map((record) => ({ ...record })),
    functionalReferenceSequence: sample.functionalReferenceSequence ? { ...sample.functionalReferenceSequence } : undefined,
  })) };
}

export function restorePipelineConfig(config: ResultConfig, resources: DownstreamResources | undefined,
  contaminationReferences: ResultBundle["contaminationReferences"]): PipelineConfig {
  if (!resources) throw new Error("This project predates resumable downstream analysis and does not contain its panel/reference sequences. Re-run from the FASTQ to compute deferred stages.");
  const bySample = new Map(resources.samples.map((sample) => [sample.name, sample]));
  return {
    dataset: config.dataset,
    contaminationPanel: config.contaminationPanel,
    contaminationPanelSequences: (contaminationReferences ?? []).map((record) => ({ ...record })),
    parameters: { ...config.parameters, deterministicSeed: BigInt(config.parameters.deterministicSeed) },
    samples: config.samples.map((sample) => {
      const stored = bySample.get(sample.name);
      if (!stored) throw new Error(`The project is missing downstream reference data for sample ${sample.name}.`);
      return {
        ...sample,
        panelSequences: stored.panelSequences.map((record) => ({ ...record })),
        functionalReferenceSequence: stored.functionalReferenceSequence ? { ...stored.functionalReferenceSequence } : undefined,
      };
    }),
  };
}

export function statusRecord(state: OptionalStageStatus["state"], detail: string, updatedUtc = new Date().toISOString()): OptionalStageStatus {
  return { state, detail, updatedUtc };
}
