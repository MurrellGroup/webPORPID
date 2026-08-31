export type PanelFilterMode = "mafft-batch" | "independent-query";

export interface PipelineParameters {
  errorRate: number;
  minLength: number;
  maxLength: number;
  primerTolerance: number;
  primerWindow: number;
  primerChop: number;
  maxReadsPerSample: number;
  familySizeThreshold: number;
  ldaThreshold: number;
  contaminationClusterThreshold: number;
  contaminationProportionThreshold: number;
  contaminationDistanceThreshold: number;
  contaminationFilter: boolean;
  agreementThreshold: number;
  artefactFraction: number;
  outlierQuantile: number;
  panelThreshold: number;
  /** Reference-panel alignment strategy; absent in projects through 0.3.8. */
  panelFilterMode: PanelFilterMode;
  functionalMatchThreshold: number;
  spoolPartitions: number;
  deterministicSeed: bigint;
}

export interface NamedSequence { name: string; sequence: string }

export interface SampleConfig {
  name: string;
  /** Optional biological donor grouping. Samples sharing a donor are self for contamination checks. */
  donorId?: string;
  cdnaPrimer: string;
  secondStrandPrimer: string;
  panel: string;
  functionalReference?: string;
  panelSequences: NamedSequence[];
  functionalReferenceSequence?: NamedSequence;
  familySizeOverride?: number;
  artefactFractionOverride?: number;
  outlierQuantileOverride?: number;
  agreementOverride?: number;
  functionalMatchOverride?: number;
}

export interface PipelineConfig {
  dataset: string;
  samples: SampleConfig[];
  contaminationPanel: string;
  contaminationPanelSequences: NamedSequence[];
  parameters: PipelineParameters;
}

export interface QualityStats {
  totalReads: number;
  qualityReads: number;
  badReads: number;
  shortReads: number;
  longReads: number;
  primerRejects: number;
  idRejects: number;
  demultiplexedReads: number;
  bpbRejects: number;
  malformedRecords: number;
  downsampledReads: number;
  perSample: number[];
}

export type FamilyDisposition =
  | "likely_real" | "BPB-rejects" | "heteroduplex" | "LDA-rejects"
  | "UMI_len != 8" | "family-size-reject";

export interface UmiFamily {
  sample: string;
  sampleIndex: number;
  umi: string;
  familySize: number;
  mostLikelyParent: string;
  posteriorProbability: number;
  logOffspringProbability: number;
  disposition: FamilyDisposition;
  minimumAgreement?: number;
}

export interface LowAgreementSite {
  position: number;
  agreement: number;
  modalReadBase: string;
  modalRunLength: number;
}

export interface ConsensusRecord {
  id: string;
  sample: string;
  sampleIndex: number;
  umi: string;
  familySize: number;
  minimumAgreement: number;
  sequence: string;
  lowAgreementSites: LowAgreementSite[];
}

export interface ContaminationCall {
  sample: string;
  sequenceId: string;
  nearestNonselfVariant: string;
  nearestNonselfDistance: number;
  flagged: boolean;
  discarded: boolean;
  suspectOnly: boolean;
}

export interface ApobecResult {
  posteriorMeanGaMultiplier: number;
  posteriorGaInflated: number;
  posteriorMeanMutationRate: number;
  gaMutations: number;
  totalMutations: number;
}

export interface PostprocRecord {
  id: string;
  sample: string;
  umi: string;
  familySize: number;
  minimumAgreement: number;
  consensusNt: string;
  alignedNt?: string;
  trimmedNt?: string;
  trimmedAa?: string;
  panelScore: number;
  artefactPass: boolean;
  agreementPass: boolean;
  contaminationPass: boolean;
  panelPass: boolean;
  functionalPass?: boolean;
  rejectionReasons: string[];
  apobec?: ApobecResult;
}

export interface SampleSummary {
  sample: string;
  demultiplexedReads: number;
  /** Reads retained for this sample after deterministic subsampling. */
  selectedReads?: number;
  /** Demultiplexed reads omitted by deterministic subsampling. */
  downsampledReads?: number;
  observedUmis: number;
  likelyRealUmis: number;
  consensusSequences: number;
  /** Absent when contamination checks have not been computed. */
  contaminationPassed?: number;
  /** Absent when downstream filtering has not been computed. */
  postprocPassed?: number;
  /** Number of distinct retained nucleotide haplotypes after family-level collapse. */
  collapsedSequences?: number;
  functionalPassed?: number;
  /** Absent when downstream filtering has not been computed. */
  artefactCutoff?: number;
}

export interface Provenance {
  webporpidVersion: string;
  createdUtc: string;
  engine: string;
  workers: number;
  inputName: string;
  inputSha256: string;
  configSha256: string;
  deterministicSeed: string;
  upstreamBranch: "nanopore";
  upstreamCommit: string;
}

export interface ResultConfig {
  dataset: string;
  samples: Array<Omit<SampleConfig, "panelSequences" | "functionalReferenceSequence">>;
  contaminationPanel: string;
  parameters: Omit<PipelineParameters, "deterministicSeed"> & { deterministicSeed: string };
}

export interface PipelineTiming {
  stage: "setup" | "preprocessing" | "umi" | "consensus" | "contamination" | "postprocessing" | "collapse" | "tree" | "analysis-total";
  seconds: number;
  workItems?: number;
}

export interface AlignmentEdit {
  fasta: string;
  frameOffset: 0 | 1 | 2;
  baselineFingerprint: string;
  editedFingerprint: string;
  source: string;
  savedUtc: string;
  treeNewick?: string;
  /** Fingerprint of the alignment used for treeNewick. */
  treeFingerprint?: string;
  /** True when the alignment changed after the most recent tree inference. */
  treeStale?: boolean;
  warnings?: string[];
  /** Current edited alignment compared with the immutable pipeline alignment. */
  changes?: AlignmentChangeSummary;
}

export interface AlignmentChangeSummary {
  rowsBefore: number;
  rowsAfter: number;
  columnsBefore: number;
  columnsAfter: number;
  rowOrderChanged: boolean;
  rowOrderBefore: string[];
  rowOrderAfter: string[];
  removedRows: string[];
  addedRows: string[];
  changedRows: string[];
  rowChanges: AlignmentRowChange[];
  removedNucleotides: number;
  insertedNucleotides: number;
  substitutedNucleotides: number;
}

export interface AlignmentRowChange {
  name: string;
  substitutedNucleotides: number;
  insertedNucleotides: number;
  removedNucleotides: number;
  gapPlacementChanged: boolean;
}

export interface AlignmentAuditEntry {
  alignmentKey: string;
  action: "alignment-edit" | "frame-change" | "tree-recalculation" | "edit-reset";
  timestamp: string;
  source: string;
  details: string[];
  beforeFingerprint?: string;
  afterFingerprint?: string;
}

export interface CollapseGroup {
  sample: string;
  representativeId: string;
  memberIds: string[];
  /** One count per retained UMI family; deliberately not a read count. */
  familyCount: number;
  /**
   * Legacy (<=0.3.2) conservative summary: the lowest minimumAgreement among
   * member UMI families. New results omit it because agreement is a family
   * property, not a collapsed-haplotype property.
   */
  minimumAgreement?: number;
}

export interface InputFileMapping {
  slot: string;
  role: "reads" | "configuration" | "panel" | "functional-reference" | "contamination-panel";
  expectedName?: string;
  uploadedName: string;
  uploadedSize: number;
}

export type OptionalStageName = "contamination" | "postprocessing" | "collapse" | "tree";
export type OptionalStageState = "completed" | "deferred" | "skipped";

export interface OptionalStageStatus {
  state: OptionalStageState;
  detail: string;
  updatedUtc: string;
}

/**
 * Reference data required to resume optional work from a consensus-only
 * project. These inputs are small compared with the FASTQ and are therefore
 * retained in the project without retaining raw reads.
 */
export interface DownstreamResources {
  samples: Array<{
    name: string;
    panelSequences: NamedSequence[];
    functionalReferenceSequence?: NamedSequence;
  }>;
}

export interface RunOptions {
  deferPhylogeny: boolean;
  deferContamination?: boolean;
  deferPostprocessing?: boolean;
  deferCollapse?: boolean;
  /** Browser read-spool location; absent in older results and CLI runs. */
  spoolStorage?: "automatic" | "external-directory";
  /** Pause at inspectable UMI and consensus-filter decision checkpoints. */
  interactiveFiltering?: boolean;
}

export type ThresholdReviewPhase = "umi" | "consensus-filters";

export interface ThresholdReviewSample {
  sample: string;
  donorId?: string;
  totalFamilies: number;
  /** Exact family-size frequency table; pairs are [family size, family count]. */
  familySizeCounts: Array<[number, number]>;
  /** Exact 201-bin posterior distribution on [0,1], including both endpoints and excluding the aggregate BPB bucket. */
  posteriorBins?: number[];
  /** Exact 101-bin minimum-agreement distribution on [0,1]. */
  agreementBins?: number[];
  /** A deterministic display-only sample. Thresholds are always applied to every family. */
  displayPoints?: Array<{ familySize: number; posteriorProbability?: number; minimumAgreement?: number; disposition: FamilyDisposition }>;
  current: {
    familySizeThreshold?: number;
    artefactFraction?: number;
    outlierQuantile?: number;
    agreementThreshold?: number;
  };
  usesGlobal: {
    familySizeThreshold?: boolean;
    artefactFraction?: boolean;
    outlierQuantile?: boolean;
    agreementThreshold?: boolean;
  };
}

export interface ThresholdReview {
  id: string;
  phase: ThresholdReviewPhase;
  title: string;
  detail: string;
  current: {
    ldaThreshold?: number;
    familySizeThreshold?: number;
    artefactFraction?: number;
    outlierQuantile?: number;
    agreementThreshold?: number;
  };
  samples: ThresholdReviewSample[];
}

export interface ThresholdSelection {
  id: string;
  phase: ThresholdReviewPhase;
  parameters: Partial<Pick<PipelineParameters, "ldaThreshold" | "familySizeThreshold" | "artefactFraction" | "outlierQuantile" | "agreementThreshold">>;
  samples: Array<{
    sample: string;
    familySizeOverride?: number;
    artefactFractionOverride?: number;
    outlierQuantileOverride?: number;
    agreementOverride?: number;
  }>;
}

export interface ThresholdSelectionRecord extends ThresholdSelection {
  acceptedUtc: string;
  changes: string[];
}

export interface ResultBundle {
  schema: "webporpid-results/1";
  provenance: Provenance;
  config: ResultConfig;
  quality: QualityStats;
  summaries: SampleSummary[];
  umiFamilies: UmiFamily[];
  consensuses: ConsensusRecord[];
  contamination: ContaminationCall[];
  /** External contamination-panel records retained for on-demand reference/donor phylogenies. */
  contaminationReferences?: NamedSequence[];
  /** Small panel/reference inputs needed to resume deferred downstream stages. */
  downstreamResources?: DownstreamResources;
  records: PostprocRecord[];
  alignments: Record<string, string>;
  trees: Record<string, string>;
  /** Reference rows projected into each stored nucleotide alignment. */
  referenceAlignments?: Record<string, string>;
  /** Collapsed haplotypes and the UMI-family members represented by each tip. */
  collapseGroups?: Record<string, CollapseGroup[]>;
  /** Exact user-file to configuration-slot assignments used for this run. */
  inputMappings?: InputFileMapping[];
  runOptions?: RunOptions;
  /** Explicit state for every optional stage; absent in results through 0.3.5. */
  optionalStages?: Record<OptionalStageName, OptionalStageStatus>;
  /** Whether stored post-processing decisions applied computed contamination calls or deliberately retained every consensus at that gate. */
  postprocessingContaminationMode?: "applied" | "bypassed";
  /** Manual Alivibe/import corrections, keyed like `sample/nucleotide`. */
  alignmentEdits?: Record<string, AlignmentEdit>;
  /** Append-only record of interactive alignment and tree operations. */
  alignmentEditHistory?: AlignmentAuditEntry[];
  /** Optional so results written by webPORPID 0.1.x remain loadable. */
  timings?: PipelineTiming[];
  /** Accepted interactive decision checkpoints, including exact parameter changes. */
  thresholdSelections?: ThresholdSelectionRecord[];
  log: string[];
}

export interface PipelineProgress {
  stage: "preprocessing" | "umi" | "consensus" | "contamination" | "postprocessing" | "collapse" | "tree" | "complete";
  fraction: number;
  detail: string;
  reads?: number;
  /** Live demultiplexed-read counts in configuration sample order. */
  sampleAssignments?: Array<{ sample: string; reads: number }>;
}
