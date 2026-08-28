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
  functionalMatchThreshold: number;
  spoolPartitions: number;
  deterministicSeed: bigint;
}

export interface NamedSequence { name: string; sequence: string }

export interface SampleConfig {
  name: string;
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
  observedUmis: number;
  likelyRealUmis: number;
  consensusSequences: number;
  contaminationPassed: number;
  postprocPassed: number;
  /** Number of distinct retained nucleotide haplotypes after family-level collapse. */
  collapsedSequences?: number;
  functionalPassed?: number;
  artefactCutoff: number;
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

export interface RunOptions {
  deferPhylogeny: boolean;
  /** Browser read-spool location; absent in older results and CLI runs. */
  spoolStorage?: "automatic" | "external-directory";
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
  /** Manual Alivibe/import corrections, keyed like `sample/nucleotide`. */
  alignmentEdits?: Record<string, AlignmentEdit>;
  /** Optional so results written by webPORPID 0.1.x remain loadable. */
  timings?: PipelineTiming[];
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
