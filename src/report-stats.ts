import type { FamilyDisposition, PostprocRecord, ResultBundle, UmiFamily } from "./types";

export interface DualCountStat {
  key: string;
  label: string;
  families: number;
  familyPercent: number;
  reads: number;
  readPercent: number;
  note?: string;
}

export interface CountStat {
  key: string;
  label: string;
  count: number;
  percent: number;
  note?: string;
}

export interface SampleOverviewStat {
  sample: string;
  demultiplexedReads: number;
  selectedReads: number;
  downsampledReads: number;
  downsampledPercent: number;
  observedFamilies: number;
  consensusFamilies: number;
  retainedFamilies: number;
  functionalConfigured: boolean;
  functionalEvaluatedFamilies: number;
  functionalPassedFamilies: number;
  collapsedHaplotypes: number;
  bpbReadPercent: number;
  umiLengthReadPercent: number;
  familySizeReadPercent: number;
  ldaReadPercent: number;
  heteroduplexReadPercent: number;
  artefactReadPercent: number;
  agreementReadPercent: number;
  contaminationReadPercent: number;
  panelReadPercent: number;
  functionalReadPercent: number;
}

export interface ParameterSettingRow {
  scope: "run" | "sample";
  sample: string;
  parameter: string;
  value: string;
}

export const FAMILY_CALLS: ReadonlyArray<{ disposition: FamilyDisposition; label: string; note: string }> = [
  { disposition: "likely_real", label: "Accepted UMI families", note: "Final UMI-family call after all PORPID family-level checks." },
  { disposition: "BPB-rejects", label: "BPB rejects", note: "Reads for which the barcode/primer boundary could not be accepted." },
  { disposition: "UMI_len != 8", label: "UMI-length rejects", note: "Observed UMI length did not match the configured eight-base model." },
  { disposition: "family-size-reject", label: "Family-size rejects", note: "UMI family did not satisfy the family-size requirement." },
  { disposition: "LDA-rejects", label: "LDA rejects", note: "Probabilistic offspring classifier rejected the UMI family." },
  { disposition: "heteroduplex", label: "Heteroduplex rejects", note: "Consensus evidence was consistent with a heteroduplex family." },
];

export const FUNCTIONAL_REASONS = [
  { key: "ambiguous", label: "Ambiguous-symbol rejects", matches: (reason: string) => reason.startsWith("ambiguousSymbols-reject") },
  { key: "frameshift", label: "Frameshift rejects", matches: (reason: string) => reason.startsWith("frameshift-reject") },
  { key: "late-start", label: "Late-start rejects", matches: (reason: string) => reason.startsWith("lateStart-reject") },
  { key: "early-stop", label: "Early-stop rejects", matches: (reason: string) => reason.startsWith("earlyStop-reject") },
  { key: "bad-match", label: "Functional-reference match rejects", matches: (reason: string) => reason.startsWith("badMatch-reject") },
] as const;

const percent = (value: number, denominator: number) => denominator ? value / denominator * 100 : 0;
const readCount = <T extends { familySize: number }>(rows: readonly T[]) => rows.reduce((sum, row) => sum + row.familySize, 0);

function dualRow<T extends { familySize: number }>(key: string, label: string, selected: readonly T[], all: readonly T[], note?: string): DualCountStat {
  const reads = readCount(selected), denominatorReads = readCount(all);
  return { key, label, families: selected.length, familyPercent: percent(selected.length, all.length), reads,
    readPercent: percent(reads, denominatorReads), note };
}

export function selectedReadCount(bundle: ResultBundle, sample: string): number {
  const summary = bundle.summaries.find((row) => row.sample === sample);
  if (summary?.selectedReads != null) return summary.selectedReads;
  const represented = readCount(bundle.umiFamilies.filter((row) => row.sample === sample));
  return represented || summary?.demultiplexedReads || 0;
}

export function porpidCallStats(bundle: ResultBundle, sample: string): DualCountStat[] {
  const entries = bundle.umiFamilies.filter((row) => row.sample === sample), families = entries.filter((row) => row.disposition !== "BPB-rejects");
  const denominatorReads = readCount(entries);
  return [
    { key: "all", label: "All observed UMI families", families: families.length, familyPercent: families.length ? 100 : 0,
      reads: denominatorReads, readPercent: denominatorReads ? 100 : 0,
      note: "UMI percentages exclude the aggregate BPB bucket; read / CCS percentages include every selected read." },
    ...FAMILY_CALLS.map(({ disposition, label, note }) => {
      const selected = entries.filter((row) => row.disposition === disposition), reads = readCount(selected);
      const count = disposition === "BPB-rejects" ? 0 : selected.length;
      return { key: disposition, label, families: count, familyPercent: percent(count, families.length), reads,
        readPercent: percent(reads, denominatorReads), note: disposition === "BPB-rejects" ? `${note} These reads have no accepted UMI, so no artificial UMI-family count is reported.` : note };
    }),
  ];
}

const passedPostproc = (record: PostprocRecord) => record.artefactPass && record.agreementPass && record.contaminationPass && record.panelPass;

export function postprocFilterStats(bundle: ResultBundle, sample: string): DualCountStat[] {
  const records = bundle.records.filter((row) => row.sample === sample);
  const matching = (predicate: (row: PostprocRecord) => boolean) => records.filter(predicate);
  const contaminationBypassed = bundle.postprocessingContaminationMode === "bypassed";
  return [
    dualRow("all", "All consensus families", records, records, "Denominator for both percentage columns; filter rows can overlap."),
    dualRow("artefact", "Artefact-cutoff rejects", matching((row) => !row.artefactPass), records),
    dualRow("agreement", "Minimum-agreement rejects", matching((row) => !row.agreementPass), records),
    dualRow("contamination", contaminationBypassed ? "Contamination filter bypassed" : "Contamination rejects",
      matching((row) => !row.contaminationPass), records, contaminationBypassed
        ? "This optional check was skipped, deferred, or disabled when downstream outputs were built; zero sequences were discarded at this gate."
        : undefined),
    dualRow("panel", "Reference-panel rejects", matching((row) => !row.panelPass), records),
    dualRow("retained", "Passed non-functional filters", matching(passedPostproc), records),
  ];
}

export function functionalWasEvaluated(record: PostprocRecord): boolean {
  if (record.functionalPass === true || record.trimmedNt != null || record.trimmedAa != null) return true;
  return record.rejectionReasons.some((reason) => FUNCTIONAL_REASONS.some((category) => category.matches(reason)));
}

export function functionalFilterStats(bundle: ResultBundle, sample: string): DualCountStat[] {
  const records = bundle.records.filter((row) => row.sample === sample);
  const configured = Boolean(bundle.config.samples.find((row) => row.name === sample)?.functionalReference);
  if (!configured) return [dualRow("not-configured", "No functional filter configured", records, records,
    "Nucleotide and amino-acid alignment views remain available without a functional reference.")];
  const evaluated = records.filter(functionalWasEvaluated), notEvaluated = records.filter((row) => !functionalWasEvaluated(row));
  return [
    dualRow("all", "All consensus families", records, records, "Denominator for both percentage columns; reason rows can overlap."),
    dualRow("not-evaluated", "Not evaluated after upstream rejection", notEvaluated, records),
    dualRow("evaluated", "Evaluated against functional reference", evaluated, records),
    dualRow("passed", "Functional-filter passes", evaluated.filter((row) => row.functionalPass === true), records),
    dualRow("rejected", "Any functional-filter reject", evaluated.filter((row) => row.functionalPass === false), records),
    ...FUNCTIONAL_REASONS.map((category) => dualRow(category.key, category.label,
      evaluated.filter((row) => row.rejectionReasons.some(category.matches)), records)),
  ];
}

function dispositionReadPercent(families: readonly UmiFamily[], disposition: FamilyDisposition): number {
  return percent(readCount(families.filter((row) => row.disposition === disposition)), readCount(families));
}

function rejectionReadPercent(records: readonly PostprocRecord[], predicate: (row: PostprocRecord) => boolean): number {
  return percent(readCount(records.filter(predicate)), readCount(records));
}

export function sampleOverviewStats(bundle: ResultBundle): SampleOverviewStat[] {
  return bundle.config.samples.map((configured) => {
    const sample = configured.name, summary = bundle.summaries.find((row) => row.sample === sample);
    const families = bundle.umiFamilies.filter((row) => row.sample === sample), records = bundle.records.filter((row) => row.sample === sample);
    const consensuses = bundle.consensuses.filter((row) => row.sample === sample);
    const contaminantIds = new Set(bundle.contamination.filter((row) => row.sample === sample && row.discarded).map((row) => row.sequenceId));
    const demultiplexedReads = summary?.demultiplexedReads ?? 0, selectedReads = selectedReadCount(bundle, sample);
    const downsampledReads = summary?.downsampledReads ?? Math.max(0, demultiplexedReads - selectedReads);
    const evaluated = records.filter(functionalWasEvaluated), functionalConfigured = Boolean(configured.functionalReference);
    return {
      sample, demultiplexedReads, selectedReads, downsampledReads,
      downsampledPercent: percent(downsampledReads, demultiplexedReads),
      observedFamilies: families.filter((row) => row.disposition !== "BPB-rejects").length,
      consensusFamilies: summary?.consensusSequences ?? consensuses.length,
      retainedFamilies: records.filter(passedPostproc).length,
      functionalConfigured,
      functionalEvaluatedFamilies: evaluated.length,
      functionalPassedFamilies: evaluated.filter((row) => row.functionalPass === true).length,
      collapsedHaplotypes: summary?.collapsedSequences ?? bundle.collapseGroups?.[sample]?.length ?? 0,
      bpbReadPercent: dispositionReadPercent(families, "BPB-rejects"),
      umiLengthReadPercent: dispositionReadPercent(families, "UMI_len != 8"),
      familySizeReadPercent: dispositionReadPercent(families, "family-size-reject"),
      ldaReadPercent: dispositionReadPercent(families, "LDA-rejects"),
      heteroduplexReadPercent: dispositionReadPercent(families, "heteroduplex"),
      artefactReadPercent: rejectionReadPercent(records, (row) => !row.artefactPass),
      agreementReadPercent: rejectionReadPercent(records, (row) => !row.agreementPass),
      contaminationReadPercent: records.length && bundle.postprocessingContaminationMode !== "bypassed" ? rejectionReadPercent(records, (row) => !row.contaminationPass)
        : percent(readCount(consensuses.filter((row) => contaminantIds.has(row.id))), readCount(consensuses)),
      panelReadPercent: rejectionReadPercent(records, (row) => !row.panelPass),
      functionalReadPercent: percent(readCount(evaluated.filter((row) => row.functionalPass === false)), readCount(records)),
    };
  });
}

export function inputFilterStats(bundle: ResultBundle): CountStat[] {
  const quality = bundle.quality, total = quality.totalReads + quality.malformedRecords;
  const row = (key: string, label: string, count: number, note?: string): CountStat => ({ key, label, count, percent: percent(count, total), note });
  return [
    row("all", "FASTQ records encountered", total, "Parsed records plus malformed records; denominator for this run-wide table."),
    row("parsed", "Parsed FASTQ records", quality.totalReads),
    row("malformed", "Malformed FASTQ records", quality.malformedRecords),
    row("quality-passed", "Passed expected-error filter", quality.qualityReads),
    row("bad-quality", "Quality-filter rejects", quality.badReads),
    row("short", "Too-short rejects", quality.shortReads),
    row("long", "Too-long rejects", quality.longReads),
    row("primer", "Primer-filter rejects", quality.primerRejects),
    row("identifier", "Read-identifier rejects", quality.idRejects),
    row("demultiplexed", "Assigned to a configured sample", quality.demultiplexedReads),
    row("unassigned", "Not assigned after preprocessing", Math.max(0, total - quality.demultiplexedReads), "Aggregate outcome; it overlaps malformed and specific upstream rejection rows above."),
    row("subsampled", "Demultiplexed reads omitted by subsampling", quality.downsampledReads),
  ];
}

export function parameterSettings(bundle: ResultBundle): ParameterSettingRow[] {
  const rows: ParameterSettingRow[] = [];
  const addRun = (parameter: string, value: unknown) => rows.push({ scope: "run", sample: "—", parameter, value: String(value ?? "—") });
  addRun("dataset", bundle.config.dataset); addRun("contaminationPanel", bundle.config.contaminationPanel);
  for (const [parameter, value] of Object.entries(bundle.config.parameters)) addRun(parameter, value);
  addRun("workers", bundle.provenance.workers); addRun("engine", bundle.provenance.engine);
  addRun("temporaryReadStorage", bundle.runOptions?.spoolStorage ?? "not recorded");
  addRun("postprocessingContaminationMode", bundle.postprocessingContaminationMode ?? "not recorded (legacy result)");
  for (const [parameter, value] of Object.entries(bundle.runOptions ?? {})) if (parameter !== "spoolStorage") addRun(parameter, value);
  for (const sample of bundle.config.samples) {
    const add = (parameter: string, value: unknown) => rows.push({ scope: "sample", sample: sample.name, parameter, value: String(value ?? "—") });
    add("cDNA primer", sample.cdnaPrimer); add("second-strand primer", sample.secondStrandPrimer);
    add("panel", sample.panel); add("functional reference", sample.functionalReference ?? "not configured");
    add("familySizeThreshold (effective)", sample.familySizeOverride ?? bundle.config.parameters.familySizeThreshold);
    add("artefactFraction (effective)", sample.artefactFractionOverride ?? bundle.config.parameters.artefactFraction);
    add("outlierQuantile (effective)", sample.outlierQuantileOverride ?? bundle.config.parameters.outlierQuantile);
    add("agreementThreshold (effective)", sample.agreementOverride ?? bundle.config.parameters.agreementThreshold);
    add("functionalMatchThreshold (effective)", sample.functionalMatchOverride ?? bundle.config.parameters.functionalMatchThreshold);
  }
  return rows;
}
