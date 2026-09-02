import type {
  ConsensusRecord, PipelineConfig, ThresholdReview, ThresholdSelection, ThresholdSelectionRecord, UmiFamily,
} from "./types";

const DISPLAY_POINT_LIMIT = 4_000;

function frequencyTable(values: readonly number[]): Array<[number, number]> {
  const counts = new Map<number, number>();
  for (const value of values) counts.set(value, (counts.get(value) ?? 0) + 1);
  return [...counts].sort((left, right) => left[0] - right[0]);
}

function stableDisplaySample<T>(rows: readonly T[], limit = DISPLAY_POINT_LIMIT): T[] {
  if (rows.length <= limit) return [...rows];
  const output: T[] = [];
  for (let index = 0; index < limit; index++) output.push(rows[Math.floor(index * rows.length / limit)]);
  return output;
}

function bins(values: readonly number[], divisions: number): number[] {
  const output = Array(divisions + 1).fill(0) as number[];
  for (const value of values) output[Math.max(0, Math.min(divisions, Math.round(value * divisions)))]++;
  return output;
}

function reviewId(phase: ThresholdReview["phase"]) {
  return `${phase}-${Date.now()}-${Math.random().toString(36).slice(2)}`;
}

export function buildUmiThresholdReview(families: readonly UmiFamily[], config: PipelineConfig): ThresholdReview {
  const rowsBySample = Array.from({ length: config.samples.length }, () => [] as UmiFamily[]);
  for (const family of families) if (family.disposition !== "BPB-rejects" && rowsBySample[family.sampleIndex])
    rowsBySample[family.sampleIndex].push(family);
  return {
    id: reviewId("umi"), phase: "umi", title: "Review UMI-family thresholds",
    detail: "The offspring model is complete. Adjust the probability and family-size decisions before consensus calling. Display points may be thinned, but accepted thresholds are applied to every observed family.",
    current: { ldaThreshold: config.parameters.ldaThreshold, familySizeThreshold: config.parameters.familySizeThreshold },
    samples: config.samples.map((sample, sampleIndex) => {
      const rows = rowsBySample[sampleIndex];
      return { sample: sample.name, donorId: sample.donorId, totalFamilies: rows.length,
        familySizeCounts: frequencyTable(rows.map((row) => row.familySize)), posteriorBins: bins(rows.map((row) => row.posteriorProbability), 200),
        displayPoints: stableDisplaySample(rows).map((row) => ({ umi: row.umi, familySize: row.familySize,
          posteriorProbability: row.posteriorProbability, disposition: row.disposition })),
        current: { familySizeThreshold: sample.familySizeOverride ?? config.parameters.familySizeThreshold },
        usesGlobal: { familySizeThreshold: sample.familySizeOverride === undefined } };
    }),
  };
}

export function buildConsensusThresholdReview(consensuses: readonly ConsensusRecord[], discardedIds: ReadonlySet<string>,
  config: PipelineConfig, families?: readonly UmiFamily[]): ThresholdReview {
  const rowsBySample = Array.from({ length: config.samples.length }, () => [] as ConsensusRecord[]);
  const abundanceRowsBySample = Array.from({ length: config.samples.length }, () => [] as ConsensusRecord[]);
  for (const row of consensuses) if (rowsBySample[row.sampleIndex]) {
    rowsBySample[row.sampleIndex].push(row);
    if (!discardedIds.has(row.id)) abundanceRowsBySample[row.sampleIndex].push(row);
  }
  return {
    id: reviewId("consensus-filters"), phase: "consensus-filters", title: "Review consensus-family filters",
    detail: "Consensus agreement and family abundance are complete. Adjust the agreement, outlier-quantile, and artefact-fraction thresholds before reference-panel screening and retained-family alignment.",
    current: { artefactFraction: config.parameters.artefactFraction, outlierQuantile: config.parameters.outlierQuantile,
      agreementThreshold: config.parameters.agreementThreshold },
    samples: config.samples.map((sample, sampleIndex) => {
      const allRows = rowsBySample[sampleIndex], abundanceRows = abundanceRowsBySample[sampleIndex];
      const agreementByUmi = new Map(allRows.map((row) => [row.umi, row.minimumAgreement]));
      const decisionFamilies = families?.filter((row) => row.sampleIndex === sampleIndex
        && (row.disposition === "likely_real" || row.disposition === "family-size-reject"));
      const displayPoints = families
        ? stableDisplaySample(decisionFamilies ?? []).map((row) => ({
            umi: row.umi, familySize: row.familySize, minimumAgreement: agreementByUmi.get(row.umi), disposition: row.disposition,
          }))
        : stableDisplaySample(allRows).map((row) => ({ umi: row.umi, familySize: row.familySize,
          minimumAgreement: row.minimumAgreement, disposition: "likely_real" as const }));
      return { sample: sample.name, donorId: sample.donorId, totalFamilies: decisionFamilies?.length ?? allRows.length,
        familySizeCounts: frequencyTable(abundanceRows.map((row) => row.familySize)),
        agreementBins: bins(allRows.map((row) => row.minimumAgreement), 100),
        displayPoints,
        current: { artefactFraction: sample.artefactFractionOverride ?? config.parameters.artefactFraction,
          outlierQuantile: sample.outlierQuantileOverride ?? config.parameters.outlierQuantile,
          agreementThreshold: sample.agreementOverride ?? config.parameters.agreementThreshold },
        usesGlobal: { artefactFraction: sample.artefactFractionOverride === undefined,
          outlierQuantile: sample.outlierQuantileOverride === undefined,
          agreementThreshold: sample.agreementOverride === undefined } };
    }),
  };
}

function finite(value: number | undefined, label: string): number | undefined {
  if (value === undefined) return undefined;
  if (!Number.isFinite(value)) throw new Error(`${label} must be a finite number.`);
  return value;
}

function familyThreshold(value: number | undefined, label: string): number | undefined {
  value = finite(value, label);
  if (value !== undefined && (!Number.isSafeInteger(value) || value < 0))
    throw new Error(`${label} must be a non-negative integer. It may exceed the slider range.`);
  return value;
}

function setChange(changes: string[], label: string, before: number | undefined, after: number | undefined) {
  if (before !== after) changes.push(`${label}: ${before ?? "global"} -> ${after ?? "global"}`);
}

/** Validate, apply, and describe one interactive decision without re-fitting the offspring model. */
export function applyThresholdSelection(config: PipelineConfig, families: UmiFamily[], selection: ThresholdSelection): ThresholdSelectionRecord {
  const changes: string[] = [], bySample = new Map(selection.samples.map((row) => [row.sample, row]));
  if (bySample.size !== selection.samples.length) throw new Error("An interactive threshold selection contains duplicate sample rows.");
  for (const name of bySample.keys()) if (!config.samples.some((sample) => sample.name === name))
    throw new Error(`Interactive thresholds reference unknown sample ${name}.`);

  if (selection.phase === "umi") {
    const lda = finite(selection.parameters.ldaThreshold, "Offspring probability threshold") ?? config.parameters.ldaThreshold;
    const globalFamily = familyThreshold(selection.parameters.familySizeThreshold, "Global family-size threshold") ?? config.parameters.familySizeThreshold;
    setChange(changes, "ldaThreshold", config.parameters.ldaThreshold, lda);
    setChange(changes, "familySizeThreshold", config.parameters.familySizeThreshold, globalFamily);
    config.parameters.ldaThreshold = lda; config.parameters.familySizeThreshold = globalFamily;
    config.samples.forEach((sample) => {
      const selected = bySample.get(sample.name);
      if (!selected) return;
      const threshold = familyThreshold(selected.familySizeOverride, `${sample.name} family-size threshold`);
      setChange(changes, `${sample.name}.familySizeOverride`, sample.familySizeOverride, threshold);
      sample.familySizeOverride = threshold;
    });
    for (const family of families) {
      if (family.disposition === "BPB-rejects") continue;
      const threshold = config.samples[family.sampleIndex]?.familySizeOverride ?? globalFamily;
      family.disposition = family.posteriorProbability < lda ? "LDA-rejects"
        : family.umi.length !== 8 ? "UMI_len != 8"
          : family.familySize < threshold ? "family-size-reject" : "likely_real";
    }
  } else {
    const artefact = finite(selection.parameters.artefactFraction, "Artefact fraction") ?? config.parameters.artefactFraction;
    const quantile = finite(selection.parameters.outlierQuantile, "Outlier quantile") ?? config.parameters.outlierQuantile;
    const agreement = finite(selection.parameters.agreementThreshold, "Minimum agreement") ?? config.parameters.agreementThreshold;
    setChange(changes, "artefactFraction", config.parameters.artefactFraction, artefact);
    setChange(changes, "outlierQuantile", config.parameters.outlierQuantile, quantile);
    setChange(changes, "agreementThreshold", config.parameters.agreementThreshold, agreement);
    config.parameters.artefactFraction = artefact; config.parameters.outlierQuantile = quantile; config.parameters.agreementThreshold = agreement;
    config.samples.forEach((sample) => {
      const selected = bySample.get(sample.name); if (!selected) return;
      const nextArtefact = finite(selected.artefactFractionOverride, `${sample.name} artefact fraction`);
      const nextQuantile = finite(selected.outlierQuantileOverride, `${sample.name} outlier quantile`);
      const nextAgreement = finite(selected.agreementOverride, `${sample.name} minimum agreement`);
      setChange(changes, `${sample.name}.artefactFractionOverride`, sample.artefactFractionOverride, nextArtefact);
      setChange(changes, `${sample.name}.outlierQuantileOverride`, sample.outlierQuantileOverride, nextQuantile);
      setChange(changes, `${sample.name}.agreementOverride`, sample.agreementOverride, nextAgreement);
      sample.artefactFractionOverride = nextArtefact; sample.outlierQuantileOverride = nextQuantile; sample.agreementOverride = nextAgreement;
    });
  }
  return { ...selection, acceptedUtc: new Date().toISOString(), changes: changes.length ? changes : ["Accepted the configured thresholds without changes."] };
}
