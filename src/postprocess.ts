import { runAlivibeMsa } from "./alivibe-msa-runtime.ts";
import { runIndependentPanelFilter } from "./independent-panel-filter-runtime.ts";
import type { PanelFilterResult } from "./independent-panel-filter.ts";
import { runMafftFftnsMsa } from "./mafft-msa-runtime.ts";
import { inspectAlignment, translateAlignedNucleotides } from "./alignment-utils.ts";
import { collapseAlignment } from "./collapse.ts";
import { extractAndScorePanel } from "./panel-profile.ts";
import { runScalableMsa } from "./scalable-msa.ts";
import { FUNCTIONAL_REFERENCE_NAME, functionalSequenceName, uncollapsedSequenceName } from "./sequence-names.ts";
import { addSequenceToProfile } from "./independent-panel-filter.ts";
import type { ApobecResult, ConsensusRecord, ContaminationCall, PipelineConfig, PostprocRecord, SampleSummary } from "./types.ts";

export interface PostprocessOutput {
  records: PostprocRecord[];
  summaries: SampleSummary[];
  alignments: Record<string, string>;
  referenceAlignments: Record<string, string>;
  collapseGroups: Record<string, ReturnType<typeof collapseAlignment>["groups"]>;
  /** A functional-filter failure is isolated to its sample instead of aborting the run. */
  functionalFilterErrors: Record<string, string>;
  collapseSeconds: number;
}

export type MsaRunner = (sequences: readonly string[], signal?: AbortSignal, iterations?: number,
  scoringMode?: "literal" | "nucleotide" | "amino-acid", onProgress?: (progress: { detail: string }) => void) => Promise<string[]>;

export type PanelFilterRunner = (sequences: readonly string[], panelRows: readonly string[], signal?: AbortSignal,
  workers?: number, onProgress?: (progress: { completed: number; total: number }) => void) => Promise<PanelFilterResult>;

export interface PostprocessProgress { fraction: number; detail: string }
export interface PostprocessOptions {
  /** Keep false when collapse is run as its own cancellable pipeline stage. */
  collapse?: boolean;
  onCollapseProgress?: (progress: PostprocessProgress) => void;
  /** Batch panel MSA. The browser default is the bundled MAFFT FFT-NS-2 runtime. */
  panelMsa?: MsaRunner;
  /** Independent query-to-panel implementation, injected by the CLI. */
  panelFilter?: PanelFilterRunner;
  /** CPU workers available to each active sample's independent panel filter. */
  panelWorkers?: number;
  /** Restrict recomputation to these samples while preserving original sample indices. */
  sampleNames?: ReadonlySet<string>;
}

const degap = (sequence: string) => sequence.replaceAll("-", "").toUpperCase();
const fasta = (rows: Array<{ name: string; sequence: string }>) => rows.map((row) => `>${row.name}\n${row.sequence.match(/.{1,80}/g)?.join("\n") ?? ""}`).join("\n") + (rows.length ? "\n" : "");

function quantile(values: number[], probability: number) {
  if (!values.length) return 0; const sorted = [...values].sort((a, b) => a - b); if (sorted.length === 1) return sorted[0];
  const position = Math.max(0, Math.min(1, probability)) * (sorted.length - 1), lower = Math.floor(position), fraction = position - lower;
  return sorted[lower] * (1 - fraction) + sorted[Math.min(lower + 1, sorted.length - 1)] * fraction;
}

function alignmentConsensus(rows: string[]) {
  if (!rows.length) return ""; let output = "";
  const counts = new Uint32Array(128), touched = new Uint8Array(128);
  for (let position = 0; position < rows[0].length; position++) {
    let touchedCount = 0;
    for (const row of rows) {
      let code = row.charCodeAt(position); if (code >= 97 && code <= 122) code -= 32;
      if (code >= counts.length) code = 63;
      if (counts[code] === 0) touched[touchedCount++] = code;
      counts[code]++;
    }
    let bestCode = touched[0] ?? 45, bestCount = counts[bestCode];
    // Touched symbols are in first-observed order; strict greater therefore
    // preserves the original modal tie rule after final counts are known.
    for (let index = 1; index < touchedCount; index++) if (counts[touched[index]] > bestCount) {
      bestCode = touched[index]; bestCount = counts[bestCode];
    }
    output += String.fromCharCode(bestCode);
    for (let index = 0; index < touchedCount; index++) counts[touched[index]] = 0;
  }
  return output;
}

function matrixMultiply(left: number[], right: number[]) {
  const output = Array(16).fill(0);
  for (let i = 0; i < 4; i++) for (let k = 0; k < 4; k++) for (let j = 0; j < 4; j++) output[i * 4 + j] += left[i * 4 + k] * right[k * 4 + j];
  return output;
}
const identity = () => [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1];
function matrixExp(matrix: number[]) {
  const norm = Math.max(...[0, 1, 2, 3].map((row) => matrix.slice(row * 4, row * 4 + 4).reduce((sum, value) => sum + Math.abs(value), 0)));
  const scale = Math.max(0, Math.ceil(Math.log2(Math.max(norm, 1)))), divisor = 2 ** scale;
  matrix = matrix.map((value) => value / divisor); let result = identity(), term = identity();
  for (let order = 1; order <= 24; order++) { term = matrixMultiply(term, matrix).map((value) => value / order); result = result.map((value, index) => value + term[index]); }
  for (let index = 0; index < scale; index++) result = matrixMultiply(result, result); return result;
}
function normalPdf(value: number, mean: number, sd: number) { return Math.exp(-((value - mean) ** 2) / (2 * sd * sd)) / (sd * Math.sqrt(2 * Math.PI)); }
interface GridPoint { t: number; ga: number; logTransitions: number[]; prior: number }
let gridCache: GridPoint[] | undefined;
function apobecGrid() {
  if (gridCache) return gridCache; const output: GridPoint[] = [];
  for (let ti = 0; ti < 23; ti++) for (let gi = 0; gi < 121; gi++) {
    const t = -12 + ti * 0.5, ga = -1 + gi * 0.05, mu = Math.exp(t), multiplier = Math.exp(ga), tv = 4.5;
    const q = [-(1 + tv + 1), 1, tv, 1, 1, -(1 + 1 + tv), 1, tv,
      tv * multiplier, 1, -(tv * multiplier + 1 + 1), 1, 1, tv, 1, -(1 + tv + 1)].map((value) => value * mu);
    const prior = Math.log(normalPdf(t, -5, 1)) + Math.log(0.99 * normalPdf(ga, 0, 0.1) + 0.01 * normalPdf(ga, 0, 1));
    const transitions = matrixExp(q);
    output.push({ t, ga, logTransitions: transitions.map((value) => Math.log(Math.max(value, 1e-300))), prior });
  }
  return (gridCache = output);
}

const baseIndex: Record<string, number> = { A: 0, C: 1, G: 2, T: 3 };
function apobec(consensus: string, sequence: string): ApobecResult {
  const counts = Array(16).fill(0);
  for (let position = 0; position < Math.min(consensus.length, sequence.length); position++) {
    const left = baseIndex[consensus[position].toUpperCase()], right = baseIndex[sequence[position].toUpperCase()];
    if (left !== undefined && right !== undefined) counts[left * 4 + right]++;
  }
  const observed: number[] = []; for (let index = 0; index < 16; index++) if (counts[index]) observed.push(index);
  const grid = apobecGrid(), logWeight = (point: GridPoint) => {
    let value = point.prior; for (const index of observed) value += counts[index] * point.logTransitions[index]; return value;
  };
  const weights = grid.map(logWeight), maximum = Math.max(...weights);
  const scaled = weights.map((value) => Math.exp(value - maximum)), total = scaled.reduce((sum, value) => sum + value, 0);
  let t = 0, ga = 0, inflated = 0;
  grid.forEach((point, index) => { const probability = scaled[index] / total; t += point.t * probability; ga += point.ga * probability; if (point.ga > 0) inflated += probability; });
  let totalMutations = 0; for (let a = 0; a < 4; a++) for (let b = 0; b < 4; b++) if (a !== b) totalMutations += counts[a * 4 + b];
  return { posteriorMeanGaMultiplier: Math.exp(ga), posteriorGaInflated: inflated, posteriorMeanMutationRate: Math.exp(t),
    gaMutations: counts[2 * 4], totalMutations };
}

const CODONS: Record<string, string> = {
  TTT: "F", TTC: "F", TTA: "L", TTG: "L", CTT: "L", CTC: "L", CTA: "L", CTG: "L",
  ATT: "I", ATC: "I", ATA: "I", ATG: "M", GTT: "V", GTC: "V", GTA: "V", GTG: "V",
  TCT: "S", TCC: "S", TCA: "S", TCG: "S", AGT: "S", AGC: "S", CCT: "P", CCC: "P", CCA: "P", CCG: "P",
  ACT: "T", ACC: "T", ACA: "T", ACG: "T", GCT: "A", GCC: "A", GCA: "A", GCG: "A",
  TAT: "Y", TAC: "Y", TAA: "*", TAG: "*", TGA: "*", CAT: "H", CAC: "H", CAA: "Q", CAG: "Q",
  AAT: "N", AAC: "N", AAA: "K", AAG: "K", GAT: "D", GAC: "D", GAA: "E", GAG: "E",
  TGT: "C", TGC: "C", TGG: "W", CGT: "R", CGC: "R", CGA: "R", CGG: "R", AGA: "R", AGG: "R",
  GGT: "G", GGC: "G", GGA: "G", GGG: "G",
};
function translate(sequence: string) { let output = ""; for (let index = 0; index + 2 < sequence.length; index += 3) output += CODONS[sequence.slice(index, index + 3)] ?? "X"; return output; }
function longestOrf(sequence: string) {
  let best: { length: number; start: number; end: number } | undefined;
  for (let frame = 0; frame < 3; frame++) {
    const aa = translate(sequence.slice(frame, frame + Math.floor((sequence.length - frame) / 3) * 3));
    for (let start = 0; start < aa.length; start++) if (aa[start] === "M") {
      const offset = aa.indexOf("*", start); if (offset >= 0 && (!best || offset + 1 - start > best.length)) best = { length: offset + 1 - start, start: frame + start * 3, end: frame + (offset + 1) * 3 };
    }
  }
  return best ? sequence.slice(best.start, best.end) : undefined;
}

function pairwise(reference: string, query: string) {
  const columns = query.length + 1, trace = new Uint8Array((reference.length + 1) * columns);
  let previous = new Int32Array(columns), current = new Int32Array(columns);
  for (let column = 1; column < columns; column++) { previous[column] = previous[column - 1] - 99; trace[column] = 2; }
  for (let row = 1; row <= reference.length; row++) {
    current[0] = previous[0] - 99; trace[row * columns] = 3;
    for (let column = 1; column < columns; column++) {
      const diagonal = previous[column - 1] + (reference[row - 1] === query[column - 1] ? 100 : -100);
      const left = current[column - 1] - (row === reference.length ? 99 : 100), up = previous[column] - (column === query.length ? 99 : 100);
      if (diagonal >= left && diagonal >= up) { current[column] = diagonal; trace[row * columns + column] = 1; }
      else if (left >= up) { current[column] = left; trace[row * columns + column] = 2; } else { current[column] = up; trace[row * columns + column] = 3; }
    }
    [previous, current] = [current, previous];
  }
  let left = "", right = "", row = reference.length, column = query.length;
  while (row || column) { const op = trace[row * columns + column]; if (op === 1) { left += reference[--row]; right += query[--column]; } else if (op === 2) { left += "-"; right += query[--column]; } else { left += reference[--row]; right += "-"; } }
  return { reference: [...left].reverse().join(""), query: [...right].reverse().join("") };
}

export function functionalFilter(reference: string, sequence: string, threshold: number) {
  sequence = degap(sequence); const reasons: string[] = [];
  if (/[^ACGT]/.test(sequence)) return { passed: false, reasons: ["ambiguousSymbols-reject"] };
  const coding = longestOrf(sequence); if (!coding) return { passed: false, reasons: ["noORF-reject"] };
  if (coding.length % 3) return { passed: false, reasons: ["frameshift-reject"] };
  const aligned = pairwise(degap(reference), coding), first = aligned.reference.search(/[^-]/);
  let last = aligned.reference.length - 1; while (last >= first && aligned.reference[last] === "-") last--;
  const referenceRegion = aligned.reference.slice(first, last + 1), queryRegion = aligned.query.slice(first, last + 1);
  const frameshift = [referenceRegion, queryRegion].some((row) => {
    for (const match of row.matchAll(/-+/g)) {
      const start = match.index, end = start + match[0].length;
      if (start > 0 && end < row.length && match[0].length % 3 !== 0) return true;
    }
    return false;
  });
  if (frameshift) reasons.push("frameshift-reject");
  const trimmed = degap(queryRegion), aa = translate(trimmed);
  if (queryRegion.slice(0, 3) !== "ATG") reasons.push("lateStart-reject");
  if (queryRegion.slice(-3) === "---") reasons.push("earlyStop-reject");
  let matches = 0; for (let index = 0; index < referenceRegion.length; index++) if (referenceRegion[index] === queryRegion[index]) matches++;
  const rawRatio = matches / Math.max(1, referenceRegion.length);
  const digits = rawRatio === 0 ? 3 : 3 - Math.floor(Math.log10(Math.abs(rawRatio))) - 1;
  const ratio = Number(rawRatio.toFixed(Math.max(0, digits)));
  if (ratio < threshold) reasons.push(`badMatch-reject (match=${ratio})`);
  return { passed: !reasons.length, reasons, nt: trimmed, aa, referenceMatch: Number(rawRatio.toFixed(2)) };
}

interface FunctionalOutcome {
  passed: boolean;
  reasons: string[];
  nt?: string;
  aa?: string;
  /** Exact-position identity to the clipped aligned reference, rounded to two decimal places. */
  referenceMatch?: number;
  /** Codon-aware alignment, clipped to the first/last reference residue. */
  alignedNt?: string;
  alignedAa?: string;
}

interface FunctionalBatchOutcome {
  outcomes: FunctionalOutcome[];
  referenceName: string;
  referenceNt: string;
  referenceAa: string;
}

function backtranslate(alignedAminoAcids: string, codingNucleotides: string): string {
  let output = "", offset = 0;
  for (const residue of alignedAminoAcids) {
    if (residue === "-") output += "---";
    else { output += codingNucleotides.slice(offset, offset + 3); offset += 3; }
  }
  if (offset !== codingNucleotides.length) throw new Error("The amino-acid alignment did not preserve its nucleotide coding sequence.");
  return output;
}

interface AnchorParts { insertions: string[]; residues: string[] }

function splitAgainstAnchor(anchor: string, alignedAnchor: string, alignedSequence: string): AnchorParts {
  if (alignedAnchor.length !== alignedSequence.length || degap(alignedAnchor) !== anchor)
    throw new Error("The sequence-preserving amino-acid fallback returned an invalid anchor.");
  const insertions = Array.from({ length: anchor.length + 1 }, () => ""), residues = Array<string>(anchor.length);
  let position = 0;
  for (let column = 0; column < alignedAnchor.length; column++) {
    if (alignedAnchor[column] === "-") insertions[position] += alignedSequence[column];
    else residues[position++] = alignedSequence[column];
  }
  if (position !== anchor.length || degap(alignedSequence).length === 0)
    throw new Error("The sequence-preserving amino-acid fallback returned an invalid row.");
  return { insertions, residues };
}

/**
 * Correctness fallback for a structurally corrupt POA result. Every query is
 * aligned independently to the first (abundance-leading) protein, then the
 * insertion slots are padded into one rectangular alignment. It is linear in
 * the number of rows and cannot delete or invent a residue.
 */
function anchorProteinAlignment(sequences: readonly string[]): string[] {
  if (sequences.length < 2) return [...sequences];
  const anchor = sequences[0], rows: AnchorParts[] = Array(sequences.length);
  rows[0] = { insertions: Array.from({ length: anchor.length + 1 }, () => ""), residues: [...anchor] };
  const widths = new Uint32Array(anchor.length + 1);
  for (let index = 1; index < sequences.length; index++) {
    const aligned = addSequenceToProfile([anchor], sequences[index]);
    const parts = splitAgainstAnchor(anchor, aligned.profileRows[0], aligned.sequence); rows[index] = parts;
    parts.insertions.forEach((value, slot) => { widths[slot] = Math.max(widths[slot], value.length); });
  }
  const result = rows.map((row) => {
    let output = row.insertions[0].padStart(widths[0], "-");
    for (let position = 0; position < anchor.length; position++)
      output += row.residues[position] + row.insertions[position + 1].padEnd(widths[position + 1], "-");
    return output;
  });
  result.forEach((row, index) => {
    if (degap(row) !== sequences[index]) throw new Error("The sequence-preserving amino-acid fallback did not preserve an input protein.");
  });
  return result;
}

/**
 * An MSA is allowed to choose gaps, never biological residues. Restore benign
 * alphabet/case normalization onto an unchanged gap path; use a strict
 * sequence-preserving fallback if the aligner changed a row's residue count.
 */
function preserveProteinSequences(inputs: readonly string[], aligned: readonly string[]): string[] {
  if (aligned.length !== inputs.length || aligned.some((row) => row.length !== aligned[0]?.length))
    return anchorProteinAlignment(inputs);
  let structuralChange = false;
  const restored = aligned.map((row, index) => {
    const residues = [...row].filter((residue) => residue !== "-");
    if (residues.length !== inputs[index].length) { structuralChange = true; return row; }
    let cursor = 0;
    return [...row].map((residue) => residue === "-" ? "-" : inputs[index][cursor++]).join("");
  });
  if (structuralChange) return anchorProteinAlignment(inputs);
  restored.forEach((row, index) => {
    if (degap(row) !== inputs[index]) throw new Error("The amino-acid alignment did not preserve an input protein.");
  });
  return restored;
}

/**
 * Align the eligible sample ORFs in amino-acid space first, freeze that sample
 * profile, and add the translated reference without moving existing residues.
 * Projecting every resulting gap back as a codon triplet prevents an aligner
 * from manufacturing one- or two-base frameshifts.
 */
async function functionalFilterBatch(
  referenceName: string, reference: string, sequences: string[], threshold: number, runMsa: MsaRunner, signal?: AbortSignal,
): Promise<FunctionalBatchOutcome> {
  const outcomes: FunctionalOutcome[] = Array(sequences.length), coding: Array<{ index: number; sequence: string }> = [];
  sequences.forEach((raw, index) => {
    const sequence = degap(raw);
    if (/[^ACGT]/.test(sequence)) outcomes[index] = { passed: false, reasons: ["ambiguousSymbols-reject"] };
    else {
      const orf = longestOrf(sequence);
      if (orf) coding.push({ index, sequence: orf });
      else outcomes[index] = { passed: false, reasons: ["noORF-reject"] };
    }
  });
  const referenceCoding = degap(reference).slice(0, Math.floor(degap(reference).length / 3) * 3);
  const resolvedReferenceName = referenceName.trim() || FUNCTIONAL_REFERENCE_NAME;
  if (!coding.length) return { outcomes, referenceName: resolvedReferenceName, referenceNt: referenceCoding, referenceAa: translate(referenceCoding) };

  // SeededAlignment first aligns the biological sequences to one another, then
  // adds the functional reference to that *fixed* profile. A joint MSA can move
  // gaps between sample sequences merely because a reference was introduced.
  const sampleAminoInputs = coding.map((row) => translate(row.sequence));
  const rawSampleAminoAlignment = sampleAminoInputs.length > 1
    ? await runScalableMsa(sampleAminoInputs, runMsa, signal, 3, "amino-acid")
    : sampleAminoInputs;
  const sampleAminoAlignment = preserveProteinSequences(sampleAminoInputs, rawSampleAminoAlignment);
  const withReference = addSequenceToProfile(sampleAminoAlignment, translate(referenceCoding));
  const aminoAlignment = [withReference.sequence, ...withReference.profileRows];
  const nucleotideAlignment = aminoAlignment.map((row, index) => backtranslate(row, index ? coding[index - 1].sequence : referenceCoding));
  const alignedReference = nucleotideAlignment[0], first = alignedReference.search(/[^-]/);
  let last = alignedReference.length - 1; while (last >= first && alignedReference[last] === "-") last--;
  const referenceRegion = alignedReference.slice(first, last + 1);
  coding.forEach((entry, position) => {
    const queryRegion = nucleotideAlignment[position + 1].slice(first, last + 1), reasons: string[] = [];
    const trimmed = degap(queryRegion), aa = translate(trimmed);
    if (queryRegion.slice(0, 3) !== "ATG") reasons.push("lateStart-reject");
    if (queryRegion.slice(-3) === "---") reasons.push("earlyStop-reject");
    let matches = 0; for (let column = 0; column < referenceRegion.length; column += 1) if (referenceRegion[column] === queryRegion[column]) matches++;
    const rawRatio = matches / Math.max(1, referenceRegion.length), digits = rawRatio === 0 ? 3 : 3 - Math.floor(Math.log10(Math.abs(rawRatio))) - 1;
    const ratio = Number(rawRatio.toFixed(Math.max(0, digits)));
    if (ratio < threshold) reasons.push(`badMatch-reject (match=${ratio})`);
    outcomes[entry.index] = { passed: !reasons.length, reasons, nt: trimmed, aa,
      alignedNt: queryRegion, alignedAa: translateAlignedNucleotides(queryRegion, 0), referenceMatch: Number(rawRatio.toFixed(2)) };
  });
  return { outcomes, referenceName: resolvedReferenceName, referenceNt: referenceRegion, referenceAa: translateAlignedNucleotides(referenceRegion, 0) };
}

export async function postprocess(
  consensuses: ConsensusRecord[], contamination: ContaminationCall[], config: PipelineConfig, signal?: AbortSignal,
  runMsa: MsaRunner = runAlivibeMsa, sampleConcurrency = 1, onProgress?: (progress: PostprocessProgress) => void,
  options: PostprocessOptions = {},
): Promise<PostprocessOutput> {
  const panelMsa = options.panelMsa ?? runMafftFftnsMsa;
  const independentPanelFilter = options.panelFilter ?? runIndependentPanelFilter;
  const discarded = new Set(contamination.filter((call) => call.discarded).map((call) => call.sequenceId));
  const activeSampleIndices = config.samples.flatMap((sample, index) =>
    !options.sampleNames || options.sampleNames.has(sample.name) ? [index] : []);
  const outputs: Array<PostprocessOutput | undefined> = Array(config.samples.length); let cursor = 0;
  const sourceBySample = Array.from({ length: config.samples.length }, () => [] as ConsensusRecord[]);
  const sampleIndexByName = new Map(config.samples.map((sample, index) => [sample.name, index]));
  for (const record of consensuses) {
    const sampleIndex = config.samples[record.sampleIndex]?.name === record.sample ? record.sampleIndex : sampleIndexByName.get(record.sample);
    if (sampleIndex !== undefined) sourceBySample[sampleIndex].push(record);
  }
  const sampleProgress = new Map(activeSampleIndices.map((index) => [index, 0]));
  const report = (sampleIndex: number, fraction: number, detail: string) => {
    sampleProgress.set(sampleIndex, Math.max(sampleProgress.get(sampleIndex) ?? 0, Math.max(0, Math.min(1, fraction))));
    onProgress?.({ fraction: [...sampleProgress.values()].reduce((sum, value) => sum + value, 0) / Math.max(1, sampleProgress.size), detail });
  };
  await Promise.all(Array.from({ length: Math.min(activeSampleIndices.length, Math.max(1, Math.floor(sampleConcurrency))) }, async () => {
    while (true) {
      const activeIndex = cursor++; if (activeIndex >= activeSampleIndices.length) return;
      const sampleIndex = activeSampleIndices[activeIndex];
      const sample = config.samples[sampleIndex], records: PostprocRecord[] = [], summaries: SampleSummary[] = [], alignments: Record<string, string> = {};
      const referenceAlignments: Record<string, string> = {}, collapseGroups: PostprocessOutput["collapseGroups"] = {};
    if (signal?.aborted) throw new DOMException("Analysis cancelled.", "AbortError");
    report(sampleIndex, .03, `Preparing filters for sample ${sample.name}`);
    const source = sourceBySample[sampleIndex], sizes: number[] = []; let discardedInSample = 0;
    for (const record of source) { if (discarded.has(record.id)) discardedInSample++; else sizes.push(record.familySize); }
    const artefactCutoff = Math.ceil(quantile(sizes, sample.outlierQuantileOverride ?? config.parameters.outlierQuantile)
      * (sample.artefactFractionOverride ?? config.parameters.artefactFraction));
    const agreementThreshold = sample.agreementOverride ?? config.parameters.agreementThreshold;
    const preliminary = source.map((record, index) => ({ record, index })).filter(({ record }) => record.familySize >= artefactCutoff
      && record.minimumAgreement >= agreementThreshold && !discarded.has(record.id));
    const scores = Array(source.length).fill(0), panelPass = Array(source.length).fill(true), extracted = new Map<number, string>();
    if (preliminary.length && sample.panelSequences.length) {
      const panelRows = sample.panelSequences.map((record) => record.sequence), rawCandidates = preliminary.map(({ record }) => degap(record.sequence));
      const panelMode = config.parameters.panelFilterMode ?? "mafft-batch";
      report(sampleIndex, .18, panelMode === "independent-query"
        ? `Aligning ${preliminary.length.toLocaleString()} candidate families independently to the reference panel for ${sample.name}`
        : `Building the batch MAFFT FFT-NS-2 alignment for ${preliminary.length.toLocaleString()} candidate families in ${sample.name}`);
      let panelResult: PanelFilterResult;
      if (panelMode === "independent-query") {
        const activeSamples = Math.max(1, Math.min(activeSampleIndices.length, Math.max(1, Math.floor(sampleConcurrency))));
        const panelWorkers = options.panelWorkers ?? Math.max(1, Math.floor(Math.max(1, sampleConcurrency) / activeSamples));
        panelResult = await independentPanelFilter(rawCandidates, panelRows, signal, panelWorkers, ({ completed, total }) =>
          report(sampleIndex, .18 + .2 * completed / Math.max(1, total),
            `Reference-panel alignment for ${sample.name}: ${completed.toLocaleString()} of ${total.toLocaleString()} families`));
      } else {
        const candidates = rawCandidates.length > 1
          ? await panelMsa(rawCandidates, signal, 0, "nucleotide", ({ detail }) =>
            report(sampleIndex, .24, `MAFFT FFT-NS-2 for ${sample.name}: ${detail}`))
          : rawCandidates;
        report(sampleIndex, .34, `Projecting the ${sample.name} candidate alignment onto its fixed reference panel`);
        panelResult = extractAndScorePanel(candidates, panelRows);
      }
      preliminary.forEach(({ index }, candidate) => {
        extracted.set(index, degap(panelResult.sequences[candidate])); scores[index] = panelResult.scores[candidate];
        panelPass[index] = scores[index] < config.parameters.panelThreshold;
      });
    } else {
      preliminary.forEach(({ record, index }) => extracted.set(index, degap(record.sequence)));
    }
    report(sampleIndex, .4, `Reference-panel screening complete for sample ${sample.name}`);
    const accepted = preliminary.filter(({ index }) => panelPass[index]);
    const displayReference = sample.functionalReferenceSequence?.sequence ?? sample.panelSequences[0]?.sequence;
    let acceptedAlignment: string[] = [], alignedReference = "";
    if (accepted.length) {
      report(sampleIndex, .48, `Building the retained-sequence alignment for sample ${sample.name}`);
      const inputs = [...(displayReference ? [degap(displayReference)] : []), ...accepted.map(({ index }) => extracted.get(index)!)];
      const aligned = inputs.length > 1 ? await runScalableMsa(inputs, runMsa, signal, 3, "nucleotide") : inputs;
      if (displayReference) { alignedReference = aligned[0]; acceptedAlignment = aligned.slice(1); }
      else { acceptedAlignment = aligned; alignedReference = alignmentConsensus(aligned); }
    }
    report(sampleIndex, .66, `Retained-sequence alignment complete for sample ${sample.name}`);
    const alignmentByIndex = new Map(accepted.map(({ index }, position) => [index, acceptedAlignment[position]]));
    const consensus = alignmentConsensus(acceptedAlignment), nucleotideRows: Array<{ name: string; sequence: string }> = [];
    report(sampleIndex, .72, `Calculating family-level annotations and filter decisions for sample ${sample.name}`);
    for (const [index, record] of source.entries()) {
      if (signal?.aborted) throw new DOMException("Downstream filtering skipped.", "AbortError");
      const artefactPass = record.familySize >= artefactCutoff, agreementPass = record.minimumAgreement >= agreementThreshold;
      const contaminationPass = !discarded.has(record.id), acceptedRow = alignmentByIndex.get(index), rejectionReasons: string[] = [];
      if (!artefactPass) rejectionReasons.push(`ccs_count < artefact cutoff (${artefactCutoff})`);
      if (!agreementPass) rejectionReasons.push(`minimum_agreement < ${agreementThreshold}`);
      if (!contaminationPass) rejectionReasons.push("contamination filter"); if (!panelPass[index]) rejectionReasons.push(`distance_from_panel >= ${config.parameters.panelThreshold}`);
      if (acceptedRow) nucleotideRows.push({ name: uncollapsedSequenceName(record), sequence: acceptedRow });
      records.push({ id: record.id, sample: sample.name, umi: record.umi, familySize: record.familySize,
        minimumAgreement: record.minimumAgreement, consensusNt: record.sequence, alignedNt: acceptedRow,
        panelScore: scores[index], artefactPass, agreementPass, contaminationPass, panelPass: panelPass[index],
        rejectionReasons, apobec: acceptedRow ? apobec(consensus, acceptedRow) : undefined });
      if ((index & 7) === 7 || index + 1 === source.length) {
        report(sampleIndex, .84 + .12 * (index + 1) / Math.max(1, source.length),
          `Annotated ${index + 1} of ${source.length} consensus-family records for ${sample.name}`);
        await new Promise<void>((resolvePromise) => setTimeout(resolvePromise, 0));
      }
    }
    if (nucleotideRows.length) {
      const uncollapsed = fasta(nucleotideRows);
      alignments[`${sample.name}/uncollapsed-nucleotide`] = uncollapsed;
      // Protein exploration is a direct view of every retained nucleotide row,
      // independent of the optional functional filter.  A complete gap codon
      // becomes '-', while mixed/ambiguous/incomplete codons become 'X'.
      alignments[`${sample.name}/uncollapsed-protein`] = fasta(nucleotideRows.map((row) => ({ ...row, sequence: translateAlignedNucleotides(row.sequence, 0) })));
      referenceAlignments[`${sample.name}/uncollapsed-nucleotide`] = fasta([{ name: "reference", sequence: alignedReference }]);
    }
    summaries.push({ sample: sample.name, demultiplexedReads: 0, observedUmis: 0, likelyRealUmis: 0,
      consensusSequences: source.length, contaminationPassed: source.length - discardedInSample,
      postprocPassed: accepted.length, artefactCutoff });
      outputs[sampleIndex] = { records, summaries, alignments, referenceAlignments, collapseGroups, functionalFilterErrors: {}, collapseSeconds: 0 };
      report(sampleIndex, 1, `Downstream processing complete for sample ${sample.name}`);
    }
  }));
  const completed = outputs.filter((output): output is PostprocessOutput => Boolean(output));
  const combined = { records: completed.flatMap((output) => output.records), summaries: completed.flatMap((output) => output.summaries),
    alignments: Object.assign({}, ...completed.map((output) => output.alignments)),
    referenceAlignments: Object.assign({}, ...completed.map((output) => output.referenceAlignments)),
    collapseGroups: Object.assign({}, ...completed.map((output) => output.collapseGroups)),
    functionalFilterErrors: Object.assign({}, ...completed.map((output) => output.functionalFilterErrors)), collapseSeconds: 0 };
  return options.collapse === false ? combined : collapsePostprocess(combined, config, signal, options.onCollapseProgress, runMsa);
}

/** Run family-count-preserving haplotype collapse as its own resumable stage. */
export async function collapsePostprocess(
  output: PostprocessOutput, config: PipelineConfig, signal?: AbortSignal,
  onProgress?: (progress: PostprocessProgress) => void,
  runMsa: MsaRunner = runAlivibeMsa,
  sampleNames?: ReadonlySet<string>,
): Promise<PostprocessOutput> {
  const started = performance.now(), alignments = { ...output.alignments }, referenceAlignments = { ...output.referenceAlignments };
  const collapseGroups = { ...output.collapseGroups }, summaries = output.summaries.map((summary) => ({ ...summary }));
  const functionalFilterErrors = { ...output.functionalFilterErrors };
  const activeSamples = config.samples.filter((sample) => !sampleNames || sampleNames.has(sample.name));
  for (const [index, sample] of activeSamples.entries()) {
    if (signal?.aborted) throw new DOMException("Haplotype collapse skipped.", "AbortError");
    onProgress?.({ fraction: index / Math.max(1, activeSamples.length), detail: `Collapsing identical retained UMI-family sequences for ${sample.name}` });
    await new Promise<void>((resolvePromise) => setTimeout(resolvePromise, 0));
    const uncollapsed = alignments[`${sample.name}/uncollapsed-nucleotide`];
    let collapsedCount = 0, functionalPassed: number | undefined;
    collapseGroups[sample.name] = [];
    delete functionalFilterErrors[sample.name];
    delete alignments[`${sample.name}/functional-nucleotide`]; delete alignments[`${sample.name}/functional-protein`];
    delete referenceAlignments[`${sample.name}/functional-nucleotide`];
    if (uncollapsed) {
      const collapsed = collapseAlignment(uncollapsed, sample.name); collapsedCount = collapsed.groups.length;
      alignments[`${sample.name}/nucleotide`] = collapsed.fasta;
      collapseGroups[sample.name] = collapsed.groups;
      const collapsedRows = inspectAlignment(collapsed.fasta, 1).records;
      alignments[`${sample.name}/protein`] = fasta(collapsedRows.map((row) => ({ name: row.name,
        sequence: translateAlignedNucleotides(row.sequence, 0) })));
      const reference = referenceAlignments[`${sample.name}/uncollapsed-nucleotide`];
      if (reference) referenceAlignments[`${sample.name}/nucleotide`] = reference;

      // Functional filtering belongs to biological variants, not individual
      // UMI families. It therefore runs once per collapsed haplotype and the
      // abundance-ranked collapsed identifiers are retained in every output.
      if (sample.functionalReferenceSequence && collapsedRows.length) {
        onProgress?.({ fraction: (index + .45) / Math.max(1, activeSamples.length),
          detail: `Checking ${collapsedRows.length.toLocaleString()} collapsed variants against the functional reference for ${sample.name}` });
        let batch: FunctionalBatchOutcome;
        try {
          batch = await functionalFilterBatch(sample.functionalReferenceSequence.name, sample.functionalReferenceSequence.sequence,
            collapsedRows.map((row) => row.sequence), sample.functionalMatchOverride ?? config.parameters.functionalMatchThreshold, runMsa, signal);
        } catch (cause) {
          if (signal?.aborted || (cause instanceof DOMException && cause.name === "AbortError")) throw cause;
          const message = cause instanceof Error ? cause.message : String(cause);
          functionalFilterErrors[sample.name] = message;
          onProgress?.({ fraction: (index + .8) / Math.max(1, activeSamples.length),
            detail: `Functional filtering could not be completed for ${sample.name}; recording the sample-specific error and continuing` });
          const summary = summaries.find((row) => row.sample === sample.name);
          if (summary) { summary.collapsedSequences = collapsedCount; delete summary.functionalPassed; }
          onProgress?.({ fraction: (index + 1) / Math.max(1, activeSamples.length),
            detail: `Collapsed ${sample.name} into ${collapsedCount.toLocaleString()} variants; functional filtering failed only for this sample` });
          continue;
        }
        const functionalNucleotideRows: Array<{ name: string; sequence: string }> = [
          { name: batch.referenceName, sequence: batch.referenceNt },
        ];
        const functionalProteinRows: Array<{ name: string; sequence: string }> = [
          { name: batch.referenceName, sequence: batch.referenceAa },
        ];
        let passedCount = 0;
        collapsed.groups.forEach((group, position) => {
          const outcome = batch.outcomes[position];
          group.functionalPass = outcome.passed; group.trimmedNt = outcome.nt; group.trimmedAa = outcome.aa;
          group.referenceMatch = outcome.referenceMatch;
          group.functionalRejectionReasons = outcome.reasons;
          if (outcome.passed) {
            passedCount++;
            if (outcome.alignedNt && outcome.alignedAa) {
              const name = functionalSequenceName(group);
              functionalNucleotideRows.push({ name, sequence: outcome.alignedNt });
              functionalProteinRows.push({ name, sequence: outcome.alignedAa });
            }
          }
        });
        functionalPassed = passedCount;
        if (passedCount) {
          alignments[`${sample.name}/functional-nucleotide`] = fasta(functionalNucleotideRows);
          alignments[`${sample.name}/functional-protein`] = fasta(functionalProteinRows);
          referenceAlignments[`${sample.name}/functional-nucleotide`] = fasta([{ name: batch.referenceName, sequence: batch.referenceNt }]);
        }
      } else if (sample.functionalReferenceSequence) functionalPassed = 0;
    }
    const summary = summaries.find((row) => row.sample === sample.name);
    if (summary) { summary.collapsedSequences = collapsedCount; summary.functionalPassed = functionalPassed; }
    onProgress?.({ fraction: (index + 1) / Math.max(1, activeSamples.length),
      detail: `Collapsed ${sample.name} into ${collapsedCount.toLocaleString()} variants${functionalPassed === undefined ? "" : `; ${functionalPassed.toLocaleString()} passed functional filtering`}; counts represent UMI families` });
  }
  return { ...output, summaries, alignments, referenceAlignments, collapseGroups, functionalFilterErrors,
    collapseSeconds: output.collapseSeconds + (performance.now() - started) / 1000 };
}
