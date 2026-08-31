import { runAlivibeMsa } from "./alivibe-msa-runtime.ts";
import { translateAlignedNucleotides } from "./alignment-utils.ts";
import { collapseAlignment } from "./collapse.ts";
import { extractAndScorePanel } from "./panel-profile.ts";
import { runScalableMsa } from "./scalable-msa.ts";
import type { ApobecResult, ConsensusRecord, ContaminationCall, PipelineConfig, PostprocRecord, SampleSummary } from "./types.ts";

export interface PostprocessOutput {
  records: PostprocRecord[];
  summaries: SampleSummary[];
  alignments: Record<string, string>;
  referenceAlignments: Record<string, string>;
  collapseGroups: Record<string, ReturnType<typeof collapseAlignment>["groups"]>;
  collapseSeconds: number;
}

export type MsaRunner = (sequences: readonly string[], signal?: AbortSignal, iterations?: number,
  scoringMode?: "literal" | "nucleotide" | "amino-acid") => Promise<string[]>;

export interface PostprocessProgress { fraction: number; detail: string }
export interface PostprocessOptions {
  /** Keep false when collapse is run as its own cancellable pipeline stage. */
  collapse?: boolean;
  onCollapseProgress?: (progress: PostprocessProgress) => void;
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
  return best ? sequence.slice(best.start, best.end) : sequence.slice(0, Math.floor(sequence.length / 3) * 3);
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
  const coding = longestOrf(sequence); if (!coding || coding.length % 3) return { passed: false, reasons: ["frameshift-reject"] };
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
  return { passed: !reasons.length, reasons, nt: trimmed, aa };
}

interface FunctionalOutcome {
  passed: boolean;
  reasons: string[];
  nt?: string;
  aa?: string;
  /** Codon-aware alignment, clipped to the first/last reference residue. */
  alignedNt?: string;
  alignedAa?: string;
}

interface FunctionalBatchOutcome {
  outcomes: FunctionalOutcome[];
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

/**
 * SeededAlignment's functional filter is codon-aware and evaluates every
 * sequence against one joint reference-anchored alignment.  Doing independent
 * nucleotide NW alignments introduced arbitrary one-base gaps and incorrectly
 * labeled almost every long-read sequence as a frameshift.  Translate first,
 * align the complete batch, and project gaps back as codon triplets instead.
 */
async function functionalFilterBatch(
  reference: string, sequences: string[], threshold: number, runMsa: MsaRunner, signal?: AbortSignal,
): Promise<FunctionalBatchOutcome> {
  const outcomes: FunctionalOutcome[] = Array(sequences.length), coding: Array<{ index: number; sequence: string }> = [];
  sequences.forEach((raw, index) => {
    const sequence = degap(raw);
    if (/[^ACGT]/.test(sequence)) outcomes[index] = { passed: false, reasons: ["ambiguousSymbols-reject"] };
    else coding.push({ index, sequence: longestOrf(sequence) });
  });
  const referenceCoding = degap(reference).slice(0, Math.floor(degap(reference).length / 3) * 3);
  if (!coding.length) return { outcomes, referenceNt: referenceCoding, referenceAa: translate(referenceCoding) };
  const aminoInputs = [translate(referenceCoding), ...coding.map((row) => translate(row.sequence))];
  const aminoAlignment = await runScalableMsa(aminoInputs, runMsa, signal, 3, "amino-acid");
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
      alignedNt: queryRegion, alignedAa: translateAlignedNucleotides(queryRegion, 0) };
  });
  return { outcomes, referenceNt: referenceRegion, referenceAa: translateAlignedNucleotides(referenceRegion, 0) };
}

export async function postprocess(
  consensuses: ConsensusRecord[], contamination: ContaminationCall[], config: PipelineConfig, signal?: AbortSignal,
  runMsa: MsaRunner = runAlivibeMsa, sampleConcurrency = 1, onProgress?: (progress: PostprocessProgress) => void,
  options: PostprocessOptions = {},
): Promise<PostprocessOutput> {
  const discarded = new Set(contamination.filter((call) => call.discarded).map((call) => call.sequenceId));
  const outputs: PostprocessOutput[] = Array(config.samples.length); let cursor = 0;
  const sourceBySample = Array.from({ length: config.samples.length }, () => [] as ConsensusRecord[]);
  const sampleIndexByName = new Map(config.samples.map((sample, index) => [sample.name, index]));
  for (const record of consensuses) {
    const sampleIndex = config.samples[record.sampleIndex]?.name === record.sample ? record.sampleIndex : sampleIndexByName.get(record.sample);
    if (sampleIndex !== undefined) sourceBySample[sampleIndex].push(record);
  }
  const sampleProgress = Array(config.samples.length).fill(0);
  const report = (sampleIndex: number, fraction: number, detail: string) => {
    sampleProgress[sampleIndex] = Math.max(sampleProgress[sampleIndex], Math.max(0, Math.min(1, fraction)));
    onProgress?.({ fraction: sampleProgress.reduce((sum, value) => sum + value, 0) / Math.max(1, sampleProgress.length), detail });
  };
  await Promise.all(Array.from({ length: Math.min(config.samples.length, Math.max(1, Math.floor(sampleConcurrency))) }, async () => {
    while (true) {
      const sampleIndex = cursor++; if (sampleIndex >= config.samples.length) return;
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
      report(sampleIndex, .18, `Screening candidate sequences against the reference panel for sample ${sample.name}`);
      const candidates = preliminary.length > 1
        ? await runScalableMsa(preliminary.map(({ record }) => degap(record.sequence)), runMsa, signal, 3, "nucleotide")
        : preliminary.map(({ record }) => degap(record.sequence));
      // The supplied panel is already a reference alignment in PORPID. Align
      // its profile to the independently aligned sample profile, matching the reference workflow.
      const panelResult = extractAndScorePanel(candidates, sample.panelSequences.map((record) => record.sequence));
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
    const functionalNucleotideRows: Array<{ name: string; sequence: string }> = [];
    const functionalProteinRows: Array<{ name: string; sequence: string }> = [];
    const functionalByIndex = new Map<number, FunctionalOutcome>();
    let functionalReferenceNt = "", functionalReferenceAa = "";
    if (sample.functionalReferenceSequence && accepted.length) {
      report(sampleIndex, .72, `Checking coding-frame and functional-reference requirements for sample ${sample.name}`);
      const batch = await functionalFilterBatch(sample.functionalReferenceSequence.sequence,
        accepted.map(({ index }) => extracted.get(index)!), sample.functionalMatchOverride ?? config.parameters.functionalMatchThreshold, runMsa, signal);
      functionalReferenceNt = batch.referenceNt; functionalReferenceAa = batch.referenceAa;
      accepted.forEach(({ index }, position) => functionalByIndex.set(index, batch.outcomes[position]));
    }
    report(sampleIndex, .84, `Calculating sequence annotations and filter decisions for sample ${sample.name}`);
    let functionalPassed = 0;
    for (const [index, record] of source.entries()) {
      if (signal?.aborted) throw new DOMException("Downstream filtering skipped.", "AbortError");
      const artefactPass = record.familySize >= artefactCutoff, agreementPass = record.minimumAgreement >= agreementThreshold;
      const contaminationPass = !discarded.has(record.id), acceptedRow = alignmentByIndex.get(index), rejectionReasons: string[] = [];
      if (!artefactPass) rejectionReasons.push(`ccs_count < artefact cutoff (${artefactCutoff})`);
      if (!agreementPass) rejectionReasons.push(`minimum_agreement < ${agreementThreshold}`);
      if (!contaminationPass) rejectionReasons.push("contamination filter"); if (!panelPass[index]) rejectionReasons.push(`distance_from_panel >= ${config.parameters.panelThreshold}`);
      let trimmedNt: string | undefined, trimmedAa: string | undefined, functionalPass: boolean | undefined;
      if (sample.functionalReferenceSequence) {
        if (acceptedRow) {
          const outcome = functionalByIndex.get(index)!;
          functionalPass = outcome.passed; trimmedNt = outcome.nt; trimmedAa = outcome.aa; rejectionReasons.push(...outcome.reasons);
          if (outcome.passed) {
            functionalPassed++;
            if (outcome.alignedNt && outcome.alignedAa) {
              functionalNucleotideRows.push({ name: record.id, sequence: outcome.alignedNt });
              functionalProteinRows.push({ name: record.id, sequence: outcome.alignedAa });
            }
          }
        }
      }
      if (acceptedRow) nucleotideRows.push({ name: record.id, sequence: acceptedRow });
      records.push({ id: record.id, sample: sample.name, umi: record.umi, familySize: record.familySize,
        minimumAgreement: record.minimumAgreement, consensusNt: record.sequence, alignedNt: acceptedRow, trimmedNt, trimmedAa,
        panelScore: scores[index], artefactPass, agreementPass, contaminationPass, panelPass: panelPass[index], functionalPass,
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
    if (functionalNucleotideRows.length) {
      // The joint codon-aware alignment is sliced at the first and last
      // non-gap reference nucleotide. Query-only flanks therefore cannot leak
      // into functional exports or the functional tree/alignment workbench.
      alignments[`${sample.name}/functional-nucleotide`] = fasta(functionalNucleotideRows);
      alignments[`${sample.name}/functional-protein`] = fasta(functionalProteinRows);
      referenceAlignments[`${sample.name}/functional-nucleotide`] = fasta([{ name: "functional_reference", sequence: functionalReferenceNt }]);
    }
    summaries.push({ sample: sample.name, demultiplexedReads: 0, observedUmis: 0, likelyRealUmis: 0,
      consensusSequences: source.length, contaminationPassed: source.length - discardedInSample,
      postprocPassed: accepted.length,
      functionalPassed: sample.functionalReferenceSequence ? functionalPassed : undefined, artefactCutoff });
      outputs[sampleIndex] = { records, summaries, alignments, referenceAlignments, collapseGroups, collapseSeconds: 0 };
      report(sampleIndex, 1, `Downstream processing complete for sample ${sample.name}`);
    }
  }));
  const combined = { records: outputs.flatMap((output) => output.records), summaries: outputs.flatMap((output) => output.summaries),
    alignments: Object.assign({}, ...outputs.map((output) => output.alignments)),
    referenceAlignments: Object.assign({}, ...outputs.map((output) => output.referenceAlignments)),
    collapseGroups: Object.assign({}, ...outputs.map((output) => output.collapseGroups)), collapseSeconds: 0 };
  return options.collapse === false ? combined : collapsePostprocess(combined, config, signal, options.onCollapseProgress);
}

/** Run family-count-preserving haplotype collapse as its own resumable stage. */
export async function collapsePostprocess(
  output: PostprocessOutput, config: PipelineConfig, signal?: AbortSignal,
  onProgress?: (progress: PostprocessProgress) => void,
): Promise<PostprocessOutput> {
  const started = performance.now(), alignments = { ...output.alignments }, referenceAlignments = { ...output.referenceAlignments };
  const collapseGroups = { ...output.collapseGroups }, summaries = output.summaries.map((summary) => ({ ...summary }));
  for (const [index, sample] of config.samples.entries()) {
    if (signal?.aborted) throw new DOMException("Haplotype collapse skipped.", "AbortError");
    onProgress?.({ fraction: index / Math.max(1, config.samples.length), detail: `Collapsing identical retained UMI-family sequences for ${sample.name}` });
    await new Promise<void>((resolvePromise) => setTimeout(resolvePromise, 0));
    const uncollapsed = alignments[`${sample.name}/uncollapsed-nucleotide`];
    let collapsedCount = 0;
    if (uncollapsed) {
      const collapsed = collapseAlignment(uncollapsed, sample.name); collapsedCount = collapsed.groups.length;
      alignments[`${sample.name}/nucleotide`] = collapsed.fasta;
      collapseGroups[sample.name] = collapsed.groups;
      const rows = new Map(uncollapsed.split(/^>/m).filter(Boolean).map((block) => {
        const [name, ...sequence] = block.trimEnd().split(/\r?\n/); return [name, sequence.join("")];
      }));
      alignments[`${sample.name}/protein`] = fasta(collapsed.groups.map((group) => ({ name: group.representativeId,
        sequence: translateAlignedNucleotides(rows.get(group.representativeId) ?? "", 0) })));
      const reference = referenceAlignments[`${sample.name}/uncollapsed-nucleotide`];
      if (reference) referenceAlignments[`${sample.name}/nucleotide`] = reference;
    }
    const summary = summaries.find((row) => row.sample === sample.name); if (summary) summary.collapsedSequences = collapsedCount;
    onProgress?.({ fraction: (index + 1) / Math.max(1, config.samples.length),
      detail: `Collapsed ${sample.name} into ${collapsedCount.toLocaleString()} haplotypes; counts represent UMI families` });
  }
  return { ...output, summaries, alignments, referenceAlignments, collapseGroups,
    collapseSeconds: output.collapseSeconds + (performance.now() - started) / 1000 };
}
