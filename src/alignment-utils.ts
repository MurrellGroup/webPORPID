import type { AlignmentChangeSummary, AlignmentEdit, AlignmentRowChange, NamedSequence, ResultBundle } from "./types";

export type AlignmentFrameOffset = 0 | 1 | 2;

export interface AlignmentInspection {
  fasta: string;
  records: NamedSequence[];
  rows: number;
  columns: number;
  fingerprint: string;
}

export interface AlignmentCorrectionInspection extends AlignmentInspection {
  removedRows: string[];
  addedRows: string[];
  changedRows: string[];
  rowChanges: AlignmentRowChange[];
  removedNucleotides: number;
  insertedNucleotides: number;
  substitutedNucleotides: number;
}

const CODONS: Record<string, string> = {
  TTT: "F", TTC: "F", TTA: "L", TTG: "L", TCT: "S", TCC: "S", TCA: "S", TCG: "S",
  TAT: "Y", TAC: "Y", TAA: "*", TAG: "*", TGT: "C", TGC: "C", TGA: "*", TGG: "W",
  CTT: "L", CTC: "L", CTA: "L", CTG: "L", CCT: "P", CCC: "P", CCA: "P", CCG: "P",
  CAT: "H", CAC: "H", CAA: "Q", CAG: "Q", CGT: "R", CGC: "R", CGA: "R", CGG: "R",
  ATT: "I", ATC: "I", ATA: "I", ATG: "M", ACT: "T", ACC: "T", ACA: "T", ACG: "T",
  AAT: "N", AAC: "N", AAA: "K", AAG: "K", AGT: "S", AGC: "S", AGA: "R", AGG: "R",
  GTT: "V", GTC: "V", GTA: "V", GTG: "V", GCT: "A", GCC: "A", GCA: "A", GCG: "A",
  GAT: "D", GAC: "D", GAA: "E", GAG: "E", GGT: "G", GGC: "G", GGA: "G", GGG: "G",
};

function exactFasta(records: NamedSequence[]): string {
  return records.map((record) => `>${record.name}\n${record.sequence}\n`).join("");
}

// Kept deliberately small and local: this module is also loaded directly by
// the Node conformance tests, whereas the application config parser pulls in
// binary serialization machinery that requires a full TypeScript transform.
function parseAlignmentFasta(source: string): NamedSequence[] {
  const records: NamedSequence[] = [];
  let name = "", sequence = "";
  for (const rawLine of source.split(/\r?\n/)) {
    const line = rawLine.trim();
    if (!line || line.startsWith(";")) continue;
    if (line.startsWith(">")) {
      if (name) records.push({ name, sequence: sequence.toUpperCase() });
      name = line.slice(1).trim() || `sequence_${records.length + 1}`;
      sequence = "";
    } else {
      if (!name) throw new Error("FASTA sequence data appeared before its first header.");
      sequence += line.replace(/\s/g, "");
    }
  }
  if (name) records.push({ name, sequence: sequence.toUpperCase() });
  if (!records.length) throw new Error("The FASTA file contains no records.");
  return records;
}

function assertFastaAlphabet(text: string): void {
  for (const line of text.split(/\r?\n/)) {
    if (!line || line.startsWith(">")) continue;
    const unsupported = line.replace(/\s/g, "").match(/[^ACGTUNRYKMSWBDHV.\-]/i)?.[0];
    if (unsupported) throw new Error(`The alignment contains the unsupported nucleotide character ${JSON.stringify(unsupported)}.`);
  }
}

function fnv1a64(value: string): string {
  let hash = 0xcbf29ce484222325n;
  const mask = 0xffffffffffffffffn;
  for (const byte of new TextEncoder().encode(value)) {
    hash ^= BigInt(byte);
    hash = (hash * 0x100000001b3n) & mask;
  }
  return hash.toString(16).padStart(16, "0");
}

export function inspectAlignment(text: string, minimumRows = 2): AlignmentInspection {
  assertFastaAlphabet(text);
  const normalizedText = text.split(/\r?\n/).map((line) => line.startsWith(">") ? line : line.replace(/U/gi, "T")).join("\n");
  const records = parseAlignmentFasta(normalizedText);
  if (records.length < minimumRows) throw new Error(`The alignment must contain at least ${minimumRows} sequences.`);
  const columns = records[0]?.sequence.length ?? 0;
  if (!columns || records.some((record) => record.sequence.length !== columns))
    throw new Error("Every alignment record must have the same non-zero aligned length.");
  const names = new Set<string>();
  for (const record of records) {
    if (names.has(record.name)) throw new Error(`The alignment contains duplicate identifier ${record.name}.`);
    names.add(record.name);
  }
  const fasta = exactFasta(records);
  return { fasta, records, rows: records.length, columns, fingerprint: fnv1a64(fasta) };
}

function ungapped(sequence: string): string { return sequence.replaceAll("-", "").replaceAll("U", "T"); }

/** Exact unit-cost Levenshtein decomposition with deterministic diagonal ties. */
function nucleotideEditCounts(original: string, replacement: string) {
  if (original === replacement) return { substitutedNucleotides: 0, insertedNucleotides: 0, removedNucleotides: 0 };
  const width = replacement.length + 1;
  let previousCost = new Uint32Array(width), previousSubstitutions = new Uint32Array(width), previousInsertions = new Uint32Array(width), previousDeletions = new Uint32Array(width);
  let currentCost = new Uint32Array(width), currentSubstitutions = new Uint32Array(width), currentInsertions = new Uint32Array(width), currentDeletions = new Uint32Array(width);
  for (let column = 0; column < width; column += 1) { previousCost[column] = column; previousInsertions[column] = column; }
  for (let row = 1; row <= original.length; row += 1) {
    currentCost[0] = row; currentSubstitutions[0] = 0; currentInsertions[0] = 0; currentDeletions[0] = row;
    for (let column = 1; column < width; column += 1) {
      const mismatch = original[row - 1] === replacement[column - 1] ? 0 : 1;
      // Start with the diagonal candidate so a tied one-base substitution is
      // reported as such instead of an arbitrary delete/insert pair.
      let cost = previousCost[column - 1] + mismatch;
      let substitutions = previousSubstitutions[column - 1] + mismatch;
      let insertions = previousInsertions[column - 1], deletions = previousDeletions[column - 1];
      const deletionCost = previousCost[column] + 1;
      if (deletionCost < cost) {
        cost = deletionCost; substitutions = previousSubstitutions[column]; insertions = previousInsertions[column]; deletions = previousDeletions[column] + 1;
      }
      const insertionCost = currentCost[column - 1] + 1;
      if (insertionCost < cost) {
        cost = insertionCost; substitutions = currentSubstitutions[column - 1]; insertions = currentInsertions[column - 1] + 1; deletions = currentDeletions[column - 1];
      }
      currentCost[column] = cost; currentSubstitutions[column] = substitutions; currentInsertions[column] = insertions; currentDeletions[column] = deletions;
    }
    [previousCost, currentCost] = [currentCost, previousCost];
    [previousSubstitutions, currentSubstitutions] = [currentSubstitutions, previousSubstitutions];
    [previousInsertions, currentInsertions] = [currentInsertions, previousInsertions];
    [previousDeletions, currentDeletions] = [currentDeletions, previousDeletions];
  }
  const column = replacement.length;
  return { substitutedNucleotides: previousSubstitutions[column], insertedNucleotides: previousInsertions[column], removedNucleotides: previousDeletions[column] };
}

/**
 * Validate the FASTA shape without restricting biological edits. Substitutions,
 * insertions, deletions, renamed rows and removed rows are all retained in the
 * edit audit so the UI can warn before saving them.
 */
export function validateCorrectedAlignment(currentText: string, correctedText: string): AlignmentCorrectionInspection {
  // A run can legitimately contain a single passing consensus.  Alivibe and
  // manual correction must still round-trip that alignment even though tree
  // inference is degenerate for one row.
  const current = inspectAlignment(currentText, 1);
  const corrected = inspectAlignment(correctedText, 1);
  const currentByName = new Map(current.records.map((record) => [record.name, record.sequence]));
  const correctedNames = new Set(corrected.records.map((record) => record.name));
  const removedRows = current.records.filter((record) => !correctedNames.has(record.name)).map((record) => record.name);
  const addedRows = corrected.records.filter((record) => !currentByName.has(record.name)).map((record) => record.name);
  const changedRows: string[] = [], rowChanges: AlignmentRowChange[] = [];
  let removedNucleotides = 0, insertedNucleotides = 0, substitutedNucleotides = 0;
  for (const record of corrected.records) {
    const originalAligned = currentByName.get(record.name);
    if (originalAligned == null) continue;
    const original = ungapped(originalAligned);
    const replacement = ungapped(record.sequence);
    if (original !== replacement || originalAligned !== record.sequence) changedRows.push(record.name);
    const edits = nucleotideEditCounts(original, replacement);
    const rowSubstitutions = edits.substitutedNucleotides, rowRemoved = edits.removedNucleotides, rowInserted = edits.insertedNucleotides;
    substitutedNucleotides += rowSubstitutions; removedNucleotides += rowRemoved; insertedNucleotides += rowInserted;
    if (original !== replacement || originalAligned !== record.sequence) rowChanges.push({ name: record.name,
      substitutedNucleotides: rowSubstitutions, insertedNucleotides: rowInserted, removedNucleotides: rowRemoved,
      gapPlacementChanged: originalAligned !== record.sequence && original === replacement });
  }
  return { ...corrected, removedRows, addedRows, changedRows, rowChanges, removedNucleotides, insertedNucleotides, substitutedNucleotides };
}

export function summarizeAlignmentChanges(currentText: string, correctedText: string): AlignmentChangeSummary {
  const current = inspectAlignment(currentText, 1), corrected = validateCorrectedAlignment(currentText, correctedText);
  const currentNames = new Set(current.records.map((record) => record.name)), correctedNames = new Set(corrected.records.map((record) => record.name));
  const sharedBefore = current.records.map((record) => record.name).filter((name) => correctedNames.has(name));
  const sharedAfter = corrected.records.map((record) => record.name).filter((name) => currentNames.has(name));
  return {
    rowsBefore: current.rows, rowsAfter: corrected.rows, columnsBefore: current.columns, columnsAfter: corrected.columns,
    rowOrderChanged: sharedBefore.join("\0") !== sharedAfter.join("\0"),
    rowOrderBefore: current.records.map((record) => record.name), rowOrderAfter: corrected.records.map((record) => record.name),
    removedRows: corrected.removedRows, addedRows: corrected.addedRows, changedRows: corrected.changedRows, rowChanges: corrected.rowChanges,
    removedNucleotides: corrected.removedNucleotides, insertedNucleotides: corrected.insertedNucleotides,
    substitutedNucleotides: corrected.substitutedNucleotides,
  };
}

export function translateAlignedNucleotides(sequence: string, frameOffset: number = 0): string {
  const offset: AlignmentFrameOffset = frameOffset === 1 || frameOffset === 2 ? frameOffset : 0;
  let result = "";
  for (let index = offset; index < sequence.length; index += 3) {
    const codon = sequence.slice(index, index + 3).toUpperCase();
    result += codon === "---" ? "-" : codon.includes("-") || /[^ACGT]/.test(codon) || codon.length < 3 ? "X" : CODONS[codon] ?? "X";
  }
  return result;
}

export function translateAlignmentFasta(fasta: string, frameOffset: AlignmentFrameOffset = 0): string {
  const inspected = inspectAlignment(fasta, 1);
  return exactFasta(inspected.records.map((record) => ({ ...record, sequence: translateAlignedNucleotides(record.sequence, frameOffset) })));
}

export type AlignmentVariant = "collapsed" | "uncollapsed" | "functional";

export function alignmentKey(sample: string, variant: AlignmentVariant = "collapsed"): string {
  return variant === "collapsed" ? `${sample}/nucleotide`
    : variant === "uncollapsed" ? `${sample}/uncollapsed-nucleotide` : `${sample}/functional-nucleotide`;
}

export function effectiveAlignment(bundle: ResultBundle, sample: string, variant: AlignmentVariant = "collapsed"): { fasta?: string; edit?: AlignmentEdit; frameOffset: AlignmentFrameOffset; key: string } {
  const key = alignmentKey(sample, variant), edit = bundle.alignmentEdits?.[key];
  const legacyFallback = variant === "uncollapsed" && !bundle.collapseGroups?.[sample] ? bundle.alignments[alignmentKey(sample)] : undefined;
  return { fasta: edit?.fasta ?? bundle.alignments[key] ?? legacyFallback, edit, frameOffset: edit?.frameOffset ?? 0, key };
}
