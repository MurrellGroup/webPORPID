import type { AlignmentEdit, NamedSequence, ResultBundle } from "./types";

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
  removedNucleotides: number;
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

function isSubsequence(candidate: string, original: string): boolean {
  let offset = 0;
  for (const base of candidate) {
    offset = original.indexOf(base, offset);
    if (offset < 0) return false;
    offset += 1;
  }
  return true;
}

/** Allow gap editing and row/base deletion, but never invented biological content. */
export function validateCorrectedAlignment(currentText: string, correctedText: string): AlignmentCorrectionInspection {
  // A run can legitimately contain a single passing consensus.  Alivibe and
  // manual correction must still round-trip that alignment even though tree
  // inference is degenerate for one row.
  const current = inspectAlignment(currentText, 1);
  const corrected = inspectAlignment(correctedText, 1);
  const currentByName = new Map(current.records.map((record) => [record.name, record.sequence]));
  const correctedNames = new Set(corrected.records.map((record) => record.name));
  const removedRows = current.records.filter((record) => !correctedNames.has(record.name)).map((record) => record.name);
  const added = corrected.records.filter((record) => !currentByName.has(record.name)).map((record) => record.name);
  if (added.length) throw new Error(`The corrected alignment contains unexpected or renamed rows: ${added.slice(0, 4).join(", ")}${added.length > 4 ? "…" : ""}.`);
  let removedNucleotides = 0;
  for (const record of corrected.records) {
    const original = ungapped(currentByName.get(record.name)!);
    const replacement = ungapped(record.sequence);
    if (!isSubsequence(replacement, original))
      throw new Error(`The ungapped sequence for ${record.name} contains a substitution, inserted base, or changed order.`);
    removedNucleotides += original.length - replacement.length;
  }
  return { ...corrected, removedRows, removedNucleotides };
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

export function alignmentKey(sample: string): string { return `${sample}/nucleotide`; }

export function effectiveAlignment(bundle: ResultBundle, sample: string): { fasta?: string; edit?: AlignmentEdit; frameOffset: AlignmentFrameOffset } {
  const key = alignmentKey(sample), edit = bundle.alignmentEdits?.[key];
  return { fasta: edit?.fasta ?? bundle.alignments[key], edit, frameOffset: edit?.frameOffset ?? 0 };
}
