import type { AlignmentFrameOffset } from "./alignment-utils";

export const ALIVIBE_BRIDGE_VERSION = 4;
export const ALIVIBE_OLDEST_COMPATIBLE_SNAPSHOT_VERSION = 3;
export const ALIVIBE_SOURCE_REVISION = "cbcd02719dd0a5f1f05d3127666f00e8579f2423";

export interface AlivibeNucleotideRecord {
  name: string;
  sequence: string;
}

export interface AlivibeNucleotideSnapshot {
  version: number;
  sourceRevision: string;
  alphabet: string;
  mode: string;
  frameOffset: number;
  records: AlivibeNucleotideRecord[];
  fasta: string;
}

export interface AlivibeSwigBridge {
  version: number;
  sourceRevision: string;
  loadNucleotideFasta: (text: string, frameOffset?: number) => unknown;
  snapshotNucleotide: () => unknown;
  installMsaRunner: (runner: AlivibeMsaRunner) => void;
  createMsaJob: (sequences: string[]) => AlivibeMsaJob | null;
  installFastTreeRunner: (runner: (fasta: string, alphabet: "nt" | "aa") => Promise<string>) => void;
  runFastTree: (fasta: string) => Promise<string> | null;
}

export interface AlivibeMsaJob {
  result: Promise<string[]>;
  cancel: () => void;
}

export type AlivibeMsaRunner = (sequences: string[], scoringMode: "nucleotide" | "amino-acid") => AlivibeMsaJob;

export interface AlivibeEditorWindow extends Window {
  swigAlivibeBridge?: AlivibeSwigBridge;
}

export interface AlivibeNucleotideTransfer {
  fasta: string;
  frameOffset: AlignmentFrameOffset;
  records: AlivibeNucleotideRecord[];
  sourceRevision: string;
}

export interface AlivibeRoundTripTarget {
  groupKey: string;
  alignmentFingerprint: string;
}

function exactFasta(records: AlivibeNucleotideRecord[]): string {
  return records.map((record) => `>${record.name}\n${record.sequence}\n`).join("");
}

function validateSnapshot(value: unknown): AlivibeNucleotideTransfer {
  if (!value || typeof value !== "object") throw new Error("Alivibe returned no nucleotide snapshot.");
  const snapshot = value as Partial<AlivibeNucleotideSnapshot>;
  const compatibleSnapshot = typeof snapshot.version === "number"
    && snapshot.version >= ALIVIBE_OLDEST_COMPATIBLE_SNAPSHOT_VERSION
    && snapshot.version <= ALIVIBE_BRIDGE_VERSION;
  if (!compatibleSnapshot || snapshot.sourceRevision !== ALIVIBE_SOURCE_REVISION) {
    throw new Error("The open Alivibe editor does not match this webPORPID release. Close it and reopen it from webPORPID.");
  }
  if (snapshot.alphabet !== "NT" || snapshot.mode !== "NT") {
    throw new Error("Alivibe did not return its nucleotide view. The alignment was not imported.");
  }
  if (snapshot.frameOffset !== 0 && snapshot.frameOffset !== 1 && snapshot.frameOffset !== 2) {
    throw new Error("Alivibe returned an invalid amino-acid frame offset.");
  }
  if (!Array.isArray(snapshot.records) || snapshot.records.length === 0) {
    throw new Error("Alivibe returned an empty nucleotide alignment.");
  }

  const seen = new Set<string>();
  const records = snapshot.records.map((record, index) => {
    if (!record || typeof record !== "object") throw new Error(`Alivibe row ${index + 1} is invalid.`);
    const name = (record as AlivibeNucleotideRecord).name;
    const sequence = (record as AlivibeNucleotideRecord).sequence;
    if (typeof name !== "string" || name.length === 0 || /[\r\n]/.test(name)) {
      throw new Error(`Alivibe row ${index + 1} has an invalid identifier.`);
    }
    if (seen.has(name)) throw new Error(`Alivibe returned duplicate identifier ${name}.`);
    seen.add(name);
    if (typeof sequence !== "string" || !/^[A-Z-]*$/.test(sequence)) {
      throw new Error(`Alivibe row ${name} contains characters outside its nucleotide alignment alphabet.`);
    }
    return { name, sequence };
  });

  const fasta = exactFasta(records);
  if (snapshot.fasta !== fasta) {
    throw new Error("Alivibe's nucleotide records and nucleotide FASTA disagree. The alignment was not imported.");
  }
  return {
    fasta,
    frameOffset: snapshot.frameOffset,
    records,
    sourceRevision: snapshot.sourceRevision,
  };
}

export function getAlivibeBridge(editor: AlivibeEditorWindow | null | undefined): AlivibeSwigBridge | undefined {
  const bridge = editor?.swigAlivibeBridge;
  if (
    bridge?.version !== ALIVIBE_BRIDGE_VERSION
    || bridge.sourceRevision !== ALIVIBE_SOURCE_REVISION
    || typeof bridge.loadNucleotideFasta !== "function"
    || typeof bridge.snapshotNucleotide !== "function"
    || typeof bridge.installMsaRunner !== "function"
    || typeof bridge.createMsaJob !== "function"
    || typeof bridge.installFastTreeRunner !== "function"
    || typeof bridge.runFastTree !== "function"
  ) return undefined;
  return bridge;
}

function getCompatibleSnapshotBridge(editor: AlivibeEditorWindow | null | undefined): Pick<AlivibeSwigBridge, "snapshotNucleotide"> | undefined {
  const bridge = editor?.swigAlivibeBridge;
  if (
    typeof bridge?.version !== "number"
    || bridge.version < ALIVIBE_OLDEST_COMPATIBLE_SNAPSHOT_VERSION
    || bridge.version > ALIVIBE_BRIDGE_VERSION
    || bridge.sourceRevision !== ALIVIBE_SOURCE_REVISION
    || typeof bridge.snapshotNucleotide !== "function"
  ) return undefined;
  return bridge;
}

export function assertAlivibeInitialLoad(expectedFasta: string, loaded: AlivibeNucleotideTransfer): void {
  if (loaded.fasta !== expectedFasta) {
    throw new Error("Alivibe did not load the exact webPORPID nucleotide alignment. The round trip was stopped.");
  }
}

export function assertAlivibeRoundTripTarget(origin: AlivibeRoundTripTarget, current: AlivibeRoundTripTarget): void {
  if (current.groupKey !== origin.groupKey) {
    throw new Error("The selected lineage changed after Alivibe was opened. Return to the originating lineage selection, or close Alivibe and open the current alignment again.");
  }
  if (current.alignmentFingerprint !== origin.alignmentFingerprint) {
    throw new Error("The Swig alignment changed after Alivibe was opened. Close Alivibe and reopen the current alignment to avoid overwriting newer work.");
  }
}

/** Load the exact nucleotide FASTA into the bundled Alivibe editor and verify its visible NT state. */
export function loadAlivibeNucleotideFasta(
  editor: AlivibeEditorWindow,
  fasta: string,
  frameOffset: AlignmentFrameOffset,
): AlivibeNucleotideTransfer {
  const bridge = getAlivibeBridge(editor);
  if (!bridge) throw new Error("The bundled Alivibe bridge is not ready.");
  return validateSnapshot(bridge.loadNucleotideFasta(fasta, frameOffset));
}

/** Read exactly the complete, ordered nucleotide rows used by Alivibe's NT viewer and NT export. */
export function readAlivibeNucleotideFasta(editor: AlivibeEditorWindow): AlivibeNucleotideTransfer {
  const bridge = getCompatibleSnapshotBridge(editor);
  if (!bridge) throw new Error("The open Alivibe editor cannot provide a compatible nucleotide snapshot.");
  return validateSnapshot(bridge.snapshotNucleotide());
}
