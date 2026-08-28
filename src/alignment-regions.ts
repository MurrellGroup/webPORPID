import { translateAlignedNucleotides, type AlignmentFrameOffset } from "./alignment-utils.ts";

export type CoordinateUnits = "nt" | "aa";

export function modalSequence(sequences: readonly string[]): string {
  const columns = sequences[0]?.length ?? 0;
  let result = "";
  for (let column = 0; column < columns; column += 1) {
    const counts = new Map<string, { count: number; first: number }>();
    sequences.forEach((sequence, index) => {
      const value = sequence[column]?.toUpperCase() ?? "-", current = counts.get(value) ?? { count: 0, first: index };
      current.count += 1; counts.set(value, current);
    });
    result += [...counts].sort((left, right) => right[1].count - left[1].count || left[1].first - right[1].first)[0]?.[0] ?? "-";
  }
  return result;
}

function coordinateMap(aligned: string): Map<number, number> {
  const result = new Map<number, number>(); let coordinate = 0;
  for (let column = 0; column < aligned.length; column += 1) {
    if (aligned[column] === "-") continue;
    coordinate += 1; result.set(coordinate, column);
  }
  return result;
}

function parseCoordinates(text: string, maximum: number): number[] {
  const selected = new Set<number>();
  for (const raw of text.split(/[\s,;]+/).filter(Boolean)) {
    const range = raw.match(/^(\d+)-(\d+)$/);
    if (range) {
      const left = Number(range[1]), right = Number(range[2]);
      for (let value = Math.max(1, Math.min(left, right)); value <= Math.min(maximum, Math.max(left, right)); value += 1) selected.add(value);
    } else {
      const value = Number(raw); if (Number.isInteger(value) && value >= 1 && value <= maximum) selected.add(value);
    }
  }
  return [...selected].sort((left, right) => left - right);
}

/** Map 1-based reference coordinates to compact display columns. */
export function referenceDisplayColumns(
  alignedReferenceNt: string,
  selection: string,
  units: CoordinateUnits,
  alphabet: "nt" | "aa",
  frameOffset: AlignmentFrameOffset,
  displayLength: number,
): number[] {
  if (!selection.trim()) return Array.from({ length: displayLength }, (_, index) => index);
  const ntMap = coordinateMap(alignedReferenceNt), aaReference = translateAlignedNucleotides(alignedReferenceNt, frameOffset), aaMap = coordinateMap(aaReference);
  const source = units === "nt" ? ntMap : aaMap, coordinates = parseCoordinates(selection, Math.max(0, ...source.keys()));
  const selected = new Set<number>();
  for (const coordinate of coordinates) {
    const sourceColumn = source.get(coordinate); if (sourceColumn == null) continue;
    if (units === alphabet) selected.add(sourceColumn);
    else if (units === "nt") {
      if (sourceColumn >= frameOffset) selected.add(Math.floor((sourceColumn - frameOffset) / 3));
    } else {
      const start = frameOffset + sourceColumn * 3;
      for (let offset = 0; offset < 3; offset += 1) selected.add(start + offset);
    }
  }
  return [...selected].filter((column) => column >= 0 && column < displayLength).sort((left, right) => left - right);
}

export function referenceCoordinateLabels(alignedReferenceNt: string, alphabet: "nt" | "aa", frameOffset: AlignmentFrameOffset): string[] {
  const aligned = alphabet === "nt" ? alignedReferenceNt : translateAlignedNucleotides(alignedReferenceNt, frameOffset);
  let coordinate = 0;
  return [...aligned].map((value) => value === "-" ? "" : String(++coordinate));
}
