import type { MsaRunner } from "./postprocess";

const MONOLITHIC_ROWS = 8_000;
const MONOLITHIC_BASES = 128 * 1024 * 1024;
const PROFILE_BATCH_ROWS = 2_000;

interface AnchoredRow { insertions: string[]; bases: string[] }

function decompose(anchor: string, alignedAnchor: string, alignedRow: string): AnchoredRow {
  if (alignedAnchor.length !== alignedRow.length || alignedAnchor.replaceAll("-", "") !== anchor)
    throw new Error("A scalable MSA batch returned an invalid anchor alignment.");
  const insertions = Array.from({ length: anchor.length + 1 }, () => ""), bases = Array<string>(anchor.length);
  let position = 0;
  for (let column = 0; column < alignedAnchor.length; column++) {
    if (alignedAnchor[column] === "-") insertions[position] += alignedRow[column];
    else {
      if (alignedAnchor[column] !== anchor[position]) throw new Error("A scalable MSA batch changed the anchor sequence.");
      bases[position++] = alignedRow[column];
    }
  }
  if (position !== anchor.length) throw new Error("A scalable MSA batch truncated the anchor sequence.");
  return { insertions, bases };
}

function padInsertion(value: string, width: number, slot: number) {
  return slot === 0 ? value.padStart(width, "-") : value.padEnd(width, "-");
}

export async function runScalableMsa(
  sequences: readonly string[], runMsa: MsaRunner, signal?: AbortSignal, iterations = 3,
  scoringMode: "literal" | "nucleotide" | "amino-acid" = "literal",
) {
  if (sequences.length < 2) return [...sequences];
  const totalBases = sequences.reduce((sum, sequence) => sum + sequence.length, 0);
  if (sequences.length <= MONOLITHIC_ROWS && totalBases <= MONOLITHIC_BASES)
    return runMsa(sequences, signal, iterations, scoringMode);

  const anchor = sequences[0], widths = new Uint32Array(anchor.length + 1);
  const decomposed: AnchoredRow[] = Array(sequences.length);
  for (let start = 0; start < sequences.length; start += PROFILE_BATCH_ROWS) {
    if (signal?.aborted) throw new DOMException("Analysis cancelled.", "AbortError");
    const end = Math.min(sequences.length, start + PROFILE_BATCH_ROWS);
    const input = start === 0 ? sequences.slice(start, end) : [anchor, ...sequences.slice(start, end)];
    const aligned = await runMsa(input, signal, iterations, scoringMode);
    if (aligned.length !== input.length || aligned.some((row) => row.length !== aligned[0].length))
      throw new Error("A scalable MSA batch returned a non-rectangular alignment.");
    const alignedAnchor = aligned[0], first = start === 0 ? 0 : 1;
    for (let row = first; row < aligned.length; row++) {
      const index = start === 0 ? row : start + row - 1;
      if (aligned[row].replaceAll("-", "") !== sequences[index]) throw new Error("A scalable MSA batch changed an input sequence.");
      const parts = decompose(anchor, alignedAnchor, aligned[row]); decomposed[index] = parts;
      parts.insertions.forEach((value, slot) => { widths[slot] = Math.max(widths[slot], value.length); });
    }
  }

  return decomposed.map((row) => {
    if (!row) throw new Error("A scalable MSA batch omitted a sequence.");
    let output = padInsertion(row.insertions[0], widths[0], 0);
    for (let position = 0; position < anchor.length; position++)
      output += row.bases[position] + padInsertion(row.insertions[position + 1], widths[position + 1], position + 1);
    return output;
  });
}
