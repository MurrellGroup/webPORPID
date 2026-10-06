import type { ProfileAddition } from "./independent-panel-filter.ts";

export type MafftScoringMode = "literal" | "nucleotide" | "amino-acid";
export interface MafftRequest { sequences: string[]; scoringMode: MafftScoringMode; profileRows?: number }

/** disttbfast flags used by MAFFT 7.520's FFT-NS-2 and --add wrappers. */
export function mafftArguments(mode: MafftScoringMode, profileRows?: number): string[] {
  return [
    "-q", "0", "-E", "2", "-V", "-1.53", "-s", "0.0", "-W", "6", "-O",
    "-C", "0-0", mode === "amino-acid" ? "-P" : "-D", "-b", "62", "-g", "0", "-f", "-1.53", "-Q", "100.0",
    "-h", "0", "-F", "-X", "0.1", "-x", "-1",
    // Existing profile rows come first; exactly one new reference comes last.
    ...(profileRows === undefined ? [] : ["-K", "-I", "1"]), "-i", "/input.fa",
  ];
}

export function numericMafftFasta(sequences: readonly string[], mode: MafftScoringMode): string {
  return sequences.map((sequence, index) => {
    // MAFFT drops '*' in protein input. Align it as X, then restore the exact
    // original residue at that position before codon backtranslation.
    const input = mode === "amino-acid" ? sequence.replaceAll("*", "X") : sequence;
    return `>${index}\n${input}\n`;
  }).join("");
}

export function parseMafftAlignment(source: string, expected: readonly string[], mode: MafftScoringMode): string[] {
  const rows = new Map<number, string>(); let index = -1, sequence = "";
  const finish = () => {
    if (index < 0) return;
    if (rows.has(index)) throw new Error("MAFFT returned a duplicate sequence identifier.");
    rows.set(index, sequence.toUpperCase());
  };
  for (const raw of source.split(/\r?\n/)) {
    const line = raw.trim(); if (!line) continue;
    if (line.startsWith(">")) {
      finish(); index = Number(line.slice(1).trim()); sequence = "";
      if (!Number.isSafeInteger(index) || index < 0 || index >= expected.length)
        throw new Error("MAFFT returned an unknown sequence identifier.");
    } else {
      if (index < 0) throw new Error("MAFFT returned sequence data before its first header.");
      sequence += line.replace(/\s/g, "");
    }
  }
  finish();
  if (rows.size !== expected.length) throw new Error("MAFFT returned the wrong number of sequences.");
  const aligned = expected.map((input, row) => {
    const original = input.replaceAll("-", "").toUpperCase(), output = rows.get(row)!;
    const normalized = mode === "amino-acid" ? original.replaceAll("*", "X") : original;
    if (output.replaceAll("-", "") !== normalized) throw new Error(`MAFFT changed sequence ${row + 1}.`);
    let cursor = 0;
    return [...output].map((residue) => residue === "-" ? "-" : original[cursor++]).join("");
  });
  const width = aligned[0]?.length ?? 0;
  if (!width || aligned.some((row) => row.length !== width)) throw new Error("MAFFT returned a non-rectangular alignment.");
  return aligned;
}

export function mafftProfileAddition(profile: readonly string[], aligned: readonly string[]): ProfileAddition {
  if (aligned.length !== profile.length + 1) throw new Error("MAFFT profile addition returned the wrong number of rows.");
  const sampleRows = aligned.slice(0, profile.length);
  const withoutGapColumns = (rows: readonly string[]) => {
    const columns = Array.from({ length: rows[0]?.length ?? 0 }, (_, column) => column)
      .filter((column) => rows.some((row) => row[column] !== "-"));
    return rows.map((row) => columns.map((column) => row[column]).join(""));
  };
  const before = withoutGapColumns(profile), after = withoutGapColumns(sampleRows);
  if (before.some((row, index) => row !== after[index]))
    throw new Error("MAFFT reference addition changed the existing sample protein alignment.");
  return { profileRows: sampleRows, sequence: aligned[profile.length] };
}
