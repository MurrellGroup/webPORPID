const MAGIC = [0x41, 0x4d, 0x53, 0x41] as const;
const encoder = new TextEncoder();
const decoder = new TextDecoder("utf-8", { fatal: true });

export function encodeAlivibeMsaSequences(sequences: readonly string[]): ArrayBuffer {
  if (sequences.length > 10_000) throw new Error("Alivibe MSA accepts at most 10,000 rows.");
  const rows = sequences.map((sequence, row) => {
    if (!/^[\x21-\x7e]*$/.test(sequence)) throw new Error(`Alivibe MSA row ${row + 1} contains a non-ASCII character.`);
    return encoder.encode(sequence);
  });
  const byteLength = 8 + rows.reduce((total, row) => total + 4 + row.byteLength, 0);
  if (byteLength > 256 * 1024 * 1024) throw new Error("Alivibe MSA input is too large.");
  const buffer = new ArrayBuffer(byteLength), bytes = new Uint8Array(buffer), view = new DataView(buffer);
  bytes.set(MAGIC); view.setUint32(4, rows.length, true); let offset = 8;
  for (const row of rows) { view.setUint32(offset, row.byteLength, true); offset += 4; bytes.set(row, offset); offset += row.byteLength; }
  return buffer;
}

export function decodeAlivibeMsaSequences(buffer: ArrayBuffer): string[] {
  const bytes = new Uint8Array(buffer);
  if (bytes.byteLength < 8 || MAGIC.some((value, index) => bytes[index] !== value)) throw new Error("Alivibe MSA returned an invalid result.");
  const view = new DataView(buffer), count = view.getUint32(4, true), rows: string[] = []; let offset = 8;
  if (count > 10_000) throw new Error("Alivibe MSA returned too many rows.");
  for (let row = 0; row < count; row++) {
    if (offset + 4 > bytes.byteLength) throw new Error("Alivibe MSA returned a truncated result.");
    const length = view.getUint32(offset, true); offset += 4;
    if (offset + length > bytes.byteLength) throw new Error("Alivibe MSA returned a truncated row.");
    rows.push(decoder.decode(bytes.subarray(offset, offset + length))); offset += length;
  }
  if (offset !== bytes.byteLength) throw new Error("Alivibe MSA returned trailing bytes."); return rows;
}

export function assertAlivibeMsaResult(input: readonly string[], aligned: readonly string[]) {
  if (aligned.length !== input.length) throw new Error("Alivibe MSA returned the wrong number of rows.");
  const width = aligned[0]?.length ?? 0;
  aligned.forEach((sequence) => { if (sequence.length !== width) throw new Error("Alivibe MSA returned a non-rectangular alignment."); });
}
