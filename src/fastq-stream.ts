import { sha256 } from "@noble/hashes/sha2.js";
import { bytesToHex } from "@noble/hashes/utils.js";
import { Gunzip } from "fflate";

export interface FastqBatch { text: string; count: number; firstOrdinal: number }
export interface FastqStreamProgress { compressedBytes: number; totalBytes: number; records: number; maxCarryBytes: number; maxBatchBytes: number }

export type StreamingHash = ReturnType<typeof sha256.create>;
export const createStreamingHash = (): StreamingHash => sha256.create();
export const finishStreamingHash = (hash: StreamingHash) => bytesToHex(hash.digest());

function gzipTransform(): TransformStream<Uint8Array, Uint8Array> {
  if ("DecompressionStream" in globalThis) return new DecompressionStream("gzip") as TransformStream<Uint8Array, Uint8Array>;
  let output: TransformStreamDefaultController<Uint8Array>;
  const gunzip = new Gunzip((chunk) => output.enqueue(chunk));
  return new TransformStream({
    start(controller) { output = controller; },
    transform(chunk) { gunzip.push(chunk, false); },
    flush() { gunzip.push(new Uint8Array(), true); },
  });
}

export async function* streamFastq(
  file: File,
  hash: StreamingHash,
  options: { recordsPerBatch?: number; maximumBatchBytes?: number; signal?: AbortSignal; onProgress?: (progress: FastqStreamProgress) => void } = {},
): AsyncGenerator<FastqBatch> {
  const recordsPerBatch = options.recordsPerBatch ?? 256;
  const maximumBatchBytes = options.maximumBatchBytes ?? 4 * 1024 * 1024;
  let compressedBytes = 0, records = 0, maxCarryBytes = 0, maxBatchBytes = 0;
  const measured = file.stream().pipeThrough(new TransformStream<Uint8Array, Uint8Array>({
    transform(chunk, controller) {
      if (options.signal?.aborted) throw new DOMException("Analysis cancelled.", "AbortError");
      hash.update(chunk); compressedBytes += chunk.byteLength; controller.enqueue(chunk);
    },
  }));
  const bytes = file.name.toLowerCase().endsWith(".gz") ? measured.pipeThrough(gzipTransform()) : measured;
  const reader = bytes.getReader(), decoder = new TextDecoder();
  let carry = "", lines: string[] = [], batch = "", batchCount = 0, firstOrdinal = 1;

  const progress = () => options.onProgress?.({ compressedBytes, totalBytes: file.size, records, maxCarryBytes, maxBatchBytes });
  const consume = (text: string): FastqBatch[] => {
    const output: FastqBatch[] = []; carry += text; let start = 0;
    for (let index = 0; index < carry.length; index++) {
      if (carry.charCodeAt(index) !== 10) continue;
      let line = carry.slice(start, index); if (line.endsWith("\r")) line = line.slice(0, -1); lines.push(line); start = index + 1;
      if (lines.length === 4) {
        if (!lines[0].startsWith("@") || !lines[2].startsWith("+") || lines[1].length !== lines[3].length)
          throw new Error(`Malformed FASTQ record ${records + 1}.`);
        batch += `${lines[0]}\n${lines[1]}\n${lines[2]}\n${lines[3]}\n`; lines = []; batchCount++; records++;
        if (batchCount >= recordsPerBatch || batch.length >= maximumBatchBytes) {
          output.push({ text: batch, count: batchCount, firstOrdinal }); firstOrdinal += batchCount;
          maxBatchBytes = Math.max(maxBatchBytes, batch.length); batch = ""; batchCount = 0;
        }
      }
    }
    carry = carry.slice(start); maxCarryBytes = Math.max(maxCarryBytes, carry.length); return output;
  };

  try {
    while (true) {
      if (options.signal?.aborted) throw new DOMException("Analysis cancelled.", "AbortError");
      const { done, value } = await reader.read(); if (done) break;
      const text = decoder.decode(value, { stream: true });
      for (const ready of consume(text)) { progress(); yield ready; }
    }
    const tail = decoder.decode(); for (const ready of consume(tail)) yield ready;
    if (carry) { let line = carry; if (line.endsWith("\r")) line = line.slice(0, -1); lines.push(line); carry = ""; }
    if (lines.length) {
      if (lines.length !== 4 || !lines[0].startsWith("@") || !lines[2].startsWith("+") || lines[1].length !== lines[3].length)
        throw new Error("The FASTQ stream is truncated at its final record.");
      batch += `${lines[0]}\n${lines[1]}\n${lines[2]}\n${lines[3]}\n`; batchCount++; records++;
    }
    if (batchCount) { maxBatchBytes = Math.max(maxBatchBytes, batch.length); yield { text: batch, count: batchCount, firstOrdinal }; }
    progress();
  } finally { reader.releaseLock(); }
}
