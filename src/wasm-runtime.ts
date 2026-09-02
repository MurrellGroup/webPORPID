import { WASI } from "@bjorn3/browser_wasi_shim";
import { BinaryReader, BinaryWriter } from "./binary";
import type { ConsensusRecord, FamilyDisposition, PipelineConfig, QualityStats, UmiFamily } from "./types";

interface CoreExports extends WebAssembly.Exports {
  memory: WebAssembly.Memory;
  wpp_alloc(size: number): number;
  wpp_free(pointer: number): void;
  wpp_init_config(pointer: number, length: number): number;
  wpp_preprocess(pointer: number, length: number, firstOrdinal: number): number;
  wpp_partition_counts(pointer: number, length: number): number;
  wpp_count_families(pointer: number, length: number, cutoffs: number, cutoffsLength: number): number;
  wpp_build_family_model(pointer: number, length: number): number;
  wpp_init_family_model(pointer: number, length: number): number;
  wpp_consensus_partition(pointer: number, length: number, cutoffs: number, cutoffsLength: number): number;
  wpp_stats(): number;
  wpp_result_ptr(): number;
  wpp_result_len(): number;
  wpp_error_ptr(): number;
  wpp_error_len(): number;
}

const encoder = new TextEncoder();
const decoder = new TextDecoder();

export class WebPorpidRuntime {
  private constructor(private core: CoreExports) {}

  static async create(module: WebAssembly.Module, compiledConfig: Uint8Array) {
    const wasi = new WASI([], [], []);
    const instance = await WebAssembly.instantiate(module, { wasi_snapshot_preview1: wasi.wasiImport });
    wasi.initialize(instance as WebAssembly.Instance & { exports: { memory: WebAssembly.Memory; _initialize?: () => unknown } });
    const runtime = new WebPorpidRuntime(instance.exports as CoreExports);
    runtime.withOne(compiledConfig, (pointer, length) => runtime.core.wpp_init_config(pointer, length), false);
    return runtime;
  }

  private errorText() {
    return decoder.decode(new Uint8Array(this.core.memory.buffer, this.core.wpp_error_ptr(), this.core.wpp_error_len())) || "The webPORPID WASM core failed.";
  }

  private resultBytes() {
    return new Uint8Array(this.core.memory.buffer, this.core.wpp_result_ptr(), this.core.wpp_result_len()).slice();
  }

  private put(bytes: Uint8Array): number {
    const pointer = this.core.wpp_alloc(bytes.byteLength);
    if (!pointer && bytes.byteLength) throw new Error("The webPORPID WASM core ran out of memory.");
    new Uint8Array(this.core.memory.buffer, pointer, bytes.byteLength).set(bytes); return pointer;
  }

  private withOne(bytes: Uint8Array, call: (pointer: number, length: number) => number, result = true): Uint8Array {
    const pointer = this.put(bytes);
    try { if (call(pointer, bytes.byteLength) < 0) throw new Error(this.errorText()); return result ? this.resultBytes() : new Uint8Array(); }
    finally { this.core.wpp_free(pointer); }
  }

  private withTwo(first: Uint8Array, second: Uint8Array, call: (a: number, al: number, b: number, bl: number) => number): Uint8Array {
    const a = this.put(first), b = this.put(second);
    try { if (call(a, first.byteLength, b, second.byteLength) < 0) throw new Error(this.errorText()); return this.resultBytes(); }
    finally { this.core.wpp_free(a); this.core.wpp_free(b); }
  }

  preprocess(fastq: string, firstOrdinal: number) {
    return this.withOne(encoder.encode(fastq), (pointer, length) => this.core.wpp_preprocess(pointer, length, firstOrdinal));
  }
  partitionCounts(partition: Uint8Array) {
    return this.withOne(partition, (pointer, length) => this.core.wpp_partition_counts(pointer, length));
  }
  countFamilies(partition: Uint8Array, cutoffs: Uint8Array) {
    return this.withTwo(partition, cutoffs, (a, al, b, bl) => this.core.wpp_count_families(a, al, b, bl));
  }
  buildFamilyModel(counts: Uint8Array) {
    return this.withOne(counts, (pointer, length) => this.core.wpp_build_family_model(pointer, length));
  }
  initFamilyModel(model: Uint8Array) {
    this.withOne(model, (pointer, length) => this.core.wpp_init_family_model(pointer, length), false);
  }
  consensus(partition: Uint8Array, cutoffs: Uint8Array) {
    return this.withTwo(partition, cutoffs, (a, al, b, bl) => this.core.wpp_consensus_partition(a, al, b, bl));
  }
  stats(): QualityStats {
    if (this.core.wpp_stats() < 0) throw new Error(this.errorText());
    const parsed = JSON.parse(decoder.decode(this.resultBytes())) as Omit<QualityStats, "downsampledReads">;
    return { ...parsed, downsampledReads: 0 };
  }
}

export function decodeBrowserFrames(bytes: Uint8Array): Array<{ partition: number; record: Uint8Array }> {
  const output: Array<{ partition: number; record: Uint8Array }> = [];
  const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength); let offset = 0;
  while (offset < bytes.byteLength) {
    if (offset + 5 > bytes.byteLength) throw new Error("Truncated spool routing frame.");
    const partition = bytes[offset++]; const length = view.getUint32(offset, true); offset += 4;
    if (offset + length > bytes.byteLength) throw new Error("Invalid spool routing frame length.");
    output.push({ partition, record: bytes.slice(offset, offset + length) }); offset += length;
  }
  return output;
}

export function decodePartitionCounts(bytes: Uint8Array): bigint[] {
  const reader = new BinaryReader(bytes); reader.magic("WPS1"); const count = reader.u32();
  const output = Array.from({ length: count }, () => reader.u64());
  if (!reader.done) throw new Error("Partition count payload has trailing bytes."); return output;
}

export function mergeSampleCounts(parts: bigint[][], samples: number): bigint[] {
  const output = Array.from({ length: samples }, () => 0n);
  for (const part of parts) for (let index = 0; index < part.length; index++) output[index] += part[index];
  return output;
}

export function makeCutoffValues(sampleCounts: bigint[], maximum: number): bigint[] {
  const maxHash = (1n << 64n) - 1n;
  return sampleCounts.map((count) => {
    const cap = maximum < 1 ? count : BigInt(maximum);
    return count === 0n || count <= cap ? maxHash : (maxHash * cap) / count;
  });
}

export function makeCutoffs(sampleCounts: bigint[], maximum: number): Uint8Array {
  const writer = new BinaryWriter(), values = makeCutoffValues(sampleCounts, maximum);
  writer.magic("WPT1"); writer.u32(values.length);
  for (const value of values) writer.u64(value);
  return writer.finish();
}

type CountEntry = { sample: number; umi: string; count: number };
export function decodeFamilyCounts(bytes: Uint8Array): CountEntry[] {
  const reader = new BinaryReader(bytes); reader.magic("WPN1"); const count = reader.u32(); const output: CountEntry[] = [];
  for (let index = 0; index < count; index++) output.push({ sample: reader.u16(), umi: reader.string(), count: reader.u32() });
  if (!reader.done) throw new Error("Family count payload has trailing bytes."); return output;
}

export function mergeFamilyCounts(parts: Uint8Array[]): Uint8Array {
  const counts = new Map<string, CountEntry>();
  for (const bytes of parts) for (const entry of decodeFamilyCounts(bytes)) {
    const key = `${entry.sample}\0${entry.umi}`, previous = counts.get(key);
    if (previous) previous.count += entry.count; else counts.set(key, { ...entry });
  }
  const entries = [...counts.values()].sort((a, b) => a.sample - b.sample || a.umi.localeCompare(b.umi));
  const writer = new BinaryWriter(); writer.magic("WPN1"); writer.u32(entries.length);
  for (const entry of entries) { writer.u16(entry.sample); writer.string(entry.umi); writer.u32(entry.count); }
  return writer.finish();
}

const DISPOSITIONS: FamilyDisposition[] = ["likely_real", "BPB-rejects", "heteroduplex", "LDA-rejects", "UMI_len != 8", "family-size-reject"];

export function decodeFamilyModel(bytes: Uint8Array, config: PipelineConfig): UmiFamily[] {
  const reader = new BinaryReader(bytes); reader.magic("WPM1"); const count = reader.u32(); const output: UmiFamily[] = [];
  for (let index = 0; index < count; index++) {
    const sampleIndex = reader.u16(), umi = reader.string(), parent = reader.string(), familySize = reader.u32();
    const probability = reader.f64(), disposition = DISPOSITIONS[reader.u8()];
    if (!disposition || !config.samples[sampleIndex]) throw new Error("Family model contains an invalid sample or disposition.");
    output.push({ sample: config.samples[sampleIndex].name, sampleIndex, umi: disposition === "BPB-rejects" ? "REJECTED" : umi, familySize, mostLikelyParent: parent,
      posteriorProbability: probability, logOffspringProbability: Math.log(1 - probability), disposition });
  }
  if (!reader.done) throw new Error("Family model has trailing bytes."); return output;
}

/** Re-encode an inspected family model after an interactive cutoff decision. */
export function encodeFamilyModel(families: readonly UmiFamily[]): Uint8Array {
  const writer = new BinaryWriter(); writer.magic("WPM1"); writer.u32(families.length);
  for (const family of families) {
    const disposition = DISPOSITIONS.indexOf(family.disposition);
    if (disposition < 0) throw new Error(`Cannot encode unknown UMI-family disposition ${family.disposition}.`);
    writer.u16(family.sampleIndex); writer.string(family.disposition === "BPB-rejects" ? "REJECTS" : family.umi);
    writer.string(family.mostLikelyParent); writer.u32(family.familySize); writer.f64(family.posteriorProbability); writer.u8(disposition);
  }
  return writer.finish();
}

export function decodeConsensusOutput(bytes: Uint8Array, config: PipelineConfig): { consensuses: ConsensusRecord[]; heteroduplexes: string[] } {
  const reader = new BinaryReader(bytes); reader.magic("WPO1"); const count = reader.u32(); const consensuses: ConsensusRecord[] = [];
  for (let index = 0; index < count; index++) {
    const sampleIndex = reader.u16(); reader.string(); // Legacy core-rendered label; normalize the public ID below.
    const umi = reader.string(), familySize = reader.u32();
    const minimumAgreement = reader.f64(), sequence = reader.string(), lowCount = reader.u32();
    const lowAgreementSites = Array.from({ length: lowCount }, () => ({ position: reader.u32(), agreement: reader.f32(),
      modalReadBase: String.fromCharCode(reader.u8()), modalRunLength: reader.u32() }));
    const sample = config.samples[sampleIndex]?.name ?? String(sampleIndex);
    consensuses.push({ id: `${sample}_${umi}`, sample, sampleIndex, umi,
      familySize, minimumAgreement, sequence, lowAgreementSites });
  }
  const heteroduplexCount = reader.u32(), heteroduplexes: string[] = [];
  for (let index = 0; index < heteroduplexCount; index++) heteroduplexes.push(`${reader.u16()}\0${reader.string()}`);
  if (!reader.done) throw new Error("Consensus payload has trailing bytes."); return { consensuses, heteroduplexes };
}

export function mergeStats(stats: QualityStats[], samples: number): QualityStats {
  const result: QualityStats = { totalReads: 0, qualityReads: 0, badReads: 0, shortReads: 0, longReads: 0,
    primerRejects: 0, idRejects: 0, demultiplexedReads: 0, bpbRejects: 0, malformedRecords: 0,
    downsampledReads: 0, perSample: Array(samples).fill(0) };
  for (const part of stats) for (const key of ["totalReads", "qualityReads", "badReads", "shortReads", "longReads", "primerRejects",
    "idRejects", "demultiplexedReads", "bpbRejects", "malformedRecords"] as const) result[key] += part[key];
  for (const part of stats) for (let index = 0; index < part.perSample.length; index++) result.perSample[index] += part.perSample[index];
  return result;
}
