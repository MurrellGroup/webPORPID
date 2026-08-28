import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import { WASI } from "node:wasi";

class Writer {
  chunks = [];
  raw(bytes) { this.chunks.push(bytes); }
  magic(text) { this.raw(new TextEncoder().encode(text)); }
  number(method, bytes, value) { const data = new Uint8Array(bytes); new DataView(data.buffer)[method](0, value, true); this.raw(data); }
  u8(value) { this.number("setUint8", 1, value); }
  u32(value) { this.number("setUint32", 4, value); }
  i32(value) { this.number("setInt32", 4, value); }
  u64(value) { this.number("setBigUint64", 8, BigInt(value)); }
  f64(value) { this.number("setFloat64", 8, value); }
  string(value) { const bytes = new TextEncoder().encode(value); this.u32(bytes.length); this.raw(bytes); }
  finish() { const size = this.chunks.reduce((sum, part) => sum + part.length, 0), output = new Uint8Array(size); let offset = 0; for (const part of this.chunks) { output.set(part, offset); offset += part.length; } return output; }
}

function configBytes() {
  const w = new Writer(); w.magic("WPC1"); w.u32(1); w.string("synthetic");
  w.f64(.05); w.u32(20); w.u32(300); w.u32(0); w.u32(150); w.u32(0); w.u32(100000); w.u32(1);
  w.f64(.995); w.f64(.015); w.f64(.2); w.f64(.015); w.u8(1); w.f64(.6); w.f64(.25); w.f64(.99); w.f64(50); w.f64(.7); w.u32(8); w.u64(0x504f52504944n);
  w.u32(1); w.string("sample_1"); w.string("CCGCTacgtaaNNNNNNNNGTCA"); w.string("TAGG"); w.string("panel.fa"); w.string(""); w.i32(-1);
  w.f64(Number.NaN); w.f64(Number.NaN); w.f64(Number.NaN); w.f64(Number.NaN); w.u32(0); w.string(""); return w.finish();
}

const complements = { A: "T", C: "G", G: "C", T: "A", N: "N" };
const reverseComplement = (text) => [...text].reverse().map((base) => complements[base] ?? "N").join("");

const bytes = await readFile(new URL("../public/webporpid.wasm", import.meta.url));
const module = await WebAssembly.compile(bytes);
const exported = new Set(WebAssembly.Module.exports(module).map((entry) => entry.name));
for (const name of ["memory", "wpp_alloc", "wpp_free", "wpp_version", "wpp_init_config", "wpp_preprocess", "wpp_consensus_partition"]) assert(exported.has(name), `missing WASM export ${name}`);

const wasi = new WASI({ version: "preview1" });
const instance = await WebAssembly.instantiate(module, { wasi_snapshot_preview1: wasi.wasiImport });
wasi.initialize(instance);
const core = instance.exports;
const versionPointer = core.wpp_version(), versionMemory = new Uint8Array(core.memory.buffer);
let versionEnd = versionPointer; while (versionMemory[versionEnd]) versionEnd++;
assert.equal(new TextDecoder().decode(versionMemory.subarray(versionPointer, versionEnd)), "0.2.0");
const put = (input) => { const pointer = core.wpp_alloc(input.length); new Uint8Array(core.memory.buffer, pointer, input.length).set(input); return pointer; };
const takeResult = () => new Uint8Array(core.memory.buffer, core.wpp_result_ptr(), core.wpp_result_len()).slice();
const one = (input, call) => { const pointer = put(input); try { assert(call(pointer, input.length) >= 0); return takeResult(); } finally { core.wpp_free(pointer); } };
const two = (first, second, call) => {
  const a = put(first), b = put(second); try { assert(call(a, first.length, b, second.length) >= 0); return takeResult(); }
  finally { core.wpp_free(a); core.wpp_free(b); }
};
const magic = (input, expected) => assert.equal(new TextDecoder().decode(input.subarray(0, expected.length)), expected);

const config = configBytes(), configPointer = put(config);
assert.equal(core.wpp_init_config(configPointer, config.length), 1); core.wpp_free(configPointer);
const payload = "ATGCCTTGGGCCATCGGACCATATGTTTACGATGGGCAGCTGACTACCGACAACCGTCAATTCGTCTCAGAGAAGTAA";
const payloads = [payload, payload, payload, payload.slice(0, 41) + "T" + payload.slice(42), payload.slice(0, 48) + "A" + payload.slice(48)];
const fastqText = payloads.map((value, index) => {
  const oriented = "ACGTAAAACCGGTTGTCA" + value, sequence = "TAGG" + reverseComplement(oriented) + reverseComplement("CCGCT");
  return `@synthetic_${index + 1}\n${sequence}\n+\n${"I".repeat(sequence.length)}\n`;
}).join("");
const routed = one(new TextEncoder().encode(fastqText), (pointer, length) => core.wpp_preprocess(pointer, length, 10));
assert(routed.length > 5, "WASM preprocessing returned no routed spool record");
assert(core.wpp_stats() > 0); const stats = JSON.parse(new TextDecoder().decode(takeResult()));
assert.deepEqual({ total: stats.totalReads, quality: stats.qualityReads, demux: stats.demultiplexedReads, bpb: stats.bpbRejects }, { total: 5, quality: 5, demux: 5, bpb: 0 });

const records = []; let offset = 0, routedView = new DataView(routed.buffer, routed.byteOffset, routed.byteLength);
while (offset < routed.length) {
  offset++; const length = routedView.getUint32(offset, true); offset += 4; records.push(routed.slice(offset, offset + length)); offset += length;
}
const partitionLength = records.reduce((sum, record) => sum + record.length, 0), partition = new Uint8Array(partitionLength);
offset = 0; for (const record of records) { partition.set(record, offset); offset += record.length; }
const cutoffWriter = new Writer(); cutoffWriter.magic("WPT1"); cutoffWriter.u32(1); cutoffWriter.u64((1n << 64n) - 1n); const cutoffs = cutoffWriter.finish();
const sampleCounts = one(partition, (pointer, length) => core.wpp_partition_counts(pointer, length)); magic(sampleCounts, "WPS1");
const familyCounts = two(partition, cutoffs, (a, al, b, bl) => core.wpp_count_families(a, al, b, bl)); magic(familyCounts, "WPN1");
const model = one(familyCounts, (pointer, length) => core.wpp_build_family_model(pointer, length)); magic(model, "WPM1");
one(model, (pointer, length) => core.wpp_init_family_model(pointer, length));
const consensus = two(partition, cutoffs, (a, al, b, bl) => core.wpp_consensus_partition(a, al, b, bl)); magic(consensus, "WPO1");
assert.equal(new DataView(consensus.buffer, consensus.byteOffset + 4).getUint32(0, true), 1, "WASM consensus pipeline returned one family");
console.log("webPORPID full WASI reactor pipeline test passed");
