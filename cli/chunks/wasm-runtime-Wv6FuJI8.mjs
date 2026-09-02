import { n as BinaryWriter, t as BinaryReader } from "./binary-68r8u1WT.mjs";
import { t as WASI } from "./dist-BFrMSSwW.mjs";
//#region src/wasm-runtime.ts
const encoder = new TextEncoder();
const decoder = new TextDecoder();
var WebPorpidRuntime = class WebPorpidRuntime {
	core;
	constructor(core) {
		this.core = core;
	}
	static async create(module, compiledConfig) {
		const wasi = new WASI([], [], []);
		const instance = await WebAssembly.instantiate(module, { wasi_snapshot_preview1: wasi.wasiImport });
		wasi.initialize(instance);
		const runtime = new WebPorpidRuntime(instance.exports);
		runtime.withOne(compiledConfig, (pointer, length) => runtime.core.wpp_init_config(pointer, length), false);
		return runtime;
	}
	errorText() {
		return decoder.decode(new Uint8Array(this.core.memory.buffer, this.core.wpp_error_ptr(), this.core.wpp_error_len())) || "The webPORPID WASM core failed.";
	}
	resultBytes() {
		return new Uint8Array(this.core.memory.buffer, this.core.wpp_result_ptr(), this.core.wpp_result_len()).slice();
	}
	put(bytes) {
		const pointer = this.core.wpp_alloc(bytes.byteLength);
		if (!pointer && bytes.byteLength) throw new Error("The webPORPID WASM core ran out of memory.");
		new Uint8Array(this.core.memory.buffer, pointer, bytes.byteLength).set(bytes);
		return pointer;
	}
	withOne(bytes, call, result = true) {
		const pointer = this.put(bytes);
		try {
			if (call(pointer, bytes.byteLength) < 0) throw new Error(this.errorText());
			return result ? this.resultBytes() : new Uint8Array();
		} finally {
			this.core.wpp_free(pointer);
		}
	}
	withTwo(first, second, call) {
		const a = this.put(first), b = this.put(second);
		try {
			if (call(a, first.byteLength, b, second.byteLength) < 0) throw new Error(this.errorText());
			return this.resultBytes();
		} finally {
			this.core.wpp_free(a);
			this.core.wpp_free(b);
		}
	}
	preprocess(fastq, firstOrdinal) {
		return this.withOne(encoder.encode(fastq), (pointer, length) => this.core.wpp_preprocess(pointer, length, firstOrdinal));
	}
	partitionCounts(partition) {
		return this.withOne(partition, (pointer, length) => this.core.wpp_partition_counts(pointer, length));
	}
	countFamilies(partition, cutoffs) {
		return this.withTwo(partition, cutoffs, (a, al, b, bl) => this.core.wpp_count_families(a, al, b, bl));
	}
	buildFamilyModel(counts) {
		return this.withOne(counts, (pointer, length) => this.core.wpp_build_family_model(pointer, length));
	}
	initFamilyModel(model) {
		this.withOne(model, (pointer, length) => this.core.wpp_init_family_model(pointer, length), false);
	}
	consensus(partition, cutoffs) {
		return this.withTwo(partition, cutoffs, (a, al, b, bl) => this.core.wpp_consensus_partition(a, al, b, bl));
	}
	stats() {
		if (this.core.wpp_stats() < 0) throw new Error(this.errorText());
		return {
			...JSON.parse(decoder.decode(this.resultBytes())),
			downsampledReads: 0
		};
	}
};
function makeCutoffValues(sampleCounts, maximum) {
	const maxHash = (1n << 64n) - 1n;
	return sampleCounts.map((count) => {
		const cap = maximum < 1 ? count : BigInt(maximum);
		return count === 0n || count <= cap ? maxHash : maxHash * cap / count;
	});
}
function makeCutoffs(sampleCounts, maximum) {
	const writer = new BinaryWriter(), values = makeCutoffValues(sampleCounts, maximum);
	writer.magic("WPT1");
	writer.u32(values.length);
	for (const value of values) writer.u64(value);
	return writer.finish();
}
function decodeFamilyCounts(bytes) {
	const reader = new BinaryReader(bytes);
	reader.magic("WPN1");
	const count = reader.u32();
	const output = [];
	for (let index = 0; index < count; index++) output.push({
		sample: reader.u16(),
		umi: reader.string(),
		count: reader.u32()
	});
	if (!reader.done) throw new Error("Family count payload has trailing bytes.");
	return output;
}
function mergeFamilyCounts(parts) {
	const counts = /* @__PURE__ */ new Map();
	for (const bytes of parts) for (const entry of decodeFamilyCounts(bytes)) {
		const key = `${entry.sample}\0${entry.umi}`, previous = counts.get(key);
		if (previous) previous.count += entry.count;
		else counts.set(key, { ...entry });
	}
	const entries = [...counts.values()].sort((a, b) => a.sample - b.sample || a.umi.localeCompare(b.umi));
	const writer = new BinaryWriter();
	writer.magic("WPN1");
	writer.u32(entries.length);
	for (const entry of entries) {
		writer.u16(entry.sample);
		writer.string(entry.umi);
		writer.u32(entry.count);
	}
	return writer.finish();
}
const DISPOSITIONS = [
	"likely_real",
	"BPB-rejects",
	"heteroduplex",
	"LDA-rejects",
	"UMI_len != 8",
	"family-size-reject"
];
function decodeFamilyModel(bytes, config) {
	const reader = new BinaryReader(bytes);
	reader.magic("WPM1");
	const count = reader.u32();
	const output = [];
	for (let index = 0; index < count; index++) {
		const sampleIndex = reader.u16(), umi = reader.string(), parent = reader.string(), familySize = reader.u32();
		const probability = reader.f64(), disposition = DISPOSITIONS[reader.u8()];
		if (!disposition || !config.samples[sampleIndex]) throw new Error("Family model contains an invalid sample or disposition.");
		output.push({
			sample: config.samples[sampleIndex].name,
			sampleIndex,
			umi: disposition === "BPB-rejects" ? "REJECTED" : umi,
			familySize,
			mostLikelyParent: parent,
			posteriorProbability: probability,
			logOffspringProbability: Math.log(1 - probability),
			disposition
		});
	}
	if (!reader.done) throw new Error("Family model has trailing bytes.");
	return output;
}
function decodeConsensusOutput(bytes, config) {
	const reader = new BinaryReader(bytes);
	reader.magic("WPO1");
	const count = reader.u32();
	const consensuses = [];
	for (let index = 0; index < count; index++) {
		const sampleIndex = reader.u16();
		reader.string();
		const umi = reader.string(), familySize = reader.u32();
		const minimumAgreement = reader.f64(), sequence = reader.string(), lowCount = reader.u32();
		const lowAgreementSites = Array.from({ length: lowCount }, () => ({
			position: reader.u32(),
			agreement: reader.f32(),
			modalReadBase: String.fromCharCode(reader.u8()),
			modalRunLength: reader.u32()
		}));
		const sample = config.samples[sampleIndex]?.name ?? String(sampleIndex);
		consensuses.push({
			id: `${sample}_${umi}`,
			sample,
			sampleIndex,
			umi,
			familySize,
			minimumAgreement,
			sequence,
			lowAgreementSites
		});
	}
	const heteroduplexCount = reader.u32(), heteroduplexes = [];
	for (let index = 0; index < heteroduplexCount; index++) heteroduplexes.push(`${reader.u16()}\0${reader.string()}`);
	if (!reader.done) throw new Error("Consensus payload has trailing bytes.");
	return {
		consensuses,
		heteroduplexes
	};
}
function mergeStats(stats, samples) {
	const result = {
		totalReads: 0,
		qualityReads: 0,
		badReads: 0,
		shortReads: 0,
		longReads: 0,
		primerRejects: 0,
		idRejects: 0,
		demultiplexedReads: 0,
		bpbRejects: 0,
		malformedRecords: 0,
		downsampledReads: 0,
		perSample: Array(samples).fill(0)
	};
	for (const part of stats) for (const key of [
		"totalReads",
		"qualityReads",
		"badReads",
		"shortReads",
		"longReads",
		"primerRejects",
		"idRejects",
		"demultiplexedReads",
		"bpbRejects",
		"malformedRecords"
	]) result[key] += part[key];
	for (const part of stats) for (let index = 0; index < part.perSample.length; index++) result.perSample[index] += part.perSample[index];
	return result;
}
//#endregion
export { makeCutoffValues as a, mergeStats as c, decodeFamilyModel as i, decodeConsensusOutput as n, makeCutoffs as o, decodeFamilyCounts as r, mergeFamilyCounts as s, WebPorpidRuntime as t };
