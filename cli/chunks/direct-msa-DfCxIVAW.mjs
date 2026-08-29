import { t as WASI } from "./dist-DmN76WoC.mjs";
import { readFile } from "node:fs/promises";
import { Worker } from "node:worker_threads";
//#region src/alivibe-msa-codec.ts
const MAGIC = [
	65,
	77,
	83,
	65
];
const encoder = new TextEncoder();
const decoder = new TextDecoder("utf-8", { fatal: true });
function encodeAlivibeMsaSequences(sequences) {
	if (sequences.length > 1e4) throw new Error("Alivibe MSA accepts at most 10,000 rows.");
	const rows = sequences.map((sequence, row) => {
		if (!/^[\x21-\x7e]*$/.test(sequence)) throw new Error(`Alivibe MSA row ${row + 1} contains a non-ASCII character.`);
		return encoder.encode(sequence);
	});
	const byteLength = 8 + rows.reduce((total, row) => total + 4 + row.byteLength, 0);
	if (byteLength > 256 * 1024 * 1024) throw new Error("Alivibe MSA input is too large.");
	const buffer = new ArrayBuffer(byteLength), bytes = new Uint8Array(buffer), view = new DataView(buffer);
	bytes.set(MAGIC);
	view.setUint32(4, rows.length, true);
	let offset = 8;
	for (const row of rows) {
		view.setUint32(offset, row.byteLength, true);
		offset += 4;
		bytes.set(row, offset);
		offset += row.byteLength;
	}
	return buffer;
}
function decodeAlivibeMsaSequences(buffer) {
	const bytes = new Uint8Array(buffer);
	if (bytes.byteLength < 8 || MAGIC.some((value, index) => bytes[index] !== value)) throw new Error("Alivibe MSA returned an invalid result.");
	const view = new DataView(buffer), count = view.getUint32(4, true), rows = [];
	let offset = 8;
	if (count > 1e4) throw new Error("Alivibe MSA returned too many rows.");
	for (let row = 0; row < count; row++) {
		if (offset + 4 > bytes.byteLength) throw new Error("Alivibe MSA returned a truncated result.");
		const length = view.getUint32(offset, true);
		offset += 4;
		if (offset + length > bytes.byteLength) throw new Error("Alivibe MSA returned a truncated row.");
		rows.push(decoder.decode(bytes.subarray(offset, offset + length)));
		offset += length;
	}
	if (offset !== bytes.byteLength) throw new Error("Alivibe MSA returned trailing bytes.");
	return rows;
}
function assertAlivibeMsaResult(input, aligned) {
	if (aligned.length !== input.length) throw new Error("Alivibe MSA returned the wrong number of rows.");
	const width = aligned[0]?.length ?? 0;
	aligned.forEach((sequence) => {
		if (sequence.length !== width) throw new Error("Alivibe MSA returned a non-rectangular alignment.");
	});
}
//#endregion
//#region cli-src/direct-msa.mjs
function createDirectMsaRunner(wasmPath) {
	let runtimePromise;
	const initialize = async () => {
		const wasi = new WASI([], [], []), module = await WebAssembly.compile(await readFile(wasmPath));
		const instance = await WebAssembly.instantiate(module, { wasi_snapshot_preview1: wasi.wasiImport });
		wasi.initialize(instance);
		return instance.exports;
	};
	return async (sequences, signal, iterations = 3, scoringMode = "nucleotide") => {
		if (signal?.aborted) throw new Error("Analysis cancelled.");
		if (sequences.length < 2) return [...sequences];
		const inputSequences = sequences.map(String), input = new Uint8Array(encodeAlivibeMsaSequences(inputSequences));
		const runtime = await (runtimePromise ??= initialize()), pointer = runtime.alivibe_msa_alloc(input.byteLength);
		if (!pointer && input.byteLength) throw new Error("Alivibe MSA ran out of WebAssembly memory.");
		try {
			new Uint8Array(runtime.memory.buffer, pointer, input.byteLength).set(input);
			if ((scoringMode === "nucleotide" ? runtime.alivibe_msa_run_nucleotide : scoringMode === "amino-acid" ? runtime.alivibe_msa_run_amino_acid : runtime.alivibe_msa_run)(pointer, input.byteLength, iterations) < 0) {
				const message = new TextDecoder().decode(new Uint8Array(runtime.memory.buffer, runtime.alivibe_msa_error_ptr(), runtime.alivibe_msa_error_len()));
				throw new Error(message || "Alivibe MSA failed.");
			}
			const aligned = decodeAlivibeMsaSequences(new Uint8Array(runtime.memory.buffer, runtime.alivibe_msa_result_ptr(), runtime.alivibe_msa_result_len()).slice().buffer);
			assertAlivibeMsaResult(inputSequences, aligned);
			return aligned;
		} finally {
			runtime.alivibe_msa_free(pointer);
		}
	};
}
var MsaWorkerClient = class {
	constructor(worker) {
		this.worker = worker;
		this.pending = /* @__PURE__ */ new Map();
		this.nextId = 1;
		this.tail = Promise.resolve();
		worker.on("message", (message) => {
			const pending = this.pending.get(message.id);
			if (!pending) return;
			this.pending.delete(message.id);
			message.error ? pending.reject(new Error(message.error)) : pending.resolve(message.result);
		});
		worker.on("error", (cause) => {
			for (const pending of this.pending.values()) pending.reject(cause);
			this.pending.clear();
		});
	}
	call(message) {
		const task = () => new Promise((resolve, reject) => {
			const id = this.nextId++;
			this.pending.set(id, {
				resolve,
				reject
			});
			this.worker.postMessage({
				id,
				...message
			});
		});
		const result = this.tail.then(task, task);
		this.tail = result.catch(() => {});
		return result;
	}
	close() {
		return this.worker.terminate();
	}
};
/** Run independent sample MSAs on separate worker-thread WASM instances. */
function createMsaRunner(wasmPath, size = 1, workerPath = new URL("../porpid-msa-worker.mjs", import.meta.url)) {
	const count = Math.max(1, Math.floor(size));
	if (count === 1) {
		const direct = createDirectMsaRunner(wasmPath);
		direct.close = async () => {};
		return direct;
	}
	const clients = Array.from({ length: count }, () => new MsaWorkerClient(new Worker(workerPath)));
	let cursor = 0;
	const run = async (sequences, signal, iterations = 3, scoringMode = "nucleotide") => {
		if (signal?.aborted) throw new Error("Analysis cancelled.");
		return clients[cursor++ % clients.length].call({
			wasmPath,
			sequences: sequences.map(String),
			iterations,
			scoringMode
		});
	};
	run.close = async () => {
		await Promise.all(clients.map((client) => client.close()));
	};
	return run;
}
//#endregion
export { encodeAlivibeMsaSequences as a, decodeAlivibeMsaSequences as i, createMsaRunner as n, assertAlivibeMsaResult as r, createDirectMsaRunner as t };
