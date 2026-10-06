import { readFile } from "node:fs/promises";
import { Worker } from "node:worker_threads";
import { pathToFileURL } from "node:url";
//#region src/mafft-msa-codec.ts
/** disttbfast flags used by MAFFT 7.520's FFT-NS-2 and --add wrappers. */
function mafftArguments(mode, profileRows) {
	return [
		"-q",
		"0",
		"-E",
		"2",
		"-V",
		"-1.53",
		"-s",
		"0.0",
		"-W",
		"6",
		"-O",
		"-C",
		"0-0",
		mode === "amino-acid" ? "-P" : "-D",
		"-b",
		"62",
		"-g",
		"0",
		"-f",
		"-1.53",
		"-Q",
		"100.0",
		"-h",
		"0",
		"-F",
		"-X",
		"0.1",
		"-x",
		"-1",
		...profileRows === void 0 ? [] : [
			"-K",
			"-I",
			"1"
		],
		"-i",
		"/input.fa"
	];
}
function numericMafftFasta(sequences, mode) {
	return sequences.map((sequence, index) => {
		return `>${index}\n${mode === "amino-acid" ? sequence.replaceAll("*", "X") : sequence}\n`;
	}).join("");
}
function parseMafftAlignment(source, expected, mode) {
	const rows = /* @__PURE__ */ new Map();
	let index = -1, sequence = "";
	const finish = () => {
		if (index < 0) return;
		if (rows.has(index)) throw new Error("MAFFT returned a duplicate sequence identifier.");
		rows.set(index, sequence.toUpperCase());
	};
	for (const raw of source.split(/\r?\n/)) {
		const line = raw.trim();
		if (!line) continue;
		if (line.startsWith(">")) {
			finish();
			index = Number(line.slice(1).trim());
			sequence = "";
			if (!Number.isSafeInteger(index) || index < 0 || index >= expected.length) throw new Error("MAFFT returned an unknown sequence identifier.");
		} else {
			if (index < 0) throw new Error("MAFFT returned sequence data before its first header.");
			sequence += line.replace(/\s/g, "");
		}
	}
	finish();
	if (rows.size !== expected.length) throw new Error("MAFFT returned the wrong number of sequences.");
	const aligned = expected.map((input, row) => {
		const original = input.replaceAll("-", "").toUpperCase(), output = rows.get(row);
		const normalized = mode === "amino-acid" ? original.replaceAll("*", "X") : original;
		if (output.replaceAll("-", "") !== normalized) throw new Error(`MAFFT changed sequence ${row + 1}.`);
		let cursor = 0;
		return [...output].map((residue) => residue === "-" ? "-" : original[cursor++]).join("");
	});
	const width = aligned[0]?.length ?? 0;
	if (!width || aligned.some((row) => row.length !== width)) throw new Error("MAFFT returned a non-rectangular alignment.");
	return aligned;
}
function mafftProfileAddition(profile, aligned) {
	if (aligned.length !== profile.length + 1) throw new Error("MAFFT profile addition returned the wrong number of rows.");
	const sampleRows = aligned.slice(0, profile.length);
	const withoutGapColumns = (rows) => {
		const columns = Array.from({ length: rows[0]?.length ?? 0 }, (_, column) => column).filter((column) => rows.some((row) => row[column] !== "-"));
		return rows.map((row) => columns.map((column) => row[column]).join(""));
	};
	const before = withoutGapColumns(profile), after = withoutGapColumns(sampleRows);
	if (before.some((row, index) => row !== after[index])) throw new Error("MAFFT reference addition changed the existing sample protein alignment.");
	return {
		profileRows: sampleRows,
		sequence: aligned[profile.length]
	};
}
//#endregion
//#region cli-src/direct-mafft.mjs
function createDirectMafftRunner(javascriptPath, wasmPath) {
	const factoryPromise = import(pathToFileURL(javascriptPath).href).then((module) => module.default);
	const wasmPromise = readFile(wasmPath);
	const run = async (sequences, signal, _iterations = 0, scoringMode = "nucleotide", onProgress, profileRows) => {
		if (signal?.aborted) throw new Error("Analysis cancelled.");
		if (sequences.length < 2) return [...sequences];
		const input = sequences.map((sequence) => String(sequence).toUpperCase()), stdout = [];
		if (input.some((sequence) => !sequence.length || /[^A-Z?*.-]/.test(sequence))) throw new Error("MAFFT input contains an unsupported symbol.");
		const [factory, wasmBinary] = await Promise.all([factoryPromise, wasmPromise]);
		let lastProgress = 0;
		const runtime = await factory({
			wasmBinary,
			noInitialRun: true,
			print: (line) => stdout.push(String(line)),
			printErr: (line) => {
				const detail = String(line).replace(/[\b\r]+/g, " ").replace(/\s+/g, " ").trim(), now = performance.now();
				if (detail && (now - lastProgress >= 200 || /done|Progressive alignment|distance matrix/i.test(detail))) {
					lastProgress = now;
					onProgress?.({ detail });
				}
			}
		});
		runtime.FS.writeFile("/input.fa", new TextEncoder().encode(numericMafftFasta(input, scoringMode)));
		try {
			const status = runtime.callMain(mafftArguments(scoringMode, profileRows));
			if (status) throw new Error(`MAFFT exited with status ${status}.`);
		} finally {
			try {
				runtime.FS.unlink("/input.fa");
			} catch {}
		}
		return parseMafftAlignment(stdout.join("\n"), input, scoringMode);
	};
	run.addToProfile = async (profileRows, sequence, signal) => mafftProfileAddition(profileRows, await run([...profileRows, sequence], signal, 0, "amino-acid", void 0, profileRows.length));
	return run;
}
var MafftWorkerClient = class {
	constructor(worker) {
		this.worker = worker;
		this.pending = /* @__PURE__ */ new Map();
		this.nextId = 1;
		this.tail = Promise.resolve();
		worker.on("message", (message) => {
			const pending = this.pending.get(message.id);
			if (!pending) return;
			if (message.progress) {
				pending.onProgress?.({ detail: message.progress });
				return;
			}
			this.pending.delete(message.id);
			message.error ? pending.reject(new Error(message.error)) : pending.resolve(message.result);
		});
		worker.on("error", (cause) => {
			for (const pending of this.pending.values()) pending.reject(cause);
			this.pending.clear();
		});
	}
	call(message, onProgress) {
		const task = () => new Promise((resolve, reject) => {
			const id = this.nextId++;
			this.pending.set(id, {
				resolve,
				reject,
				onProgress
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
/** Run each sample's MAFFT alignment on an isolated CPU worker. */
function createMafftRunner(javascriptPath, wasmPath, size = 1, workerPath = new URL("../porpid-mafft-worker.mjs", import.meta.url)) {
	const count = Math.max(1, Math.floor(size));
	const clients = Array.from({ length: count }, () => new MafftWorkerClient(new Worker(workerPath)));
	let cursor = 0;
	const run = async (sequences, signal, iterations = 0, scoringMode = "nucleotide", onProgress) => {
		if (signal?.aborted) throw new Error("Analysis cancelled.");
		return clients[cursor++ % clients.length].call({
			javascriptPath,
			wasmPath,
			sequences: sequences.map(String),
			iterations,
			scoringMode
		}, onProgress);
	};
	run.addToProfile = async (profileRows, sequence, signal) => {
		if (signal?.aborted) throw new Error("Analysis cancelled.");
		return clients[cursor++ % clients.length].call({
			javascriptPath,
			wasmPath,
			profileRows: profileRows.map(String),
			sequence
		});
	};
	run.close = async () => {
		await Promise.all(clients.map((client) => client.close()));
	};
	return run;
}
//#endregion
export { createMafftRunner as n, mafftProfileAddition as r, createDirectMafftRunner as t };
