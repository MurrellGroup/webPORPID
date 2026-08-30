import { t as createDirectMsaRunner } from "./chunks/direct-msa-DNMUfcBa.mjs";
import { parentPort } from "node:worker_threads";
//#region cli-src/porpid-msa-worker.mjs
let wasmPath, runner;
parentPort.on("message", async (message) => {
	try {
		if (!runner || wasmPath !== message.wasmPath) {
			wasmPath = message.wasmPath;
			runner = createDirectMsaRunner(wasmPath);
		}
		const result = await runner(message.sequences, void 0, message.iterations, message.scoringMode);
		parentPort.postMessage({
			id: message.id,
			result
		});
	} catch (cause) {
		parentPort.postMessage({
			id: message.id,
			error: cause instanceof Error ? cause.message : String(cause)
		});
	}
});
//#endregion
