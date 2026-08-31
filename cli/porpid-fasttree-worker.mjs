import { t as createDirectFastTreeRunner } from "./chunks/direct-fasttree-8Kw3hjMk.mjs";
import { parentPort } from "node:worker_threads";
//#region cli-src/porpid-fasttree-worker.mjs
let javascriptPath, wasmPath, runner;
parentPort.on("message", async (message) => {
	try {
		if (!runner || javascriptPath !== message.javascriptPath || wasmPath !== message.wasmPath) {
			javascriptPath = message.javascriptPath;
			wasmPath = message.wasmPath;
			runner = createDirectFastTreeRunner(javascriptPath, wasmPath);
		}
		parentPort.postMessage({
			id: message.id,
			result: await runner(message.alignment)
		});
	} catch (cause) {
		parentPort.postMessage({
			id: message.id,
			error: cause instanceof Error ? cause.message : String(cause)
		});
	}
});
//#endregion
