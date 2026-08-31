import { t as createDirectMafftRunner } from "./chunks/direct-mafft-IVivIuqW.mjs";
import { parentPort } from "node:worker_threads";
//#region cli-src/porpid-mafft-worker.mjs
let key = "", runner;
parentPort.on("message", async (message) => {
	try {
		const nextKey = `${message.javascriptPath}\0${message.wasmPath}`;
		if (!runner || key !== nextKey) {
			key = nextKey;
			runner = createDirectMafftRunner(message.javascriptPath, message.wasmPath);
		}
		const result = await runner(message.sequences, void 0, message.iterations, message.scoringMode, ({ detail }) => parentPort.postMessage({
			id: message.id,
			progress: detail
		}));
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
