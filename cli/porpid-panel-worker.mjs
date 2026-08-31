import { t as filterQueriesAgainstPanel } from "./chunks/independent-panel-filter-DmDC1Ts9.mjs";
import { parentPort } from "node:worker_threads";
//#region cli-src/porpid-panel-worker.mjs
parentPort.on("message", (message) => {
	try {
		const result = filterQueriesAgainstPanel(message.sequences, message.panelRows, (progress) => parentPort.postMessage({ progress }));
		parentPort.postMessage({ result });
	} catch (cause) {
		parentPort.postMessage({ error: cause instanceof Error ? cause.message : String(cause) });
	}
});
//#endregion
