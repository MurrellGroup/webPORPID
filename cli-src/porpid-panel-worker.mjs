import { parentPort } from "node:worker_threads";
import { filterQueriesAgainstPanel } from "../src/independent-panel-filter.ts";

parentPort.on("message", (message) => {
  try {
    const result = filterQueriesAgainstPanel(message.sequences, message.panelRows,
      (progress) => parentPort.postMessage({ progress }));
    parentPort.postMessage({ result });
  } catch (cause) { parentPort.postMessage({ error: cause instanceof Error ? cause.message : String(cause) }); }
});
