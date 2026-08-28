import { parentPort } from "node:worker_threads";
import { createDirectFastTreeRunner } from "./direct-fasttree.mjs";

let javascriptPath, wasmPath, runner;
parentPort.on("message", async (message) => {
  try {
    if (!runner || javascriptPath !== message.javascriptPath || wasmPath !== message.wasmPath) {
      javascriptPath = message.javascriptPath; wasmPath = message.wasmPath;
      runner = createDirectFastTreeRunner(javascriptPath, wasmPath);
    }
    parentPort.postMessage({ id: message.id, result: await runner(message.alignment) });
  } catch (cause) { parentPort.postMessage({ id: message.id, error: cause instanceof Error ? cause.message : String(cause) }); }
});
