import { readFile } from "node:fs/promises";
import { parentPort } from "node:worker_threads";
import { WebPorpidRuntime } from "../src/wasm-runtime.ts";

let runtime;
async function handle(message) {
  if (message.type === "init") {
    const module = await WebAssembly.compile(await readFile(message.wasmPath));
    runtime = await WebPorpidRuntime.create(module, new Uint8Array(message.config)); return { ready: true };
  }
  if (!runtime) throw new Error("The porpid-cli worker is not initialized.");
  if (message.type === "preprocess") return runtime.preprocess(message.text, message.firstOrdinal);
  if (message.type === "countFamilies") return runtime.countFamilies(new Uint8Array(message.bytes), new Uint8Array(message.cutoffs));
  if (message.type === "buildModel") return runtime.buildFamilyModel(new Uint8Array(message.bytes));
  if (message.type === "initModel") { runtime.initFamilyModel(new Uint8Array(message.bytes)); return { ready: true }; }
  if (message.type === "consensus") {
    if (message.model) runtime.initFamilyModel(new Uint8Array(message.model));
    return runtime.consensus(new Uint8Array(message.bytes), new Uint8Array(message.cutoffs));
  }
  if (message.type === "stats") return runtime.stats();
  throw new Error(`Unknown porpid-cli worker request ${message.type}.`);
}

const send = (message, transfer = []) => parentPort ? parentPort.postMessage(message, transfer) : globalThis.postMessage(message, transfer);
const receive = (message) => void handle(message).then((result) => {
  if (result instanceof Uint8Array) {
    const buffer = result.buffer.slice(result.byteOffset, result.byteOffset + result.byteLength); send({ id: message.id, result: buffer }, [buffer]);
  } else send({ id: message.id, result });
}).catch((cause) => send({ id: message.id, error: cause instanceof Error ? (cause.stack || cause.message) : String(cause) }));

if (parentPort) parentPort.on("message", receive);
else globalThis.addEventListener("message", (event) => receive(event.data));
