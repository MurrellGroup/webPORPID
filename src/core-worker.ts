/// <reference lib="webworker" />

import { WebPorpidRuntime } from "./wasm-runtime";

type Request =
  | { id: number; type: "initialize"; module: WebAssembly.Module; config: ArrayBuffer }
  | { id: number; type: "shutdown" }
  | { id: number; type: "preprocess"; text: string; firstOrdinal: number }
  | { id: number; type: "partitionCounts"; bytes: ArrayBuffer }
  | { id: number; type: "countFamilies"; bytes: ArrayBuffer; cutoffs: ArrayBuffer }
  | { id: number; type: "buildModel"; bytes: ArrayBuffer }
  | { id: number; type: "initModel"; bytes: ArrayBuffer }
  | { id: number; type: "consensus"; bytes: ArrayBuffer; cutoffs: ArrayBuffer; model?: ArrayBuffer }
  | { id: number; type: "stats" };

let runtime: WebPorpidRuntime | null = null;

async function handle(request: Request): Promise<unknown> {
  if (request.type === "initialize") {
    runtime = await WebPorpidRuntime.create(request.module, new Uint8Array(request.config)); return { ready: true };
  }
  if (request.type === "shutdown") { runtime = null; return { closed: true }; }
  if (!runtime) throw new Error("The webPORPID compute worker is not initialized.");
  switch (request.type) {
    case "preprocess": return runtime.preprocess(request.text, request.firstOrdinal);
    case "partitionCounts": return runtime.partitionCounts(new Uint8Array(request.bytes));
    case "countFamilies": return runtime.countFamilies(new Uint8Array(request.bytes), new Uint8Array(request.cutoffs));
    case "buildModel": return runtime.buildFamilyModel(new Uint8Array(request.bytes));
    case "initModel": runtime.initFamilyModel(new Uint8Array(request.bytes)); return { ready: true };
    case "consensus":
      if (request.model) runtime.initFamilyModel(new Uint8Array(request.model));
      return runtime.consensus(new Uint8Array(request.bytes), new Uint8Array(request.cutoffs));
    case "stats": return runtime.stats();
  }
}

self.addEventListener("message", (event: MessageEvent<Request>) => {
  const request = event.data;
  Promise.resolve(handle(request)).then((result) => {
    if (result instanceof Uint8Array) {
      const buffer = result.buffer.slice(result.byteOffset, result.byteOffset + result.byteLength);
      self.postMessage({ id: request.id, result: buffer }, { transfer: [buffer] });
    } else self.postMessage({ id: request.id, result });
    if (request.type === "shutdown") self.close();
  }).catch((cause) => {
    const raw = cause instanceof Error ? cause.message : String(cause);
    const error = request.type === "consensus" && /unreachable|out of memory|memory access out of bounds/i.test(raw)
      ? `A WebAssembly consensus worker stopped inside its current read block (usually a per-block memory limit). Increase the number of temporary read blocks or lower CPU workers, then retry. Internal error: ${raw}`
      : raw;
    self.postMessage({ id: request.id, error });
  });
});
