/// <reference lib="webworker" />

import { WebPorpidRuntime } from "./wasm-runtime";

type Request =
  | { id: number; type: "initialize"; module: WebAssembly.Module; config: ArrayBuffer }
  | { id: number; type: "preprocess"; text: string; firstOrdinal: number }
  | { id: number; type: "partitionCounts"; bytes: ArrayBuffer }
  | { id: number; type: "countFamilies"; bytes: ArrayBuffer; cutoffs: ArrayBuffer }
  | { id: number; type: "buildModel"; bytes: ArrayBuffer }
  | { id: number; type: "initModel"; bytes: ArrayBuffer }
  | { id: number; type: "consensus"; bytes: ArrayBuffer; cutoffs: ArrayBuffer }
  | { id: number; type: "stats" };

let runtime: WebPorpidRuntime | null = null;

async function handle(request: Request): Promise<unknown> {
  if (request.type === "initialize") {
    runtime = await WebPorpidRuntime.create(request.module, new Uint8Array(request.config)); return { ready: true };
  }
  if (!runtime) throw new Error("The webPORPID compute worker is not initialized.");
  switch (request.type) {
    case "preprocess": return runtime.preprocess(request.text, request.firstOrdinal);
    case "partitionCounts": return runtime.partitionCounts(new Uint8Array(request.bytes));
    case "countFamilies": return runtime.countFamilies(new Uint8Array(request.bytes), new Uint8Array(request.cutoffs));
    case "buildModel": return runtime.buildFamilyModel(new Uint8Array(request.bytes));
    case "initModel": runtime.initFamilyModel(new Uint8Array(request.bytes)); return { ready: true };
    case "consensus": return runtime.consensus(new Uint8Array(request.bytes), new Uint8Array(request.cutoffs));
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
  }).catch((cause) => self.postMessage({ id: request.id, error: cause instanceof Error ? cause.message : String(cause) }));
});
