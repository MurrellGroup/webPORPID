/// <reference lib="webworker" />
import { WASI } from "@bjorn3/browser_wasi_shim";
import alivibeWasmUrl from "/alivibe-msa.wasm?url";

interface AlivibeExports extends WebAssembly.Exports {
  memory: WebAssembly.Memory;
  alivibe_msa_alloc(size: number): number;
  alivibe_msa_free(pointer: number): void;
  alivibe_msa_run(pointer: number, length: number, iterations: number): number;
  alivibe_msa_run_nucleotide(pointer: number, length: number, iterations: number): number;
  alivibe_msa_run_amino_acid(pointer: number, length: number, iterations: number): number;
  alivibe_msa_result_ptr(): number;
  alivibe_msa_result_len(): number;
  alivibe_msa_error_ptr(): number;
  alivibe_msa_error_len(): number;
}

type Request = { input: ArrayBuffer; iterations: number; scoringMode: "literal" | "nucleotide" | "amino-acid" };
const decoder = new TextDecoder();

async function align(request: Request) {
  const url = alivibeWasmUrl;
  let module: WebAssembly.Module;
  try { module = await WebAssembly.compileStreaming(fetch(url)); }
  catch { const response = await fetch(url); if (!response.ok) throw new Error("Alivibe MSA could not be loaded."); module = await WebAssembly.compile(await response.arrayBuffer()); }
  const wasi = new WASI([], [], []), instance = await WebAssembly.instantiate(module, { wasi_snapshot_preview1: wasi.wasiImport });
  wasi.initialize(instance as WebAssembly.Instance & { exports: { memory: WebAssembly.Memory; _initialize?: () => unknown } });
  const runtime = instance.exports as AlivibeExports, input = new Uint8Array(request.input), pointer = runtime.alivibe_msa_alloc(input.byteLength);
  if (!pointer && input.byteLength) throw new Error("Alivibe MSA ran out of WebAssembly memory.");
  try {
    new Uint8Array(runtime.memory.buffer, pointer, input.byteLength).set(input);
    const run = request.scoringMode === "nucleotide" ? runtime.alivibe_msa_run_nucleotide
      : request.scoringMode === "amino-acid" ? runtime.alivibe_msa_run_amino_acid : runtime.alivibe_msa_run;
    if (run(pointer, input.byteLength, request.iterations) < 0) {
      const message = decoder.decode(new Uint8Array(runtime.memory.buffer, runtime.alivibe_msa_error_ptr(), runtime.alivibe_msa_error_len()));
      throw new Error(message || "Alivibe MSA failed.");
    }
    const result = new Uint8Array(runtime.memory.buffer, runtime.alivibe_msa_result_ptr(), runtime.alivibe_msa_result_len()).slice();
    self.postMessage({ type: "result", result: result.buffer }, { transfer: [result.buffer] });
  } finally { runtime.alivibe_msa_free(pointer); }
}

self.addEventListener("message", (event: MessageEvent<Request>) => void align(event.data).catch((error) => {
  self.postMessage({ type: "error", message: error instanceof Error ? error.message : String(error) });
}));
