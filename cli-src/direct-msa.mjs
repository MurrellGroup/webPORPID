import { readFile } from "node:fs/promises";
import { Worker as NodeWorker } from "node:worker_threads";
import { WASI } from "@bjorn3/browser_wasi_shim";
import { assertAlivibeMsaResult, decodeAlivibeMsaSequences, encodeAlivibeMsaSequences } from "../src/alivibe-msa-codec.ts";

export function createDirectMsaRunner(wasmPath) {
  let runtimePromise;
  const initialize = async () => {
    const wasi = new WASI([], [], []), module = await WebAssembly.compile(await readFile(wasmPath));
    const instance = await WebAssembly.instantiate(module, { wasi_snapshot_preview1: wasi.wasiImport });
    wasi.initialize(instance); return instance.exports;
  };
  return async (sequences, signal, iterations = 3, scoringMode = "nucleotide") => {
    if (signal?.aborted) throw new Error("Analysis cancelled.");
    if (sequences.length < 2) return [...sequences];
    const inputSequences = sequences.map(String), input = new Uint8Array(encodeAlivibeMsaSequences(inputSequences));
    const runtime = await (runtimePromise ??= initialize()), pointer = runtime.alivibe_msa_alloc(input.byteLength);
    if (!pointer && input.byteLength) throw new Error("Alivibe MSA ran out of WebAssembly memory.");
    try {
      new Uint8Array(runtime.memory.buffer, pointer, input.byteLength).set(input);
      const run = scoringMode === "nucleotide" ? runtime.alivibe_msa_run_nucleotide
        : scoringMode === "amino-acid" ? runtime.alivibe_msa_run_amino_acid : runtime.alivibe_msa_run;
      if (run(pointer, input.byteLength, iterations) < 0) {
        const message = new TextDecoder().decode(new Uint8Array(runtime.memory.buffer,
          runtime.alivibe_msa_error_ptr(), runtime.alivibe_msa_error_len()));
        throw new Error(message || "Alivibe MSA failed.");
      }
      const output = new Uint8Array(runtime.memory.buffer, runtime.alivibe_msa_result_ptr(), runtime.alivibe_msa_result_len()).slice();
      const aligned = decodeAlivibeMsaSequences(output.buffer); assertAlivibeMsaResult(inputSequences, aligned); return aligned;
    } finally { runtime.alivibe_msa_free(pointer); }
  };
}

class MsaWorkerClient {
  constructor(worker) {
    this.worker = worker; this.pending = new Map(); this.nextId = 1; this.tail = Promise.resolve();
    worker.on("message", (message) => { const pending = this.pending.get(message.id); if (!pending) return; this.pending.delete(message.id); message.error ? pending.reject(new Error(message.error)) : pending.resolve(message.result); });
    worker.on("error", (cause) => { for (const pending of this.pending.values()) pending.reject(cause); this.pending.clear(); });
  }
  call(message) {
    const task = () => new Promise((resolve, reject) => { const id = this.nextId++; this.pending.set(id, { resolve, reject }); this.worker.postMessage({ id, ...message }); });
    const result = this.tail.then(task, task); this.tail = result.catch(() => {}); return result;
  }
  close() { return this.worker.terminate(); }
}

/** Run independent sample MSAs on separate worker-thread WASM instances. */
export function createMsaRunner(wasmPath, size = 1, workerPath = new URL("../porpid-msa-worker.mjs", import.meta.url)) {
  const count = Math.max(1, Math.floor(size));
  if (count === 1) {
    const direct = createDirectMsaRunner(wasmPath); direct.close = async () => {}; return direct;
  }
  const clients = Array.from({ length: count }, () => new MsaWorkerClient(new NodeWorker(workerPath)));
  let cursor = 0;
  const run = async (sequences, signal, iterations = 3, scoringMode = "nucleotide") => {
    if (signal?.aborted) throw new Error("Analysis cancelled.");
    return clients[cursor++ % clients.length].call({ wasmPath, sequences: sequences.map(String), iterations, scoringMode });
  };
  run.close = async () => { await Promise.all(clients.map((client) => client.close())); };
  return run;
}
