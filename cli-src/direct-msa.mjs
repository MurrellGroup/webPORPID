import { readFile } from "node:fs/promises";
import { WASI } from "@bjorn3/browser_wasi_shim";
import { assertAlivibeMsaResult, decodeAlivibeMsaSequences, encodeAlivibeMsaSequences } from "../src/alivibe-msa-codec.ts";

export function createMsaRunner(wasmPath) {
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
