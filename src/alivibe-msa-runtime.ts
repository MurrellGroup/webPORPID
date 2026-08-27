import { assertAlivibeMsaResult, decodeAlivibeMsaSequences, encodeAlivibeMsaSequences } from "./alivibe-msa-codec";

export async function runAlivibeMsa(
  sequences: readonly string[],
  signal?: AbortSignal,
  iterations = 3,
  scoringMode: "literal" | "nucleotide" | "amino-acid" = "nucleotide",
): Promise<string[]> {
  if (sequences.length < 2) return [...sequences];
  const inputSequences = sequences.map(String), input = encodeAlivibeMsaSequences(inputSequences);
  const worker = new Worker(new URL("./alivibe-msa-worker.ts", import.meta.url), { type: "module" });
  return new Promise<string[]>((resolve, reject) => {
    const abort = () => { worker.terminate(); reject(new DOMException("MSA alignment cancelled.", "AbortError")); };
    if (signal?.aborted) return abort(); signal?.addEventListener("abort", abort, { once: true });
    const finish = () => { signal?.removeEventListener("abort", abort); worker.terminate(); };
    worker.onmessage = (event: MessageEvent<{ type: "result" | "error"; result?: ArrayBuffer; message?: string }>) => {
      if (event.data.type === "error") { finish(); reject(new Error(event.data.message || "Alivibe MSA failed.")); return; }
      try { const aligned = decodeAlivibeMsaSequences(event.data.result ?? new ArrayBuffer(0)); assertAlivibeMsaResult(inputSequences, aligned); finish(); resolve(aligned); }
      catch (error) { finish(); reject(error); }
    };
    worker.onerror = (event) => { finish(); reject(new Error(event.message || "The Alivibe MSA worker stopped unexpectedly.")); };
    worker.postMessage({ input, iterations, scoringMode }, [input]);
  });
}
