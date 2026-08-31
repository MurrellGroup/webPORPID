export interface MafftProgress { detail: string }

/** Run MAFFT 7.520 FFT-NS-2 in an isolated worker with growable WASM memory. */
export async function runMafftFftnsMsa(
  sequences: readonly string[], signal?: AbortSignal, _iterations = 0,
  scoringMode: "literal" | "nucleotide" | "amino-acid" = "nucleotide",
  onProgress?: (progress: MafftProgress) => void,
): Promise<string[]> {
  if (sequences.length < 2) return [...sequences];
  if (scoringMode === "amino-acid") throw new Error("The reference-panel MAFFT runner expects nucleotide sequences.");
  const worker = new Worker(new URL("./mafft-msa-worker.ts", import.meta.url), { type: "module" });
  return new Promise<string[]>((resolve, reject) => {
    let settled = false;
    const finish = () => {
      if (settled) return false; settled = true;
      signal?.removeEventListener("abort", abort); worker.terminate(); return true;
    };
    const abort = () => { if (finish()) reject(new DOMException("MAFFT panel alignment cancelled.", "AbortError")); };
    worker.onmessage = (event: MessageEvent<{ type: "result" | "progress" | "error"; result?: string[]; detail?: string; message?: string }>) => {
      if (event.data.type === "progress") { onProgress?.({ detail: event.data.detail ?? "MAFFT is aligning candidate sequences" }); return; }
      if (!finish()) return;
      if (event.data.type === "error") reject(new Error(event.data.message || "MAFFT failed."));
      else resolve(event.data.result ?? []);
    };
    worker.onerror = (event) => { if (finish()) reject(new Error(event.message || "The MAFFT worker stopped unexpectedly.")); };
    if (signal?.aborted) { abort(); return; }
    signal?.addEventListener("abort", abort, { once: true });
    worker.postMessage({ sequences: sequences.map(String) });
  });
}
