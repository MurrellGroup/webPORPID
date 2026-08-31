import type { PanelFilterResult } from "./independent-panel-filter.ts";

export interface PanelFilterProgress { completed: number; total: number }

/** Split independent query/profile alignments over isolated CPU workers. */
export async function runIndependentPanelFilter(
  sequences: readonly string[], panelRows: readonly string[], signal?: AbortSignal, workers = 1,
  onProgress?: (progress: PanelFilterProgress) => void,
): Promise<PanelFilterResult> {
  if (!sequences.length) return { sequences: [], scores: [] };
  if (signal?.aborted) throw new DOMException("Independent panel filtering cancelled.", "AbortError");
  const count = Math.max(1, Math.min(sequences.length, Math.floor(workers) || 1));
  const outputSequences = new Array<string>(sequences.length), scores = new Array<number>(sequences.length);
  const instances: Worker[] = [], completed = new Uint32Array(count);
  let settled = false, rejectCancellation: ((cause: DOMException) => void) | undefined;
  const terminate = () => { for (const worker of instances) worker.terminate(); };
  const cancellation = new Promise<never>((_resolve, reject) => { rejectCancellation = reject; });
  const abort = () => {
    if (!settled) { terminate(); rejectCancellation?.(new DOMException("Independent panel filtering cancelled.", "AbortError")); }
  };
  signal?.addEventListener("abort", abort, { once: true });
  try {
    const jobs = Array.from({ length: count }, (_, workerIndex) => {
      const start = Math.floor(workerIndex * sequences.length / count), end = Math.floor((workerIndex + 1) * sequences.length / count);
      const worker = new Worker(new URL("./independent-panel-filter-worker.ts", import.meta.url), { type: "module" }); instances.push(worker);
      return new Promise<void>((resolve, reject) => {
        worker.onmessage = (event: MessageEvent<{ type: "result" | "progress" | "error"; result?: PanelFilterResult;
          completed?: number; message?: string }>) => {
          if (event.data.type === "progress") {
            completed[workerIndex] = event.data.completed ?? 0;
            onProgress?.({ completed: completed.reduce((sum, value) => sum + value, 0), total: sequences.length }); return;
          }
          if (event.data.type === "error") { reject(new Error(event.data.message || "Independent panel filtering failed.")); return; }
          const result = event.data.result;
          if (!result || result.sequences.length !== end - start || result.scores.length !== end - start) {
            reject(new Error("An independent panel-filter worker returned an invalid result.")); return;
          }
          for (let index = 0; index < result.sequences.length; index++) {
            outputSequences[start + index] = result.sequences[index]; scores[start + index] = result.scores[index];
          }
          completed[workerIndex] = end - start;
          onProgress?.({ completed: completed.reduce((sum, value) => sum + value, 0), total: sequences.length }); resolve();
        };
        worker.onerror = (event) => reject(new Error(event.message || "An independent panel-filter worker stopped unexpectedly."));
        worker.postMessage({ sequences: sequences.slice(start, end).map(String), panelRows: panelRows.map(String) });
      });
    });
    await Promise.race([Promise.all(jobs), cancellation]);
    if (signal?.aborted) throw new DOMException("Independent panel filtering cancelled.", "AbortError");
    settled = true; return { sequences: outputSequences, scores };
  } finally {
    signal?.removeEventListener("abort", abort); terminate();
  }
}
