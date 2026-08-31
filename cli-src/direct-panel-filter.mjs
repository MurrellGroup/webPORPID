import { Worker as NodeWorker } from "node:worker_threads";
import { filterQueriesAgainstPanel } from "../src/independent-panel-filter.ts";

export function createIndependentPanelFilterRunner(workerPath = new URL("../porpid-panel-worker.mjs", import.meta.url)) {
  return async (sequences, panelRows, signal, workers = 1, onProgress) => {
    if (!sequences.length) return { sequences: [], scores: [] };
    if (signal?.aborted) throw new Error("Analysis cancelled.");
    const count = Math.max(1, Math.min(sequences.length, Math.floor(workers) || 1));
    if (count === 1) return filterQueriesAgainstPanel(sequences, panelRows,
      ({ completed, total }) => onProgress?.({ completed, total }));
    const outputSequences = new Array(sequences.length), scores = new Array(sequences.length), completed = new Uint32Array(count), instances = [];
    let rejectCancellation, settled = false;
    const cancellation = new Promise((_resolve, reject) => { rejectCancellation = reject; });
    const terminate = () => { for (const worker of instances) worker.terminate(); };
    const abort = () => { if (!settled) { terminate(); rejectCancellation(new Error("Analysis cancelled.")); } };
    signal?.addEventListener("abort", abort, { once: true });
    try {
      const jobs = Array.from({ length: count }, (_, workerIndex) => {
        const start = Math.floor(workerIndex * sequences.length / count), end = Math.floor((workerIndex + 1) * sequences.length / count);
        const worker = new NodeWorker(workerPath); instances.push(worker);
        return new Promise((resolve, reject) => {
          worker.on("message", (message) => {
            if (message.progress) {
              completed[workerIndex] = message.progress.completed;
              onProgress?.({ completed: completed.reduce((sum, value) => sum + value, 0), total: sequences.length }); return;
            }
            if (message.error) { reject(new Error(message.error)); return; }
            const result = message.result;
            if (!result || result.sequences.length !== end - start || result.scores.length !== end - start) {
              reject(new Error("An independent panel-filter worker returned an invalid result.")); return;
            }
            for (let index = 0; index < result.sequences.length; index++) {
              outputSequences[start + index] = result.sequences[index]; scores[start + index] = result.scores[index];
            }
            completed[workerIndex] = end - start; resolve();
          });
          worker.on("error", reject);
          worker.postMessage({ sequences: sequences.slice(start, end).map(String), panelRows: panelRows.map(String) });
        });
      });
      await Promise.race([Promise.all(jobs), cancellation]);
      if (signal?.aborted) throw new Error("Analysis cancelled."); settled = true;
      return { sequences: outputSequences, scores };
    } finally { signal?.removeEventListener("abort", abort); terminate(); }
  };
}
