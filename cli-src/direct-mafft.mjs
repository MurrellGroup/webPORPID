import { readFile } from "node:fs/promises";
import { pathToFileURL } from "node:url";
import { Worker as NodeWorker } from "node:worker_threads";

import { mafftArguments, numericMafftFasta, parseMafftAlignment, mafftProfileAddition } from "../src/mafft-msa-codec.ts";

export function createDirectMafftRunner(javascriptPath, wasmPath) {
  const factoryPromise = import(pathToFileURL(javascriptPath).href).then((module) => module.default);
  const wasmPromise = readFile(wasmPath);
  const run = async (sequences, signal, _iterations = 0, scoringMode = "nucleotide", onProgress, profileRows) => {
    if (signal?.aborted) throw new Error("Analysis cancelled.");
    if (sequences.length < 2) return [...sequences];
    const input = sequences.map((sequence) => String(sequence).toUpperCase()), stdout = [];
    if (input.some((sequence) => !sequence.length || /[^A-Z?*.-]/.test(sequence))) throw new Error("MAFFT input contains an unsupported symbol.");
    const [factory, wasmBinary] = await Promise.all([factoryPromise, wasmPromise]); let lastProgress = 0;
    const runtime = await factory({ wasmBinary, noInitialRun: true, print: (line) => stdout.push(String(line)),
      printErr: (line) => {
        const detail = String(line).replace(/[\b\r]+/g, " ").replace(/\s+/g, " ").trim(), now = performance.now();
        if (detail && (now - lastProgress >= 200 || /done|Progressive alignment|distance matrix/i.test(detail))) {
          lastProgress = now; onProgress?.({ detail });
        }
      } });
    runtime.FS.writeFile("/input.fa", new TextEncoder().encode(numericMafftFasta(input, scoringMode)));
    try { const status = runtime.callMain(mafftArguments(scoringMode, profileRows)); if (status) throw new Error(`MAFFT exited with status ${status}.`); }
    finally { try { runtime.FS.unlink("/input.fa"); } catch { /* best effort */ } }
    return parseMafftAlignment(stdout.join("\n"), input, scoringMode);
  };
  run.addToProfile = async (profileRows, sequence, signal) => mafftProfileAddition(profileRows,
    await run([...profileRows, sequence], signal, 0, "amino-acid", undefined, profileRows.length));
  return run;
}

class MafftWorkerClient {
  constructor(worker) {
    this.worker = worker; this.pending = new Map(); this.nextId = 1; this.tail = Promise.resolve();
    worker.on("message", (message) => {
      const pending = this.pending.get(message.id); if (!pending) return;
      if (message.progress) { pending.onProgress?.({ detail: message.progress }); return; }
      this.pending.delete(message.id); message.error ? pending.reject(new Error(message.error)) : pending.resolve(message.result);
    });
    worker.on("error", (cause) => { for (const pending of this.pending.values()) pending.reject(cause); this.pending.clear(); });
  }
  call(message, onProgress) {
    const task = () => new Promise((resolve, reject) => { const id = this.nextId++; this.pending.set(id, { resolve, reject, onProgress }); this.worker.postMessage({ id, ...message }); });
    const result = this.tail.then(task, task); this.tail = result.catch(() => {}); return result;
  }
  close() { return this.worker.terminate(); }
}

/** Run each sample's MAFFT alignment on an isolated CPU worker. */
export function createMafftRunner(javascriptPath, wasmPath, size = 1, workerPath = new URL("../porpid-mafft-worker.mjs", import.meta.url)) {
  const count = Math.max(1, Math.floor(size));
  const clients = Array.from({ length: count }, () => new MafftWorkerClient(new NodeWorker(workerPath))); let cursor = 0;
  const run = async (sequences, signal, iterations = 0, scoringMode = "nucleotide", onProgress) => {
    if (signal?.aborted) throw new Error("Analysis cancelled.");
    return clients[cursor++ % clients.length].call({ javascriptPath, wasmPath, sequences: sequences.map(String), iterations, scoringMode }, onProgress);
  };
  run.addToProfile = async (profileRows, sequence, signal) => {
    if (signal?.aborted) throw new Error("Analysis cancelled.");
    return clients[cursor++ % clients.length].call({ javascriptPath, wasmPath, profileRows: profileRows.map(String), sequence });
  };
  run.close = async () => { await Promise.all(clients.map((client) => client.close())); };
  return run;
}
