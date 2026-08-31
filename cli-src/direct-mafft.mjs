import { readFile } from "node:fs/promises";
import { pathToFileURL } from "node:url";
import { Worker as NodeWorker } from "node:worker_threads";

const FFT_NS_2_ARGUMENTS = [
  "-q", "0", "-E", "2", "-V", "-1.53", "-s", "0.0", "-W", "6", "-O",
  "-C", "0-0", "-D", "-b", "62", "-g", "0", "-f", "-1.53", "-Q", "100.0",
  "-h", "0", "-F", "-X", "0.1", "-x", "-1", "-i", "/input.fa",
];

const numericFasta = (sequences) => sequences.map((sequence, index) => `>${index}\n${sequence}\n`).join("");

function parseNumericAlignment(source, expected) {
  const rows = new Map(); let index = -1, sequence = "";
  const finish = () => { if (index >= 0) { if (rows.has(index)) throw new Error("MAFFT returned a duplicate sequence identifier."); rows.set(index, sequence.toUpperCase()); } };
  for (const raw of source.split(/\r?\n/)) {
    const line = raw.trim(); if (!line) continue;
    if (line.startsWith(">")) {
      finish(); index = Number(line.slice(1).trim()); sequence = "";
      if (!Number.isSafeInteger(index) || index < 0 || index >= expected.length) throw new Error("MAFFT returned an unknown sequence identifier.");
    } else { if (index < 0) throw new Error("MAFFT returned sequence data before its first header."); sequence += line.replace(/\s/g, ""); }
  }
  finish();
  if (rows.size !== expected.length) throw new Error("MAFFT returned the wrong number of sequences.");
  const aligned = expected.map((input, row) => {
    const output = rows.get(row); if (output.replaceAll("-", "") !== input.toUpperCase()) throw new Error(`MAFFT changed sequence ${row + 1}.`); return output;
  });
  const width = aligned[0]?.length ?? 0;
  if (!width || aligned.some((row) => row.length !== width)) throw new Error("MAFFT returned a non-rectangular alignment.");
  return aligned;
}

export function createDirectMafftRunner(javascriptPath, wasmPath) {
  const factoryPromise = import(pathToFileURL(javascriptPath).href).then((module) => module.default);
  const wasmPromise = readFile(wasmPath);
  return async (sequences, signal, _iterations = 0, scoringMode = "nucleotide", onProgress) => {
    if (signal?.aborted) throw new Error("Analysis cancelled.");
    if (scoringMode === "amino-acid") throw new Error("The reference-panel MAFFT runner expects nucleotide sequences.");
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
    runtime.FS.writeFile("/input.fa", new TextEncoder().encode(numericFasta(input)));
    try { const status = runtime.callMain([...FFT_NS_2_ARGUMENTS]); if (status) throw new Error(`MAFFT exited with status ${status}.`); }
    finally { try { runtime.FS.unlink("/input.fa"); } catch { /* best effort */ } }
    return parseNumericAlignment(stdout.join("\n"), input);
  };
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

/** Run each sample's batch MAFFT alignment on an isolated CPU worker. */
export function createMafftRunner(javascriptPath, wasmPath, size = 1, workerPath = new URL("../porpid-mafft-worker.mjs", import.meta.url)) {
  const count = Math.max(1, Math.floor(size));
  const clients = Array.from({ length: count }, () => new MafftWorkerClient(new NodeWorker(workerPath))); let cursor = 0;
  const run = async (sequences, signal, iterations = 0, scoringMode = "nucleotide", onProgress) => {
    if (signal?.aborted) throw new Error("Analysis cancelled.");
    return clients[cursor++ % clients.length].call({ javascriptPath, wasmPath, sequences: sequences.map(String), iterations, scoringMode }, onProgress);
  };
  run.close = async () => { await Promise.all(clients.map((client) => client.close())); };
  return run;
}
