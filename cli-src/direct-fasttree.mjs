import { readFile } from "node:fs/promises";
import { Worker as NodeWorker } from "node:worker_threads";
import { parseFasta } from "../src/config.ts";
import { treeTipNames } from "../src/tree-names.ts";

function completeNewick(output) {
  const candidates = output.split(/\r?\n/).filter((line) => line.includes("(") && line.includes(";")).reverse();
  for (const candidate of candidates) {
    const start = candidate.indexOf("("), end = candidate.lastIndexOf(";");
    if (start >= 0 && end > start) return candidate.slice(start, end + 1);
  }
  throw new Error("FastTree did not return a complete Newick tree.");
}

export function createDirectFastTreeRunner(javascriptPath, wasmPath) {
  let output = [], errors = [], invocation = 0;
  const modulePromise = Promise.all([readFile(javascriptPath, "utf8"), readFile(wasmPath)]).then(([source, wasmBinary]) => {
    const commonJs = { exports: {} };
    Function("module", "exports", source)(commonJs, commonJs.exports);
    if (typeof commonJs.exports !== "function") throw new Error("FastTree JavaScript runtime did not export a module factory.");
    return commonJs.exports({ wasmBinary,
      print: (line) => output.push(String(line)), printErr: (line) => errors.push(String(line)) });
  });
  return async (alignedFasta) => {
    const records = parseFasta(alignedFasta);
    const safeNames = treeTipNames(records.map((record) => record.name));
    if (records.length < 3) return records.length === 2
      ? `(${safeNames[0]}:0.0,${safeNames[1]}:0.0);`
      : records.length === 1 ? `(${safeNames[0]}:0.0);` : ";";
    const width = records[0].sequence.length;
    if (!width || records.some((record) => record.sequence.length !== width)) throw new Error("FastTree input must be a rectangular alignment.");
    const numeric = records.map((record, index) => `>${index}\n${record.sequence}`).join("\n") + "\n";
    const runtime = await modulePromise, path = `/webporpid-fasttree-${invocation++}.fa`;
    output = []; errors = []; runtime.FS.writeFile(path, numeric);
    try { runtime.callMain(["-nosupport", "-nt", "-gtr", path]); }
    finally { try { runtime.FS.unlink(path); } catch { /* best effort */ } }
    let tree = completeNewick(output.join("\n"));
    records.forEach((record, index) => {
      tree = tree.replace(new RegExp(`([,(])${index}(?=[:),])`, "g"), `$1${safeNames[index]}`);
    });
    return tree;
  };
}

class FastTreeWorkerClient {
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

export function createFastTreeRunner(javascriptPath, wasmPath, size = 1,
  workerPath = new URL("../porpid-fasttree-worker.mjs", import.meta.url)) {
  const count = Math.max(1, Math.floor(size));
  if (count === 1) {
    const direct = createDirectFastTreeRunner(javascriptPath, wasmPath); direct.close = async () => {}; return direct;
  }
  const clients = Array.from({ length: count }, () => new FastTreeWorkerClient(new NodeWorker(workerPath)));
  let cursor = 0;
  const run = (alignment) => clients[cursor++ % clients.length].call({ javascriptPath, wasmPath, alignment });
  run.close = async () => { await Promise.all(clients.map((client) => client.close())); };
  return run;
}
