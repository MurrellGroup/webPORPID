import { readFile } from "node:fs/promises";
import { parseFasta } from "../src/config.ts";

function completeNewick(output) {
  const candidates = output.split(/\r?\n/).filter((line) => line.includes("(") && line.includes(";")).reverse();
  for (const candidate of candidates) {
    const start = candidate.indexOf("("), end = candidate.lastIndexOf(";");
    if (start >= 0 && end > start) return candidate.slice(start, end + 1);
  }
  throw new Error("FastTree did not return a complete Newick tree.");
}

export function createFastTreeRunner(javascriptPath, wasmPath) {
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
    const safeName = (name, index) => name.replace(/[^A-Za-z0-9_.|*+\-]/g, "_") || `tip_${index + 1}`;
    if (records.length < 3) return records.length === 2
      ? `(${safeName(records[0].name, 0)}:0.0,${safeName(records[1].name, 1)}:0.0);`
      : records.length === 1 ? `(${safeName(records[0].name, 0)}:0.0);` : ";";
    const width = records[0].sequence.length;
    if (!width || records.some((record) => record.sequence.length !== width)) throw new Error("FastTree input must be a rectangular alignment.");
    const numeric = records.map((record, index) => `>${index}\n${record.sequence}`).join("\n") + "\n";
    const runtime = await modulePromise, path = `/webporpid-fasttree-${invocation++}.fa`;
    output = []; errors = []; runtime.FS.writeFile(path, numeric);
    try { runtime.callMain(["-nosupport", "-nt", "-gtr", path]); }
    finally { try { runtime.FS.unlink(path); } catch { /* best effort */ } }
    let tree = completeNewick(output.join("\n"));
    records.forEach((record, index) => {
      const safe = safeName(record.name, index);
      tree = tree.replace(new RegExp(`([,(])${index}(?=[:),])`, "g"), `$1${safe}`);
    });
    return tree;
  };
}
