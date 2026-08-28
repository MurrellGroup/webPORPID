import Aioli from "@biowasm/aioli";
import fastTreeWasmUrl from "/biowasm/fasttree/fasttree.wasm?url";
import { parseFasta } from "./config";
import { treeTipNames } from "./tree-names";

interface AioliRuntime {
  write(options: { path: string; buffer: Uint8Array }): Promise<void>;
  exec(command: string): Promise<string>;
}
let runtime: Promise<AioliRuntime> | undefined;
const tools = () => {
  const urlPrefix = new URL(".", new URL(fastTreeWasmUrl, globalThis.location.href)).href.replace(/\/$/, "");
  return (runtime ??= Promise.resolve(new Aioli([{ tool: "fasttree", program: "fasttree", version: "2.1.11", urlPrefix }]) as unknown as Promise<AioliRuntime>));
};

function completeNewick(output: string) {
  const candidates = output.split(/\r?\n/).filter((line) => line.includes("(") && line.includes(";")).reverse();
  for (const candidate of candidates) {
    const start = candidate.indexOf("("), end = candidate.lastIndexOf(";");
    if (start >= 0 && end > start) return candidate.slice(start, end + 1);
  }
  throw new Error("FastTree did not return a complete Newick tree.");
}

export async function runFastTree(alignedFasta: string, alphabet: "nt" | "aa" = "nt"): Promise<string> {
  const records = parseFasta(alignedFasta);
  const safeNames = treeTipNames(records.map((record) => record.name));
  if (records.length < 3) return records.length === 2
    ? `(${safeNames[0]}:0.0,${safeNames[1]}:0.0);`
    : records.length === 1 ? `(${safeNames[0]}:0.0);` : ";";
  const width = records[0].sequence.length;
  if (!width || records.some((record) => record.sequence.length !== width)) throw new Error("FastTree input must be a rectangular alignment.");
  const numeric = records.map((record, index) => `>${index}\n${record.sequence}`).join("\n") + "\n";
  const cli = await tools(), path = "/shared/data/webporpid-fasttree.fa";
  await cli.write({ path, buffer: new TextEncoder().encode(numeric) });
  let tree = completeNewick(await cli.exec(alphabet === "nt" ? `fasttree -nosupport -nt -gtr ${path}` : `fasttree -nosupport ${path}`));
  records.forEach((record, index) => {
    tree = tree.replace(new RegExp(`([,(])${index}(?=[:),])`, "g"), `$1${safeNames[index]}`);
  });
  return tree;
}

/** Isolate FastTree so independent samples can use separate CPU cores safely. */
export function runFastTreeIsolated(alignedFasta: string): Promise<string> {
  const worker = new Worker(new URL("./fasttree-worker.ts", import.meta.url), { type: "module" });
  return new Promise((resolve, reject) => {
    const finish = () => worker.terminate();
    worker.onmessage = (event: MessageEvent<{ result?: string; error?: string }>) => {
      finish(); event.data.error ? reject(new Error(event.data.error)) : resolve(event.data.result ?? "");
    };
    worker.onerror = (event) => { finish(); reject(new Error(event.message || "FastTree worker failed.")); };
    worker.postMessage({ fasta: alignedFasta });
  });
}
