import Aioli from "@biowasm/aioli";
import fastTreeWasmUrl from "/biowasm/fasttree/fasttree.wasm?url";
import { parseFasta } from "./config";

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

export async function runFastTree(alignedFasta: string): Promise<string> {
  const records = parseFasta(alignedFasta);
  const safeName = (name: string, index: number) => name.replace(/[^A-Za-z0-9_.|*+\-]/g, "_") || `tip_${index + 1}`;
  if (records.length < 3) return records.length === 2
    ? `(${safeName(records[0].name, 0)}:0.0,${safeName(records[1].name, 1)}:0.0);`
    : records.length === 1 ? `(${safeName(records[0].name, 0)}:0.0);` : ";";
  const width = records[0].sequence.length;
  if (!width || records.some((record) => record.sequence.length !== width)) throw new Error("FastTree input must be a rectangular alignment.");
  const numeric = records.map((record, index) => `>${index}\n${record.sequence}`).join("\n") + "\n";
  const cli = await tools(), path = "/shared/data/webporpid-fasttree.fa";
  await cli.write({ path, buffer: new TextEncoder().encode(numeric) });
  let tree = completeNewick(await cli.exec(`fasttree -nosupport -nt -gtr ${path}`));
  records.forEach((record, index) => {
    const safe = safeName(record.name, index);
    tree = tree.replace(new RegExp(`([,(])${index}(?=[:),])`, "g"), `$1${safe}`);
  });
  return tree;
}
