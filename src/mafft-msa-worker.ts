/// <reference lib="webworker" />

import mafftJavascriptUrl from "/biowasm/mafft/disttbfast.mjs?url";
import mafftWasmUrl from "/biowasm/mafft/disttbfast.wasm?url";

interface MafftRuntime {
  FS: {
    writeFile(path: string, contents: Uint8Array): void;
    unlink(path: string): void;
  };
  callMain(arguments_: string[]): number;
}

type MafftFactory = (options: {
  wasmBinary: ArrayBuffer;
  noInitialRun: boolean;
  print(line: string): void;
  printErr(line: string): void;
}) => Promise<MafftRuntime>;

interface MafftRequest { sequences: string[] }

const encoder = new TextEncoder();

function numericFasta(sequences: readonly string[]) {
  return sequences.map((sequence, index) => `>${index}\n${sequence}\n`).join("");
}

function parseNumericAlignment(source: string, expected: readonly string[]) {
  const rows = new Map<number, string>();
  let index = -1, sequence = "";
  const finish = () => {
    if (index < 0) return;
    if (rows.has(index)) throw new Error("MAFFT returned a duplicate sequence identifier.");
    rows.set(index, sequence.toUpperCase());
  };
  for (const raw of source.split(/\r?\n/)) {
    const line = raw.trim();
    if (!line) continue;
    if (line.startsWith(">")) {
      finish(); index = Number(line.slice(1).trim()); sequence = "";
      if (!Number.isSafeInteger(index) || index < 0 || index >= expected.length)
        throw new Error("MAFFT returned an unknown sequence identifier.");
    } else {
      if (index < 0) throw new Error("MAFFT returned sequence data before its first header.");
      sequence += line.replace(/\s/g, "");
    }
  }
  finish();
  if (rows.size !== expected.length) throw new Error("MAFFT returned the wrong number of sequences.");
  const aligned = expected.map((input, row) => {
    const output = rows.get(row)!;
    if (output.replaceAll("-", "") !== input.toUpperCase()) throw new Error(`MAFFT changed sequence ${row + 1}.`);
    return output;
  });
  const width = aligned[0]?.length ?? 0;
  if (!width || aligned.some((row) => row.length !== width)) throw new Error("MAFFT returned a non-rectangular alignment.");
  return aligned;
}

const FFT_NS_2_ARGUMENTS = [
  "-q", "0", "-E", "2", "-V", "-1.53", "-s", "0.0", "-W", "6", "-O",
  "-C", "0-0", "-D", "-b", "62", "-g", "0", "-f", "-1.53", "-Q", "100.0",
  "-h", "0", "-F", "-X", "0.1", "-x", "-1", "-i", "/input.fa",
] as const;

self.addEventListener("message", (event: MessageEvent<MafftRequest>) => {
  void (async () => {
    const sequences = event.data.sequences.map((sequence) => sequence.toUpperCase());
    if (sequences.length < 2) { self.postMessage({ type: "result", result: sequences }); return; }
    if (sequences.some((sequence) => !sequence.length || /[^A-Z?*.-]/.test(sequence)))
      throw new Error("MAFFT input contains an unsupported symbol.");
    const [module, wasmResponse] = await Promise.all([
      import(/* @vite-ignore */ new URL(mafftJavascriptUrl, self.location.href).href),
      fetch(new URL(mafftWasmUrl, self.location.href)),
    ]);
    if (!wasmResponse.ok) throw new Error("The bundled MAFFT WebAssembly module could not be loaded.");
    const factory = module.default as MafftFactory;
    if (typeof factory !== "function") throw new Error("The bundled MAFFT JavaScript runtime is invalid.");
    const stdout: string[] = [];
    let lastProgress = 0;
    const runtime = await factory({
      wasmBinary: await wasmResponse.arrayBuffer(), noInitialRun: true,
      print: (line) => stdout.push(String(line)),
      printErr: (line) => {
        const detail = String(line).replace(/[\b\r]+/g, " ").replace(/\s+/g, " ").trim();
        const now = performance.now();
        if (detail && (now - lastProgress >= 200 || /done|Progressive alignment|distance matrix/i.test(detail))) {
          lastProgress = now; self.postMessage({ type: "progress", detail });
        }
      },
    });
    runtime.FS.writeFile("/input.fa", encoder.encode(numericFasta(sequences)));
    try {
      const status = runtime.callMain([...FFT_NS_2_ARGUMENTS]);
      if (status) throw new Error(`MAFFT exited with status ${status}.`);
    } finally { try { runtime.FS.unlink("/input.fa"); } catch { /* best effort */ } }
    self.postMessage({ type: "result", result: parseNumericAlignment(stdout.join("\n"), sequences) });
  })().catch((cause) => self.postMessage({ type: "error", message: cause instanceof Error ? cause.message : String(cause) }));
});
