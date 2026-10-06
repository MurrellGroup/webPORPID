/// <reference lib="webworker" />

import { mafftArguments, numericMafftFasta, parseMafftAlignment } from "./mafft-msa-codec.ts";
import type { MafftRequest } from "./mafft-msa-codec.ts";

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

const encoder = new TextEncoder();

self.addEventListener("message", (event: MessageEvent<MafftRequest>) => {
  void (async () => {
    const sequences = event.data.sequences.map((sequence) => sequence.toUpperCase());
    const { scoringMode, profileRows } = event.data;
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
    runtime.FS.writeFile("/input.fa", encoder.encode(numericMafftFasta(sequences, scoringMode)));
    try {
      const status = runtime.callMain(mafftArguments(scoringMode, profileRows));
      if (status) throw new Error(`MAFFT exited with status ${status}.`);
    } finally { try { runtime.FS.unlink("/input.fa"); } catch { /* best effort */ } }
    self.postMessage({ type: "result", result: parseMafftAlignment(stdout.join("\n"), sequences, scoringMode) });
  })().catch((cause) => self.postMessage({ type: "error", message: cause instanceof Error ? cause.message : String(cause) }));
});
