#!/usr/bin/env bun
import wasmPath from "../public/webporpid.wasm" with { type: "file" };
import msaPath from "../public/alivibe-msa.wasm" with { type: "file" };
import fastTreeJavascriptPath from "../public/biowasm/fasttree/fasttree.cjs" with { type: "file" };
import fastTreeWasmPath from "../public/biowasm/fasttree/fasttree.wasm" with { type: "file" };
import { runCli } from "./porpid-cli.mjs";
const msaWorkerPath = new URL("./porpid-msa-worker.mjs", import.meta.url);
const fastTreeWorkerPath = new URL("./porpid-fasttree-worker.mjs", import.meta.url);
runCli({ wasmPath, msaPath, fastTreeJavascriptPath, fastTreeWasmPath, msaWorkerPath, fastTreeWorkerPath })
  .catch((cause) => { process.stderr.write(`porpid-cli: ${cause instanceof Error ? cause.message : String(cause)}\n`); process.exitCode = 1; });
