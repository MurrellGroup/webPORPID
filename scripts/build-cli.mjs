import { chmod, copyFile, mkdir, rm } from "node:fs/promises";
import { build } from "rolldown";

const root = new URL("../", import.meta.url), output = new URL("../cli/", import.meta.url);
await rm(output, { recursive: true, force: true });
await build({
  input: {
    "porpid-cli": new URL("../cli-src/porpid-cli-node.mjs", import.meta.url).pathname,
    "porpid-worker": new URL("../cli-src/porpid-worker.mjs", import.meta.url).pathname,
    "porpid-msa-worker": new URL("../cli-src/porpid-msa-worker.mjs", import.meta.url).pathname,
    "porpid-fasttree-worker": new URL("../cli-src/porpid-fasttree-worker.mjs", import.meta.url).pathname,
  },
  external: /^node:/,
  output: { dir: output.pathname, format: "es", entryFileNames: "[name].mjs", chunkFileNames: "chunks/[name]-[hash].mjs" },
});
const assets = new URL("assets/", output); await mkdir(assets, { recursive: true });
for (const [source, target] of [
  ["public/webporpid.wasm", "webporpid.wasm"], ["public/alivibe-msa.wasm", "alivibe-msa.wasm"],
  ["public/biowasm/fasttree/fasttree.cjs", "fasttree.cjs"], ["public/biowasm/fasttree/fasttree.wasm", "fasttree.wasm"],
]) await copyFile(new URL(source, root), new URL(target, assets));
await chmod(new URL("porpid-cli.mjs", output), 0o755);
console.log("Built cli/porpid-cli.mjs with local WASM assets");
