/// <reference lib="webworker" />
import { runFastTree } from "./biowasm";

self.addEventListener("message", (event: MessageEvent<{ fasta: string }>) => {
  void runFastTree(event.data.fasta).then((result) => self.postMessage({ result }))
    .catch((cause) => self.postMessage({ error: cause instanceof Error ? cause.message : String(cause) }));
});
