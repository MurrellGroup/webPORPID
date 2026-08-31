/// <reference lib="webworker" />

import { filterQueriesAgainstPanel } from "./independent-panel-filter.ts";

interface Request { sequences: string[]; panelRows: string[] }

self.addEventListener("message", (event: MessageEvent<Request>) => {
  try {
    const result = filterQueriesAgainstPanel(event.data.sequences, event.data.panelRows,
      (progress) => self.postMessage({ type: "progress", ...progress }));
    self.postMessage({ type: "result", result });
  } catch (cause) {
    self.postMessage({ type: "error", message: cause instanceof Error ? cause.message : String(cause) });
  }
});
