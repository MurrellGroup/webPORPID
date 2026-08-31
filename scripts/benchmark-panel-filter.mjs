import { createHash } from "node:crypto";
import { performance } from "node:perf_hooks";
import { filterQueriesAgainstPanel } from "../src/independent-panel-filter.ts";

const queryCount = Math.max(1, Number(process.env.WEBPORPID_PANEL_QUERIES ?? 500));
const sequenceLength = Math.max(100, Number(process.env.WEBPORPID_PANEL_LENGTH ?? 2_400));
const panelSize = Math.max(1, Number(process.env.WEBPORPID_PANEL_ROWS ?? 12));
let state = 0x7f4a7c15;
const random = () => { state ^= state << 13; state ^= state >>> 17; state ^= state << 5; return (state >>> 0) / 0x1_0000_0000; };
const bases = "ACGT";
let reference = ""; for (let index = 0; index < sequenceLength; index++) reference += bases[Math.floor(random() * 4)];
const mutate = (source, rate) => [...source].map((base) => random() < rate ? bases[Math.floor(random() * 4)] : base).join("");
const panel = Array.from({ length: panelSize }, () => mutate(reference, 0.01));
const queries = Array.from({ length: queryCount }, (_, index) => {
  let sequence = mutate(reference, 0.015);
  if (index % 11 === 0) sequence = `${sequence.slice(0, Math.floor(sequence.length / 2))}${"A".repeat(40)}${sequence.slice(Math.floor(sequence.length / 2))}`;
  if (index % 17 === 0) sequence = `${sequence.slice(0, 800)}${sequence.slice(835)}`;
  return sequence;
});

filterQueriesAgainstPanel(queries.slice(0, Math.min(5, queries.length)), panel);
const started = performance.now(), result = filterQueriesAgainstPanel(queries, panel), seconds = (performance.now() - started) / 1000;
console.log(JSON.stringify({ runtime: process.version, queries: queryCount, sequenceLength, panelRows: panelSize, seconds,
  queriesPerSecond: queryCount / seconds, basesPerSecond: queries.reduce((sum, row) => sum + row.length, 0) / seconds,
  resultSha256: createHash("sha256").update(JSON.stringify(result)).digest("hex") }, null, 2));
