import { createHash } from "node:crypto";
import { readFile } from "node:fs/promises";
import { resolve } from "node:path";
import { performance } from "node:perf_hooks";
import { pathToFileURL } from "node:url";
import { build } from "rolldown";

const runtimeBuildDirectory = resolve(".build/benchmark-consensus-runtime");
await build({ input: { config: resolve("src/config.ts"), runtime: resolve("src/wasm-runtime.ts") },
  output: { dir: runtimeBuildDirectory, format: "es", entryFileNames: "[name].mjs" } });
const [{ compileConfig }, { decodeConsensusOutput, encodeFamilyModel, makeCutoffs, WebPorpidRuntime }] = await Promise.all([
  import(pathToFileURL(resolve(runtimeBuildDirectory, "config.mjs")).href),
  import(pathToFileURL(resolve(runtimeBuildDirectory, "runtime.mjs")).href),
]);

const familyCount = Math.max(1, Number.parseInt(process.env.WEBPORPID_BENCHMARK_FAMILIES ?? "20000", 10));
const repeats = Math.max(1, Number.parseInt(process.env.WEBPORPID_BENCHMARK_REPEATS ?? "3", 10));
const corePath = process.env.WEBPORPID_CORE_WASM ?? new URL("../public/webporpid.wasm", import.meta.url);
const sample = { name: "benchmark", cdnaPrimer: "AAAaaaNNNNNNNNTTT", secondStrandPrimer: "GGG", panel: "panel.fa",
  panelSequences: [] };
const config = { dataset: "synthetic-consensus-grouping", samples: [sample], contaminationPanel: "contam.fa",
  contaminationPanelSequences: [], parameters: { errorRate: .05, minLength: 20, maxLength: 5000, primerTolerance: 1,
    primerWindow: 200, primerChop: 0, maxReadsPerSample: 0, familySizeThreshold: 1, ldaThreshold: .995,
    contaminationClusterThreshold: .015, contaminationProportionThreshold: .2, contaminationDistanceThreshold: .015,
    contaminationFilter: false, agreementThreshold: .6, artefactFraction: .25, outlierQuantile: .99,
    panelThreshold: 50, functionalMatchThreshold: .7, spoolPartitions: 64, deterministicSeed: 1n } };

function umi(index) {
  const symbols = "ACGT"; let output = "";
  for (let position = 0; position < 8; position++) { output = symbols[index & 3] + output; index >>>= 2; }
  return output;
}

function spoolRecord(tag, name, sequence) {
  const encoder = new TextEncoder(), tagBytes = encoder.encode(tag), nameBytes = encoder.encode(name), sequenceBytes = encoder.encode(sequence);
  const bodyLength = 20 + tagBytes.length + nameBytes.length + sequenceBytes.length * 2;
  const output = new Uint8Array(bodyLength + 4), view = new DataView(output.buffer);
  view.setUint32(0, bodyLength, true); view.setUint16(4, 0, true); view.setUint16(6, tagBytes.length, true);
  view.setUint32(8, nameBytes.length, true); view.setUint32(12, sequenceBytes.length, true); view.setBigUint64(16, 0n, true);
  let offset = 24; output.set(tagBytes, offset); offset += tagBytes.length; output.set(nameBytes, offset); offset += nameBytes.length;
  output.set(sequenceBytes, offset); offset += sequenceBytes.length; output.fill(73, offset); return output;
}

const families = [], records = [];
for (let family = 0; family < familyCount; family++) {
  const tag = umi(family), sequence = `AAA${tag}TTT${"ACGT".repeat(14)}`;
  families.push({ sample: sample.name, sampleIndex: 0, umi: tag, familySize: 3, mostLikelyParent: tag,
    posteriorProbability: 1, logOffspringProbability: Number.NEGATIVE_INFINITY, disposition: "likely_real" });
  for (let replicate = 0; replicate < 3; replicate++) records.push(spoolRecord(tag, `read-${family}-${replicate}`, sequence));
}
const partitionLength = records.reduce((sum, record) => sum + record.byteLength, 0), partition = new Uint8Array(partitionLength);
let partitionOffset = 0; for (const record of records) { partition.set(record, partitionOffset); partitionOffset += record.byteLength; }

const module = await WebAssembly.compile(await readFile(corePath));
const runtime = await WebPorpidRuntime.create(module, compileConfig(config));
runtime.initFamilyModel(encodeFamilyModel(families));
const cutoffs = makeCutoffs([BigInt(records.length)], 0);
runtime.consensus(partition, cutoffs); // Compile/hydrate hot functions before timing.
const seconds = []; let output;
for (let repeat = 0; repeat < repeats; repeat++) {
  const started = performance.now(); output = runtime.consensus(partition, cutoffs); seconds.push((performance.now() - started) / 1000);
}
const decoded = decodeConsensusOutput(output, config).consensuses.sort((left, right) => left.id.localeCompare(right.id));
const canonical = decoded.map((row) => [row.id, row.sequence, row.minimumAgreement, row.familySize]);
console.log(JSON.stringify({ runtime: process.version, core: String(corePath), families: familyCount, reads: records.length,
  repeats, seconds, bestSeconds: Math.min(...seconds), outputItems: decoded.length,
  outputSha256: createHash("sha256").update(JSON.stringify(canonical)).digest("hex") }, null, 2));
