import { performance } from "node:perf_hooks";
import { createHash } from "node:crypto";
const contaminationModule = process.env.WEBPORPID_CONTAMINATION_MODULE ?? new URL("../src/contamination.ts", import.meta.url).href;
const postprocessModule = process.env.WEBPORPID_POSTPROCESS_MODULE ?? new URL("../src/postprocess.ts", import.meta.url).href;
const reportStatsModule = process.env.WEBPORPID_REPORT_STATS_MODULE ?? new URL("../src/report-stats.ts", import.meta.url).href;
const { classifyContaminationAsync } = await import(contaminationModule);
const { postprocess } = await import(postprocessModule);
const { sampleOverviewStats } = await import(reportStatsModule);

let randomState = 0x8f31a2c7;
function random() { randomState ^= randomState << 13; randomState ^= randomState >>> 17; randomState ^= randomState << 5; return (randomState >>> 0) / 0x1_0000_0000; }
const bases = "ACGT";
function sequence(length) { let output = ""; for (let index = 0; index < length; index++) output += bases[Math.floor(random() * 4)]; return output; }
function mutate(source, rate) { return [...source].map((base) => random() < rate ? bases[(bases.indexOf(base) + 1 + Math.floor(random() * 3)) % 4] : base).join(""); }
const parameters = { errorRate: .05, minLength: 20, maxLength: 5000, primerTolerance: 1, primerWindow: 200, primerChop: 0,
  maxReadsPerSample: 0, familySizeThreshold: 1, ldaThreshold: .995, contaminationClusterThreshold: .015,
  contaminationProportionThreshold: .2, contaminationDistanceThreshold: .015, contaminationFilter: true,
  agreementThreshold: .6, artefactFraction: .25, outlierQuantile: .99, panelThreshold: 50,
  functionalMatchThreshold: .7, spoolPartitions: 64, deterministicSeed: 1n };
const samples = Array.from({ length: 6 }, (_, index) => ({ name: `sample_${index + 1}`,
  cdnaPrimer: "AAAaaaNNNNNNNNTTT", secondStrandPrimer: "GGG", panel: "panel.fa", panelSequences: [] }));
const founders = Array.from({ length: 12 }, () => sequence(600));
const consensuses = Array.from({ length: 6_000 }, (_, index) => {
  const sampleIndex = index % samples.length, template = founders[(index * 7 + sampleIndex) % founders.length];
  return { id: `c${index}`, sample: samples[sampleIndex].name, sampleIndex, umi: index.toString(16).padStart(8, "0").slice(-8),
    familySize: 1 + index % 20, minimumAgreement: .55 + (index % 46) / 100, sequence: mutate(template, .004 + (index % 4) * .002), lowAgreementSites: [] };
});
const config = { dataset: "synthetic-benchmark", samples, contaminationPanel: "contam.fa",
  contaminationPanelSequences: Array.from({ length: 80 }, (_, index) => ({ name: `panel_${index}`, sequence: mutate(founders[index % founders.length], .02) })), parameters };

async function timed(name, work) {
  const started = performance.now(); const value = await work();
  return { name, seconds: (performance.now() - started) / 1000, outputItems: Array.isArray(value) ? value.length : value,
    ...(Array.isArray(value) ? { outputSha256: createHash("sha256").update(JSON.stringify(value)).digest("hex") } : {}) };
}

const contamination = await timed("contamination-6000", () => classifyContaminationAsync(consensuses, config));
const downstreamConfig = { ...config, parameters: { ...parameters, contaminationFilter: false },
  samples: [samples[0]], contaminationPanelSequences: [] };
const downstreamRows = consensuses.filter((row) => row.sampleIndex === 0).slice(0, 250).map((row) => ({ ...row, sampleIndex: 0, sample: samples[0].name }));
const postprocessing = await timed("postprocess-250", async () => {
  const output = await postprocess(downstreamRows, [], downstreamConfig, undefined, async (rows) => {
    const width = Math.max(...rows.map((row) => row.length)); return rows.map((row) => row.padEnd(width, "-"));
  }, 1, undefined, { collapse: false });
  if (process.env.WEBPORPID_SHOW_APOBEC === "1") console.error(JSON.stringify(output.records.find((row) => row.apobec)?.apobec));
  return output.records.map(({ apobec: _apobec, ...row }) => row);
});

const overviewSamples = Array.from({ length: 100 }, (_, index) => ({ ...samples[0], name: `overview_${index}`, donorId: `d${Math.floor(index / 4)}` }));
const overviewFamilies = Array.from({ length: 100_000 }, (_, index) => ({ sample: overviewSamples[index % overviewSamples.length].name,
  sampleIndex: index % overviewSamples.length, umi: index.toString(16), familySize: 1 + index % 10, mostLikelyParent: "p",
  posteriorProbability: .999, logOffspringProbability: -7, disposition: index % 13 ? "likely_real" : "LDA-rejects" }));
const overviewBundle = { config: { dataset: "overview", samples: overviewSamples, contaminationPanel: "", parameters: { ...parameters, deterministicSeed: "1" } },
  summaries: overviewSamples.map((sample) => ({ sample: sample.name, demultiplexedReads: 5_500, selectedReads: 5_500,
    downsampledReads: 0, observedUmis: 1_000, likelyRealUmis: 923, consensusSequences: 0 })),
  umiFamilies: overviewFamilies, records: [], consensuses: [], contamination: [], collapseGroups: {}, postprocessingContaminationMode: "bypassed" };
const overview = await timed("overview-100000-families", () => sampleOverviewStats(overviewBundle));

console.log(JSON.stringify({ runtime: process.version, platform: `${process.platform}/${process.arch}`, results: [contamination, postprocessing, overview] }, null, 2));
