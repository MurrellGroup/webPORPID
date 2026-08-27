import YAML from "yaml";
import { BinaryWriter } from "./binary";
import type { NamedSequence, PipelineConfig, PipelineParameters, ResultConfig, SampleConfig } from "./types";

export const DEFAULT_PARAMETERS: PipelineParameters = {
  errorRate: 0.05,
  minLength: 2100,
  maxLength: 4300,
  primerTolerance: 1,
  primerWindow: 200,
  primerChop: 0,
  maxReadsPerSample: 100_000,
  familySizeThreshold: 1,
  ldaThreshold: 0.995,
  contaminationClusterThreshold: 0.015,
  contaminationProportionThreshold: 0.2,
  contaminationDistanceThreshold: 0.015,
  contaminationFilter: true,
  agreementThreshold: 0.6,
  artefactFraction: 0.25,
  outlierQuantile: 0.99,
  panelThreshold: 50,
  functionalMatchThreshold: 0.7,
  spoolPartitions: 64,
  deterministicSeed: 0x504f_5250_4944n,
};

type Mapping = Record<string, unknown>;
const mapping = (value: unknown, context: string): Mapping => {
  if (!value || typeof value !== "object" || Array.isArray(value)) throw new Error(`${context} must be a YAML mapping.`);
  return value as Mapping;
};
const text = (value: unknown, context: string, optional = false): string => {
  if (optional && value == null) return "";
  if (typeof value !== "string" || (!optional && !value.trim())) throw new Error(`${context} must be text.`);
  return value;
};
const optionalNumber = (value: unknown, context: string): number | undefined => {
  if (value == null) return undefined;
  const parsed = Number(value); if (!Number.isFinite(parsed)) throw new Error(`${context} must be numeric.`); return parsed;
};
const pick = (object: Mapping, ...keys: string[]) => keys.find((key) => object[key] !== undefined) ? object[keys.find((key) => object[key] !== undefined)!] : undefined;

function parseParameters(source: Mapping | undefined): PipelineParameters {
  const result = { ...DEFAULT_PARAMETERS };
  if (!source) return result;
  const assignNumber = (property: keyof PipelineParameters, ...keys: string[]) => {
    const value = pick(source, ...keys); if (value !== undefined) (result[property] as number) = Number(value);
  };
  assignNumber("errorRate", "errorRate", "error_rate");
  assignNumber("minLength", "minLength", "min_length");
  assignNumber("maxLength", "maxLength", "max_length");
  assignNumber("primerTolerance", "primerTolerance", "primer_tolerance", "primer_tol");
  assignNumber("primerWindow", "primerWindow", "primer_window");
  assignNumber("primerChop", "primerChop", "primer_chop");
  assignNumber("maxReadsPerSample", "maxReadsPerSample", "max_reads_per_sample", "max_reads");
  assignNumber("familySizeThreshold", "familySizeThreshold", "family_size_threshold", "fs_thresh");
  assignNumber("ldaThreshold", "ldaThreshold", "lda_threshold", "lda_thresh");
  assignNumber("contaminationClusterThreshold", "contaminationClusterThreshold", "contamination_cluster_threshold", "cluster_thresh");
  assignNumber("contaminationProportionThreshold", "contaminationProportionThreshold", "contamination_proportion_threshold", "proportion_thresh");
  assignNumber("contaminationDistanceThreshold", "contaminationDistanceThreshold", "contamination_distance_threshold", "dist_thresh");
  assignNumber("agreementThreshold", "agreementThreshold", "agreement_threshold", "agreement_thresh");
  assignNumber("artefactFraction", "artefactFraction", "artefact_fraction", "af_thresh");
  assignNumber("outlierQuantile", "outlierQuantile", "outlier_quantile", "q_thresh");
  assignNumber("panelThreshold", "panelThreshold", "panel_threshold", "panel_thresh");
  assignNumber("functionalMatchThreshold", "functionalMatchThreshold", "functional_match_threshold", "ff_match");
  assignNumber("spoolPartitions", "spoolPartitions", "spool_partitions");
  const contamination = pick(source, "contaminationFilter", "contamination_filter", "contam_toggle");
  if (contamination !== undefined) result.contaminationFilter = contamination === true || contamination === "on" || contamination === "true";
  const seed = pick(source, "deterministicSeed", "deterministic_seed");
  if (seed !== undefined) result.deterministicSeed = BigInt(String(seed));
  for (const [key, value] of Object.entries(result)) {
    if (key === "deterministicSeed" || key === "contaminationFilter") continue;
    if (!Number.isFinite(value)) throw new Error(`Parameter ${key} must be numeric.`);
  }
  if (result.deterministicSeed < 0n || result.deterministicSeed > 0xffff_ffff_ffff_ffffn)
    throw new Error("deterministicSeed must be an unsigned 64-bit integer.");
  const probabilities: Array<keyof PipelineParameters> = ["ldaThreshold", "contaminationProportionThreshold", "agreementThreshold",
    "artefactFraction", "outlierQuantile", "functionalMatchThreshold"];
  for (const key of probabilities) if (!(Number(result[key]) >= 0 && Number(result[key]) <= 1)) throw new Error(`${key} must be between 0 and 1.`);
  const integers: Array<keyof PipelineParameters> = ["minLength", "maxLength", "primerTolerance", "primerWindow", "primerChop",
    "maxReadsPerSample", "familySizeThreshold", "spoolPartitions"];
  for (const key of integers) if (!Number.isSafeInteger(Number(result[key])) || Number(result[key]) < 0) throw new Error(`${key} must be a non-negative integer.`);
  if (!(result.errorRate >= 0 && result.errorRate < 1) || result.minLength >= result.maxLength) throw new Error("Invalid error-rate or read-length bounds.");
  if (result.familySizeThreshold < 1) throw new Error("familySizeThreshold must be at least 1.");
  if (result.contaminationClusterThreshold < 0 || result.contaminationDistanceThreshold < 0 || result.panelThreshold < 0)
    throw new Error("Distance and panel thresholds cannot be negative.");
  if (result.primerWindow < 16) throw new Error("primerWindow must be at least 16.");
  if (result.spoolPartitions < 1 || result.spoolPartitions > 256 || (result.spoolPartitions & (result.spoolPartitions - 1)) !== 0)
    throw new Error("spoolPartitions must be a power of two from 1 to 256.");
  return result;
}

function parseSamples(source: Mapping): SampleConfig[] {
  return Object.entries(source).map(([name, raw]) => {
    const value = mapping(raw, `Sample ${name}`);
    return {
      name,
      cdnaPrimer: text(pick(value, "cDNA_primer", "cdna_primer", "cdnaPrimer"), `${name}.cDNA_primer`),
      secondStrandPrimer: text(pick(value, "sec_str_primer", "secondStrandPrimer"), `${name}.sec_str_primer`),
      panel: text(value.panel, `${name}.panel`),
      functionalReference: text(pick(value, "ff_ref", "functionalReference"), `${name}.ff_ref`, true) || undefined,
      panelSequences: [],
      familySizeOverride: optionalNumber(pick(value, "fs_override", "familySizeOverride"), `${name}.fs_override`),
      artefactFractionOverride: optionalNumber(pick(value, "af_override", "artefactFractionOverride"), `${name}.af_override`),
      outlierQuantileOverride: optionalNumber(pick(value, "q_override", "outlierQuantileOverride"), `${name}.q_override`),
      agreementOverride: optionalNumber(pick(value, "ma_override", "agreementOverride"), `${name}.ma_override`),
      functionalMatchOverride: optionalNumber(pick(value, "ff_match_override", "functionalMatchOverride"), `${name}.ff_match_override`),
    };
  });
}

export function parseConfigYaml(source: string): PipelineConfig {
  const parsed = mapping(YAML.parse(source), "Configuration");
  let dataset: string, samples: Mapping, parameters: Mapping | undefined;
  if (typeof parsed.dataset === "string" && parsed.samples) {
    dataset = parsed.dataset; samples = mapping(parsed.samples, "samples");
    parameters = parsed.parameters ? mapping(parsed.parameters, "parameters") : undefined;
  } else {
    const entries = Object.entries(parsed);
    if (entries.length !== 1) throw new Error("Original PORPID YAML must contain exactly one dataset.");
    dataset = entries[0][0]; samples = mapping(entries[0][1], dataset);
  }
  const contaminationPanel = typeof parsed.contaminationPanel === "string" ? parsed.contaminationPanel
    : typeof parsed.contamination_panel === "string" ? parsed.contamination_panel : "panels/contam_panel.fasta";
  const config: PipelineConfig = { dataset, samples: parseSamples(samples), contaminationPanel,
    contaminationPanelSequences: [], parameters: parseParameters(parameters) };
  if (!config.dataset.trim() || !config.samples.length) throw new Error("The configuration needs a dataset name and at least one sample.");
  const assays = new Set<string>();
  for (const sample of config.samples) {
    if (!sample.name.trim()) throw new Error("Every sample needs a non-empty name.");
    const id = sample.cdnaPrimer.match(/[a-z]+/)?.[0].toUpperCase();
    if (!id) throw new Error(`${sample.name} cDNA_primer has no lower-case sample ID.`);
    const assay = `${sample.secondStrandPrimer.toUpperCase()}\0${sample.cdnaPrimer.toUpperCase()}`;
    if (assays.has(assay)) throw new Error(`${sample.name} duplicates another sample's complete primer pair.`); assays.add(assay);
    if (!/[Nn]+/.test(sample.cdnaPrimer)) throw new Error(`${sample.name} cDNA_primer has no N-marked UMI.`);
    if (sample.familySizeOverride !== undefined && (!Number.isSafeInteger(sample.familySizeOverride) || sample.familySizeOverride < 1))
      throw new Error(`${sample.name}.fs_override must be an integer of at least 1.`);
    for (const [label, value] of [["af_override", sample.artefactFractionOverride], ["q_override", sample.outlierQuantileOverride],
      ["ma_override", sample.agreementOverride], ["ff_match_override", sample.functionalMatchOverride]] as const)
      if (value !== undefined && !(value >= 0 && value <= 1)) throw new Error(`${sample.name}.${label} must be between 0 and 1.`);
  }
  return config;
}

export function parseFasta(source: string): NamedSequence[] {
  const records: NamedSequence[] = []; let name = "", sequence = "";
  for (const rawLine of source.split(/\r?\n/)) {
    const line = rawLine.trim(); if (!line || line.startsWith(";")) continue;
    if (line.startsWith(">")) {
      if (name) records.push({ name, sequence: sequence.toUpperCase() });
      name = line.slice(1).trim() || `sequence_${records.length + 1}`; sequence = "";
    } else {
      if (!name) throw new Error("FASTA sequence data appeared before its first header.");
      sequence += line.replace(/\s/g, "");
    }
  }
  if (name) records.push({ name, sequence: sequence.toUpperCase() });
  if (!records.length) throw new Error("The FASTA file contains no records.");
  return records;
}

const normalizePath = (value: string) => value.replaceAll("\\", "/").replace(/^\.\//, "");
const basename = (value: string) => normalizePath(value).split("/").at(-1)!;

export async function resolveReferenceFiles(config: PipelineConfig, files: Map<string, () => Promise<string>>): Promise<PipelineConfig> {
  const cache = new Map<string, NamedSequence[]>();
  const find = (requested: string) => {
    const normalized = normalizePath(requested);
    const candidates = [...files.entries()].filter(([name]) => normalizePath(name) === normalized || basename(name) === basename(normalized));
    if (candidates.length !== 1) throw new Error(candidates.length ? `Reference ${requested} is ambiguous.` : `Reference ${requested} was not supplied.`);
    return candidates[0];
  };
  const load = async (requested: string) => {
    if (cache.has(requested)) return cache.get(requested)!;
    const [, content] = find(requested); const parsed = parseFasta(await content()); cache.set(requested, parsed); return parsed;
  };
  for (const sample of config.samples) {
    sample.panelSequences = await load(sample.panel);
    if (sample.functionalReference) sample.functionalReferenceSequence = (await load(sample.functionalReference))[0];
  }
  config.contaminationPanelSequences = config.parameters.contaminationFilter ? await load(config.contaminationPanel) : [];
  return config;
}

const nan = (value: number | undefined) => value ?? Number.NaN;

export function compileConfig(config: PipelineConfig): Uint8Array {
  const writer = new BinaryWriter(); const p = config.parameters;
  writer.magic("WPC1"); writer.u32(1); writer.string(config.dataset);
  writer.f64(p.errorRate); writer.u32(p.minLength); writer.u32(p.maxLength); writer.u32(p.primerTolerance);
  writer.u32(p.primerWindow); writer.u32(p.primerChop); writer.u32(p.maxReadsPerSample); writer.u32(p.familySizeThreshold);
  writer.f64(p.ldaThreshold); writer.f64(p.contaminationClusterThreshold); writer.f64(p.contaminationProportionThreshold);
  writer.f64(p.contaminationDistanceThreshold); writer.u8(p.contaminationFilter ? 1 : 0); writer.f64(p.agreementThreshold);
  writer.f64(p.artefactFraction); writer.f64(p.outlierQuantile); writer.f64(p.panelThreshold); writer.f64(p.functionalMatchThreshold);
  writer.u32(p.spoolPartitions); writer.u64(p.deterministicSeed); writer.u32(config.samples.length);
  for (const sample of config.samples) {
    writer.string(sample.name); writer.string(sample.cdnaPrimer); writer.string(sample.secondStrandPrimer);
    writer.string(sample.panel); writer.string(sample.functionalReference ?? ""); writer.i32(sample.familySizeOverride ?? -1);
    writer.f64(nan(sample.artefactFractionOverride)); writer.f64(nan(sample.outlierQuantileOverride));
    writer.f64(nan(sample.agreementOverride)); writer.f64(nan(sample.functionalMatchOverride));
    writer.u32(sample.panelSequences.length);
    for (const sequence of sample.panelSequences) { writer.string(sequence.name); writer.string(sequence.sequence); }
    writer.string(sample.functionalReferenceSequence?.sequence ?? "");
  }
  return writer.finish();
}

export function resultConfig(config: PipelineConfig): ResultConfig {
  return {
    dataset: config.dataset,
    samples: config.samples.map(({ panelSequences: _panel, functionalReferenceSequence: _reference, ...sample }) => sample),
    contaminationPanel: config.contaminationPanel,
    parameters: { ...config.parameters, deterministicSeed: config.parameters.deterministicSeed.toString() },
  };
}

export function serializeConfigYaml(config: PipelineConfig): string {
  const samples: Record<string, Record<string, unknown>> = {};
  for (const sample of config.samples) {
    if (!sample.name.trim()) throw new Error("Every sample needs a name before YAML can be generated.");
    if (Object.hasOwn(samples, sample.name)) throw new Error(`Sample name ${sample.name} is duplicated.`);
    const row: Record<string, unknown> = {
      cDNA_primer: sample.cdnaPrimer,
      sec_str_primer: sample.secondStrandPrimer,
      panel: sample.panel,
    };
    if (sample.functionalReference) row.ff_ref = sample.functionalReference;
    if (sample.familySizeOverride !== undefined) row.fs_override = sample.familySizeOverride;
    if (sample.artefactFractionOverride !== undefined) row.af_override = sample.artefactFractionOverride;
    if (sample.outlierQuantileOverride !== undefined) row.q_override = sample.outlierQuantileOverride;
    if (sample.agreementOverride !== undefined) row.ma_override = sample.agreementOverride;
    if (sample.functionalMatchOverride !== undefined) row.ff_match_override = sample.functionalMatchOverride;
    samples[sample.name] = row;
  }
  return YAML.stringify({
    dataset: config.dataset,
    samples,
    contaminationPanel: config.contaminationPanel,
    parameters: { ...config.parameters, deterministicSeed: config.parameters.deterministicSeed.toString() },
  }, { lineWidth: 0 });
}

export function blankConfig(): PipelineConfig {
  return {
    dataset: "analysis",
    samples: [{ name: "sample_1", cdnaPrimer: "", secondStrandPrimer: "", panel: "panels/panel.fasta", panelSequences: [] }],
    contaminationPanel: "panels/contam_panel.fasta", contaminationPanelSequences: [], parameters: { ...DEFAULT_PARAMETERS },
  };
}
