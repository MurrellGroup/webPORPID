import { decode, encode } from "@msgpack/msgpack";
import { gzipSync, gunzipSync } from "fflate";
import { inspectAlignment, translateAlignmentFasta, validateCorrectedAlignment } from "./alignment-utils.ts";
import type { ResultBundle } from "./types";

const MAGIC = Uint8Array.of(0x57, 0x50, 0x52, 0x00, 0x01, 0x0d, 0x0a, 0x1a);
const MAX_COMPRESSED_BYTES = 1024 * 1024 * 1024;
const MAX_UNCOMPRESSED_BYTES = 2 * 1024 * 1024 * 1024;
const DISPOSITIONS = new Set(["likely_real", "BPB-rejects", "heteroduplex", "LDA-rejects", "UMI_len != 8", "family-size-reject"]);

type UnknownRecord = Record<string, unknown>;
const object = (value: unknown, label: string): UnknownRecord => {
  if (!value || typeof value !== "object" || Array.isArray(value)) throw new Error(`${label} must be an object.`); return value as UnknownRecord;
};
const array = (value: unknown, label: string): unknown[] => { if (!Array.isArray(value)) throw new Error(`${label} must be an array.`); return value; };
const text = (value: unknown, label: string) => { if (typeof value !== "string") throw new Error(`${label} must be text.`); return value; };
const bool = (value: unknown, label: string) => { if (typeof value !== "boolean") throw new Error(`${label} must be boolean.`); return value; };
const numeric = (value: unknown, label: string, finite = true) => {
  if (typeof value !== "number" || Number.isNaN(value) || (finite && !Number.isFinite(value))) throw new Error(`${label} must be numeric.`); return value;
};
const count = (value: unknown, label: string) => { const result = numeric(value, label); if (!Number.isSafeInteger(result) || result < 0) throw new Error(`${label} must be a non-negative integer.`); return result; };
const optionalText = (value: unknown, label: string) => { if (value != null) text(value, label); };
const optionalNumber = (value: unknown, label: string) => { if (value != null) numeric(value, label); };
const optionalBool = (value: unknown, label: string) => { if (value != null) bool(value, label); };

function validateResult(value: unknown): asserts value is ResultBundle {
  const bundle = object(value, "Results payload");
  if (bundle.schema !== "webporpid-results/1") throw new Error("Unsupported webPORPID result schema.");
  const provenance = object(bundle.provenance, "provenance");
  for (const key of ["webporpidVersion", "createdUtc", "engine", "inputName", "inputSha256", "configSha256", "deterministicSeed", "upstreamBranch", "upstreamCommit"])
    text(provenance[key], `provenance.${key}`);
  if (count(provenance.workers, "provenance.workers") < 1) throw new Error("provenance.workers must be at least one.");

  const config = object(bundle.config, "config"); text(config.dataset, "config.dataset");
  const parameters = object(config.parameters, "config.parameters");
  for (const key of ["errorRate", "minLength", "maxLength", "primerTolerance", "primerWindow", "primerChop", "maxReadsPerSample",
    "familySizeThreshold", "ldaThreshold", "contaminationClusterThreshold", "contaminationProportionThreshold",
    "contaminationDistanceThreshold", "agreementThreshold", "artefactFraction", "outlierQuantile", "panelThreshold",
    "functionalMatchThreshold", "spoolPartitions"])
    numeric(parameters[key], `config.parameters.${key}`);
  bool(parameters.contaminationFilter, "config.parameters.contaminationFilter");
  text(parameters.deterministicSeed, "config.parameters.deterministicSeed");
  const samples = array(config.samples, "config.samples").map((entry, index) => {
    const sample = object(entry, `config.samples[${index}]`); const name = text(sample.name, `config.samples[${index}].name`);
    text(sample.cdnaPrimer, `config.samples[${index}].cdnaPrimer`); text(sample.secondStrandPrimer, `config.samples[${index}].secondStrandPrimer`);
    text(sample.panel, `config.samples[${index}].panel`); optionalText(sample.functionalReference, `config.samples[${index}].functionalReference`); return name;
  });
  if (new Set(samples).size !== samples.length) throw new Error("Result sample names must be unique.");
  const sampleSet = new Set(samples);
  const sampleIndices = new Map(samples.map((name, index) => [name, index]));
  const knownSample = (sample: string, label: string) => {
    if (!sampleSet.has(sample)) throw new Error(`${label} references an unknown sample.`);
  };

  const quality = object(bundle.quality, "quality");
  for (const key of ["totalReads", "qualityReads", "badReads", "shortReads", "longReads", "primerRejects", "idRejects", "demultiplexedReads", "bpbRejects", "malformedRecords", "downsampledReads"])
    count(quality[key], `quality.${key}`);
  const perSample = array(quality.perSample, "quality.perSample"); if (perSample.length !== samples.length) throw new Error("quality.perSample has the wrong sample count.");
  perSample.forEach((entry, index) => count(entry, `quality.perSample[${index}]`));

  const summarySamples = new Set<string>();
  array(bundle.summaries, "summaries").forEach((entry, index) => {
    const row = object(entry, `summaries[${index}]`), sample = text(row.sample, `summaries[${index}].sample`);
    if (!sampleSet.has(sample) || summarySamples.has(sample)) throw new Error("Result summaries contain an unknown or duplicate sample."); summarySamples.add(sample);
    for (const key of ["demultiplexedReads", "observedUmis", "likelyRealUmis", "consensusSequences", "contaminationPassed", "postprocPassed", "artefactCutoff"])
      count(row[key], `summaries[${index}].${key}`);
    if (row.functionalPassed != null) count(row.functionalPassed, `summaries[${index}].functionalPassed`);
  });
  if (summarySamples.size !== samples.length) throw new Error("Result summaries are missing a configured sample.");

  const familyKeys = new Set<string>();
  array(bundle.umiFamilies, "umiFamilies").forEach((entry, index) => {
    const row = object(entry, `umiFamilies[${index}]`), sample = text(row.sample, `umiFamilies[${index}].sample`), sampleIndex = count(row.sampleIndex, `umiFamilies[${index}].sampleIndex`);
    knownSample(sample, `umiFamilies[${index}]`); if (sampleIndices.get(sample) !== sampleIndex) throw new Error("A UMI family has an inconsistent sample index.");
    const umi = text(row.umi, `umiFamilies[${index}].umi`), familyKey = `${sampleIndex}\0${umi}`;
    if (familyKeys.has(familyKey)) throw new Error("UMI family identifiers must be unique within a sample."); familyKeys.add(familyKey);
    count(row.familySize, `umiFamilies[${index}].familySize`); text(row.mostLikelyParent, `umiFamilies[${index}].mostLikelyParent`);
    numeric(row.posteriorProbability, `umiFamilies[${index}].posteriorProbability`); numeric(row.logOffspringProbability, `umiFamilies[${index}].logOffspringProbability`, false);
    if (!DISPOSITIONS.has(text(row.disposition, `umiFamilies[${index}].disposition`))) throw new Error("A UMI family has an unknown disposition.");
    optionalNumber(row.minimumAgreement, `umiFamilies[${index}].minimumAgreement`);
  });

  const consensusIds = new Set<string>();
  const consensusById = new Map<string, UnknownRecord>();
  array(bundle.consensuses, "consensuses").forEach((entry, index) => {
    const row = object(entry, `consensuses[${index}]`), id = text(row.id, `consensuses[${index}].id`);
    if (consensusIds.has(id)) throw new Error("Consensus identifiers must be unique."); consensusIds.add(id); consensusById.set(id, row);
    const sample = text(row.sample, `consensuses[${index}].sample`), sampleIndex = count(row.sampleIndex, `consensuses[${index}].sampleIndex`);
    knownSample(sample, `consensuses[${index}]`); if (sampleIndices.get(sample) !== sampleIndex) throw new Error("A consensus has an inconsistent sample index.");
    text(row.umi, `consensuses[${index}].umi`);
    count(row.familySize, `consensuses[${index}].familySize`); numeric(row.minimumAgreement, `consensuses[${index}].minimumAgreement`); text(row.sequence, `consensuses[${index}].sequence`);
    array(row.lowAgreementSites, `consensuses[${index}].lowAgreementSites`).forEach((site, siteIndex) => {
      const low = object(site, `consensuses[${index}].lowAgreementSites[${siteIndex}]`); count(low.position, "low-agreement position");
      numeric(low.agreement, "low-agreement value"); text(low.modalReadBase, "low-agreement modal base"); count(low.modalRunLength, "low-agreement run length");
    });
  });

  array(bundle.contamination, "contamination").forEach((entry, index) => {
    const row = object(entry, `contamination[${index}]`), sample = text(row.sample, "contamination sample"), sequenceId = text(row.sequenceId, "contamination sequence ID");
    knownSample(sample, `contamination[${index}]`);
    if (consensusById.get(sequenceId)?.sample !== sample) throw new Error("A contamination call references an unknown consensus or sample.");
    text(row.nearestNonselfVariant, "nearest non-self variant"); numeric(row.nearestNonselfDistance, "nearest non-self distance");
    bool(row.flagged, "contamination flagged"); bool(row.discarded, "contamination discarded"); bool(row.suspectOnly, "contamination suspectOnly");
  });

  const recordIds = new Set<string>();
  array(bundle.records, "records").forEach((entry, index) => {
    const row = object(entry, `records[${index}]`), id = text(row.id, `records[${index}].id`);
    if (recordIds.has(id)) throw new Error("Post-processing identifiers must be unique."); recordIds.add(id);
    const sample = text(row.sample, "postproc sample"), source = consensusById.get(id); knownSample(sample, `records[${index}]`);
    if (!source || source.sample !== sample) throw new Error("A post-processing record references an unknown consensus or sample.");
    const umi = text(row.umi, "postproc UMI"), familySize = count(row.familySize, "postproc family size"), minimumAgreement = numeric(row.minimumAgreement, "postproc agreement");
    if (umi !== source.umi || familySize !== source.familySize || minimumAgreement !== source.minimumAgreement)
      throw new Error("A post-processing record has inconsistent consensus metadata.");
    const consensusNt = text(row.consensusNt, "postproc consensus"); if (consensusNt !== source.sequence) throw new Error("A post-processing record has inconsistent consensus sequence data.");
    optionalText(row.alignedNt, "postproc aligned sequence"); optionalText(row.trimmedNt, "postproc trimmed nucleotide"); optionalText(row.trimmedAa, "postproc trimmed protein");
    numeric(row.panelScore, "postproc panel score"); for (const key of ["artefactPass", "agreementPass", "contaminationPass", "panelPass"]) bool(row[key], `postproc ${key}`);
    optionalBool(row.functionalPass, "postproc functionalPass"); array(row.rejectionReasons, "postproc rejectionReasons").forEach((reason) => text(reason, "postproc rejection reason"));
    if (row.apobec != null) { const model = object(row.apobec, "postproc APOBEC"); for (const key of ["posteriorMeanGaMultiplier", "posteriorGaInflated", "posteriorMeanMutationRate", "gaMutations", "totalMutations"]) numeric(model[key], `APOBEC ${key}`); }
  });
  if (recordIds.size !== consensusIds.size || [...consensusIds].some((id) => !recordIds.has(id)))
    throw new Error("Consensus and post-processing records are inconsistent.");

  for (const [label, entries] of [["alignments", object(bundle.alignments, "alignments")], ["trees", object(bundle.trees, "trees")]] as const)
    for (const [name, contents] of Object.entries(entries)) {
      text(name, `${label} name`); text(contents, `${label}.${name}`);
      const sample = name.split("/", 1)[0]; knownSample(sample, `${label}.${name}`);
    }
  if (bundle.alignmentEdits != null) for (const [name, rawEdit] of Object.entries(object(bundle.alignmentEdits, "alignmentEdits"))) {
    const sample = name.split("/", 1)[0]; knownSample(sample, `alignmentEdits.${name}`);
    if (name !== `${sample}/nucleotide`) throw new Error("Edited alignment keys must end in /nucleotide.");
    const edit = object(rawEdit, `alignmentEdits.${name}`), fasta = text(edit.fasta, `alignmentEdits.${name}.fasta`);
    const frameOffset = count(edit.frameOffset, `alignmentEdits.${name}.frameOffset`);
    if (frameOffset > 2) throw new Error("Edited alignment frame offsets must be 0, 1, or 2.");
    const baselineFingerprint = text(edit.baselineFingerprint, `alignmentEdits.${name}.baselineFingerprint`);
    const editedFingerprint = text(edit.editedFingerprint, `alignmentEdits.${name}.editedFingerprint`);
    text(edit.source, `alignmentEdits.${name}.source`); text(edit.savedUtc, `alignmentEdits.${name}.savedUtc`);
    optionalText(edit.treeNewick, `alignmentEdits.${name}.treeNewick`);
    const original = text(object(bundle.alignments, "alignments")[name], `alignments.${name}`);
    if (inspectAlignment(original, 1).fingerprint !== baselineFingerprint) throw new Error(`alignmentEdits.${name} has an inconsistent baseline fingerprint.`);
    if (inspectAlignment(fasta, 1).fingerprint !== editedFingerprint) throw new Error(`alignmentEdits.${name} has an inconsistent fingerprint.`);
    validateCorrectedAlignment(original, fasta);
  }
  if (bundle.timings != null) array(bundle.timings, "timings").forEach((entry, index) => {
    const timing = object(entry, `timings[${index}]`);
    const stage = text(timing.stage, `timings[${index}].stage`);
    if (!["setup", "preprocessing", "umi", "consensus", "contamination", "postprocessing", "tree", "analysis-total"].includes(stage))
      throw new Error(`timings[${index}] has an unknown stage.`);
    const seconds = numeric(timing.seconds, `timings[${index}].seconds`);
    if (seconds < 0) throw new Error(`timings[${index}].seconds must be non-negative.`);
    if (timing.workItems != null) count(timing.workItems, `timings[${index}].workItems`);
  });
  array(bundle.log, "log").forEach((entry, index) => text(entry, `log[${index}]`));
}

export function encodeResult(bundle: ResultBundle): Uint8Array {
  validateResult(bundle); const body = gzipSync(encode(bundle), { level: 9 }), output = new Uint8Array(MAGIC.byteLength + body.byteLength);
  output.set(MAGIC); output.set(body, MAGIC.byteLength); return output;
}

export function decodeResult(bytes: Uint8Array): ResultBundle {
  if (bytes.byteLength < MAGIC.byteLength || MAGIC.some((value, index) => bytes[index] !== value)) throw new Error("This is not a webPORPID results file.");
  const compressed = bytes.subarray(MAGIC.byteLength);
  if (compressed.byteLength > MAX_COMPRESSED_BYTES) throw new Error("The webPORPID results file is too large to load safely.");
  if (compressed.byteLength < 18) throw new Error("The webPORPID results payload is truncated.");
  const footer = new DataView(compressed.buffer, compressed.byteOffset + compressed.byteLength - 4, 4).getUint32(0, true);
  if (footer > MAX_UNCOMPRESSED_BYTES) throw new Error("The uncompressed webPORPID results payload is too large to load safely.");
  let unpacked: Uint8Array;
  try { unpacked = gunzipSync(compressed); } catch { throw new Error("The webPORPID results payload is corrupt or truncated."); }
  if (unpacked.byteLength > MAX_UNCOMPRESSED_BYTES) throw new Error("The uncompressed webPORPID results payload is too large to load safely.");
  let value: unknown;
  try { value = decode(unpacked); } catch { throw new Error("The webPORPID results payload contains invalid MessagePack data."); }
  validateResult(value); return value;
}

const quote = (value: unknown) => {
  const valueText = value == null ? "" : String(value); return /[",\r\n]/.test(valueText) ? `"${valueText.replaceAll('"', '""')}"` : valueText;
};
const csv = (headers: string[], rows: unknown[][]) => [headers, ...rows].map((row) => row.map(quote).join(",")).join("\n") + "\n";
const fasta = (rows: Array<{ id: string; sequence: string }>) => rows.map((row) => `>${row.id}\n${row.sequence.match(/.{1,80}/g)?.join("\n") ?? ""}`).join("\n") + (rows.length ? "\n" : "");
const passed = (record: ResultBundle["records"][number]) => record.artefactPass && record.agreementPass && record.contaminationPass && record.panelPass;

export type ExportKind = "consensus-fasta" | "passed-consensus-fasta" | "rejected-consensus-fasta" | "trimmed-nt-fasta" | "trimmed-aa-fasta"
  | "family-csv" | "low-agreement-csv" | "contamination-csv" | "postproc-csv" | "apobec-csv"
  | "nucleotide-alignment" | "protein-alignment" | "newick" | "log";

function alignmentSample(bundle: ResultBundle, sample?: string) {
  if (sample) return sample;
  if (bundle.config.samples.length === 1) return bundle.config.samples[0].name;
  throw new Error("Choose a sample when exporting a sample-specific alignment or tree.");
}

export function exportComponent(bundle: ResultBundle, kind: ExportKind, sample?: string): { extension: string; mime: string; text: string } {
  const consensuses = bundle.consensuses.filter((record) => !sample || record.sample === sample);
  const records = bundle.records.filter((record) => !sample || record.sample === sample);
  switch (kind) {
    case "consensus-fasta": return { extension: "consensus.fasta", mime: "text/x-fasta", text: fasta(consensuses.map((record) => ({ id: record.id, sequence: record.sequence }))) };
    case "passed-consensus-fasta": return { extension: "passed-consensus.fasta", mime: "text/x-fasta", text: fasta(records.filter(passed).map((record) => ({ id: record.id, sequence: record.consensusNt }))) };
    case "rejected-consensus-fasta": return { extension: "rejected-consensus.fasta", mime: "text/x-fasta", text: fasta(records.filter((record) => !passed(record)).map((record) => ({ id: record.id, sequence: record.consensusNt }))) };
    case "trimmed-nt-fasta": return { extension: "trimmed-nt.fasta", mime: "text/x-fasta", text: fasta(records.filter((record) => record.functionalPass && record.trimmedNt).map((record) => ({ id: record.id, sequence: record.trimmedNt! }))) };
    case "trimmed-aa-fasta": return { extension: "trimmed-aa.fasta", mime: "text/x-fasta", text: fasta(records.filter((record) => record.functionalPass && record.trimmedAa).map((record) => ({ id: record.id, sequence: record.trimmedAa! }))) };
    case "family-csv": return { extension: "families.csv", mime: "text/csv", text: csv(
      ["sample", "UMI", "fs", "tags", "posterior_probability", "log_offspring_probability", "minag"],
      bundle.umiFamilies.filter((row) => !sample || row.sample === sample).map((row) => [row.sample, row.umi, row.familySize, row.disposition, row.posteriorProbability, row.logOffspringProbability, row.minimumAgreement]),
    ) };
    case "low-agreement-csv": return { extension: "low-agreement.csv", mime: "text/csv", text: csv(
      ["sample", "sequence_id", "UMI", "position_from_3prime", "agreement", "modal_read_base", "modal_run_length"],
      consensuses.flatMap((record) => record.lowAgreementSites.map((site) => [record.sample, record.id, record.umi, site.position, site.agreement, site.modalReadBase, site.modalRunLength])),
    ) };
    case "contamination-csv": return { extension: "contamination.csv", mime: "text/csv", text: csv(
      ["sample", "sequence_name", "nearest_nonself_variant", "nearest_nonself_distance", "flagged", "discarded", "suspect_only"],
      bundle.contamination.filter((row) => !sample || row.sample === sample).map((row) => [row.sample, row.sequenceId, row.nearestNonselfVariant, row.nearestNonselfDistance, row.flagged, row.discarded, row.suspectOnly]),
    ) };
    case "postproc-csv": return { extension: "postproc.csv", mime: "text/csv", text: csv(
      ["sample", "id", "UMI", "fs", "minag", "panel_score", "artefact_pass", "agreement_pass", "contamination_pass", "panel_pass", "functional_pass", "rejection_reasons"],
      records.map((row) => [row.sample, row.id, row.umi, row.familySize, row.minimumAgreement, row.panelScore, row.artefactPass, row.agreementPass, row.contaminationPass, row.panelPass, row.functionalPass, row.rejectionReasons.join(";")]),
    ) };
    case "apobec-csv": return { extension: "apobec.csv", mime: "text/csv", text: csv(
      ["sample", "id", "posterior_mean_GA_multiplier", "posterior_probability_GA_inflated", "posterior_mean_mutation_rate", "GA_mutations", "total_mutations"],
      records.filter((row) => row.apobec).map((row) => [row.sample, row.id, row.apobec!.posteriorMeanGaMultiplier, row.apobec!.posteriorGaInflated, row.apobec!.posteriorMeanMutationRate, row.apobec!.gaMutations, row.apobec!.totalMutations]),
    ) };
    case "nucleotide-alignment": {
      const selected = alignmentSample(bundle, sample), key = `${selected}/nucleotide`;
      return { extension: "nucleotide-alignment.fasta", mime: "text/x-fasta", text: bundle.alignmentEdits?.[key]?.fasta ?? bundle.alignments[key] ?? "" };
    }
    case "protein-alignment": {
      const selected = alignmentSample(bundle, sample), key = `${selected}/nucleotide`, edit = bundle.alignmentEdits?.[key];
      const nucleotide = edit?.fasta ?? bundle.alignments[key];
      return { extension: "protein-alignment.fasta", mime: "text/x-fasta", text: nucleotide ? translateAlignmentFasta(nucleotide, edit?.frameOffset ?? 0) : "" };
    }
    case "newick": { const selected = alignmentSample(bundle, sample), key = `${selected}/nucleotide`, edit = bundle.alignmentEdits?.[key]; return { extension: "tree.newick", mime: "text/plain", text: edit ? edit.treeNewick ?? "" : bundle.trees[key] ?? "" }; }
    case "log": return { extension: "log.txt", mime: "text/plain", text: bundle.log.join("\n") + "\n" };
  }
}

export function safeDatasetName(value: string) { return value.replace(/[^A-Za-z0-9_.-]+/g, "_") || "webporpid"; }
