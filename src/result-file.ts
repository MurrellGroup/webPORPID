import { decode, encode } from "@msgpack/msgpack";
import { gzipSync, gunzipSync } from "fflate";
import { inspectAlignment, summarizeAlignmentChanges, translateAlignmentFasta, validateCorrectedAlignment } from "./alignment-utils.ts";
import { deduplicateContaminationCalls } from "./contamination.ts";
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
  if (parameters.panelFilterMode != null && parameters.panelFilterMode !== "mafft-batch" && parameters.panelFilterMode !== "independent-query")
    throw new Error("config.parameters.panelFilterMode must be mafft-batch or independent-query.");
  text(parameters.deterministicSeed, "config.parameters.deterministicSeed");
  const samples = array(config.samples, "config.samples").map((entry, index) => {
    const sample = object(entry, `config.samples[${index}]`); const name = text(sample.name, `config.samples[${index}].name`);
    optionalText(sample.donorId, `config.samples[${index}].donorId`);
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
    for (const key of ["demultiplexedReads", "observedUmis", "likelyRealUmis", "consensusSequences"])
      count(row[key], `summaries[${index}].${key}`);
    for (const key of ["contaminationPassed", "postprocPassed", "artefactCutoff"])
      if (row[key] != null) count(row[key], `summaries[${index}].${key}`);
    if (row.selectedReads != null) count(row.selectedReads, `summaries[${index}].selectedReads`);
    if (row.downsampledReads != null) count(row.downsampledReads, `summaries[${index}].downsampledReads`);
    if (row.selectedReads != null && row.downsampledReads != null
      && count(row.selectedReads, `summaries[${index}].selectedReads`) + count(row.downsampledReads, `summaries[${index}].downsampledReads`) !== count(row.demultiplexedReads, `summaries[${index}].demultiplexedReads`))
      throw new Error("A sample summary has inconsistent selected and subsampled read counts.");
    if (row.collapsedSequences != null) count(row.collapsedSequences, `summaries[${index}].collapsedSequences`);
    if (row.functionalPassed != null) count(row.functionalPassed, `summaries[${index}].functionalPassed`);
  });
  if (summarySamples.size !== samples.length) throw new Error("Result summaries are missing a configured sample.");

  const familyKeys = new Set<string>(), familyReadsBySample = new Map<string, number>();
  array(bundle.umiFamilies, "umiFamilies").forEach((entry, index) => {
    const row = object(entry, `umiFamilies[${index}]`), sample = text(row.sample, `umiFamilies[${index}].sample`), sampleIndex = count(row.sampleIndex, `umiFamilies[${index}].sampleIndex`);
    knownSample(sample, `umiFamilies[${index}]`); if (sampleIndices.get(sample) !== sampleIndex) throw new Error("A UMI family has an inconsistent sample index.");
    const umi = text(row.umi, `umiFamilies[${index}].umi`), familyKey = `${sampleIndex}\0${umi}`;
    if (familyKeys.has(familyKey)) throw new Error("UMI family identifiers must be unique within a sample."); familyKeys.add(familyKey);
    const familySize = count(row.familySize, `umiFamilies[${index}].familySize`); familyReadsBySample.set(sample, (familyReadsBySample.get(sample) ?? 0) + familySize);
    text(row.mostLikelyParent, `umiFamilies[${index}].mostLikelyParent`);
    numeric(row.posteriorProbability, `umiFamilies[${index}].posteriorProbability`); numeric(row.logOffspringProbability, `umiFamilies[${index}].logOffspringProbability`, false);
    if (!DISPOSITIONS.has(text(row.disposition, `umiFamilies[${index}].disposition`))) throw new Error("A UMI family has an unknown disposition.");
    optionalNumber(row.minimumAgreement, `umiFamilies[${index}].minimumAgreement`);
  });
  array(bundle.summaries, "summaries").forEach((entry, index) => {
    const row = object(entry, `summaries[${index}]`), sample = text(row.sample, `summaries[${index}].sample`);
    if (row.selectedReads != null && count(row.selectedReads, `summaries[${index}].selectedReads`) !== (familyReadsBySample.get(sample) ?? 0))
      throw new Error("A sample summary selected-read count does not match its stored family calls.");
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
  if (bundle.contaminationReferences != null) {
    array(bundle.contaminationReferences, "contaminationReferences").forEach((entry, index) => {
      const record = object(entry, `contaminationReferences[${index}]`); text(record.name, `contaminationReferences[${index}].name`);
      text(record.sequence, `contaminationReferences[${index}].sequence`);
    });
  }
  if (bundle.downstreamResources != null) {
    const resources = object(bundle.downstreamResources, "downstreamResources"), resourceSamples = new Set<string>();
    array(resources.samples, "downstreamResources.samples").forEach((rawSample, index) => {
      const resource = object(rawSample, `downstreamResources.samples[${index}]`), name = text(resource.name, "downstream resource sample");
      knownSample(name, `downstreamResources.samples[${index}]`); if (resourceSamples.has(name)) throw new Error("Downstream resource sample names must be unique.");
      resourceSamples.add(name);
      array(resource.panelSequences, "downstream panel sequences").forEach((rawRecord) => {
        const record = object(rawRecord, "downstream panel record"); text(record.name, "downstream panel name"); text(record.sequence, "downstream panel sequence");
      });
      if (resource.functionalReferenceSequence != null) {
        const record = object(resource.functionalReferenceSequence, "downstream functional reference");
        text(record.name, "downstream functional reference name"); text(record.sequence, "downstream functional reference sequence");
      }
    });
    if (resourceSamples.size !== samples.length) throw new Error("Downstream resources are missing a configured sample.");
  }

  const rawStageStatuses = bundle.optionalStages == null ? undefined : object(bundle.optionalStages, "optionalStages");
  const contaminationComplete = rawStageStatuses == null
    || object(rawStageStatuses.contamination, "optionalStages.contamination").state === "completed";
  const postprocessingComplete = rawStageStatuses == null
    || object(rawStageStatuses.postprocessing, "optionalStages.postprocessing").state === "completed";
  if (bundle.postprocessingContaminationMode != null) {
    const mode = text(bundle.postprocessingContaminationMode, "postprocessingContaminationMode");
    if (!["applied", "bypassed"].includes(mode)) throw new Error("postprocessingContaminationMode is not recognized.");
    if (!postprocessingComplete) throw new Error("An uncomputed post-processing stage cannot record a contamination mode.");
    if (mode === "applied" && !contaminationComplete) throw new Error("Post-processing cannot apply an uncomputed contamination stage.");
    if (mode === "applied" && parameters.contaminationFilter !== true) throw new Error("Post-processing cannot apply a disabled contamination filter.");
  }
  if (rawStageStatuses != null && !contaminationComplete && postprocessingComplete && bundle.postprocessingContaminationMode !== "bypassed")
    throw new Error("Post-processing completed without contamination decisions must be marked as contamination-bypassed.");

  const recordIds = new Set<string>(), recordMetadata = new Map<string, { sample: string; alignedNt?: string; minimumAgreement: number }>();
  array(bundle.records, "records").forEach((entry, index) => {
    const row = object(entry, `records[${index}]`), id = text(row.id, `records[${index}].id`);
    if (recordIds.has(id)) throw new Error("Post-processing identifiers must be unique."); recordIds.add(id);
    const sample = text(row.sample, "postproc sample"), source = consensusById.get(id); knownSample(sample, `records[${index}]`);
    if (!source || source.sample !== sample) throw new Error("A post-processing record references an unknown consensus or sample.");
    const umi = text(row.umi, "postproc UMI"), familySize = count(row.familySize, "postproc family size"), minimumAgreement = numeric(row.minimumAgreement, "postproc agreement");
    if (umi !== source.umi || familySize !== source.familySize || minimumAgreement !== source.minimumAgreement)
      throw new Error("A post-processing record has inconsistent consensus metadata.");
    const consensusNt = text(row.consensusNt, "postproc consensus"); if (consensusNt !== source.sequence) throw new Error("A post-processing record has inconsistent consensus sequence data.");
    const alignedNt = row.alignedNt == null ? undefined : text(row.alignedNt, "postproc aligned sequence");
    optionalText(row.trimmedNt, "postproc trimmed nucleotide"); optionalText(row.trimmedAa, "postproc trimmed protein");
    numeric(row.panelScore, "postproc panel score"); for (const key of ["artefactPass", "agreementPass", "contaminationPass", "panelPass"]) bool(row[key], `postproc ${key}`);
    optionalBool(row.functionalPass, "postproc functionalPass"); array(row.rejectionReasons, "postproc rejectionReasons").forEach((reason) => text(reason, "postproc rejection reason"));
    if (row.apobec != null) { const model = object(row.apobec, "postproc APOBEC"); for (const key of ["posteriorMeanGaMultiplier", "posteriorGaInflated", "posteriorMeanMutationRate", "gaMutations", "totalMutations"]) numeric(model[key], `APOBEC ${key}`); }
    recordMetadata.set(id, { sample, alignedNt, minimumAgreement });
  });
  if (postprocessingComplete && (recordIds.size !== consensusIds.size || [...consensusIds].some((id) => !recordIds.has(id))))
    throw new Error("Consensus and post-processing records are inconsistent.");
  if (!postprocessingComplete && recordIds.size) throw new Error("An uncomputed post-processing stage cannot contain partial post-processing records.");
  if (bundle.postprocessingContaminationMode === "bypassed"
    && array(bundle.records, "records").some((entry) => object(entry, "record").contaminationPass !== true))
    throw new Error("Bypassed contamination cannot reject a post-processing record.");

  for (const [label, entries] of [["alignments", object(bundle.alignments, "alignments")], ["trees", object(bundle.trees, "trees")]] as const)
    for (const [name, contents] of Object.entries(entries)) {
      text(name, `${label} name`); text(contents, `${label}.${name}`);
      const sample = name.split("/", 1)[0]; knownSample(sample, `${label}.${name}`);
    }
  if (bundle.referenceAlignments != null) for (const [name, contents] of Object.entries(object(bundle.referenceAlignments, "referenceAlignments"))) {
    const sample = name.split("/", 1)[0]; knownSample(sample, `referenceAlignments.${name}`);
    const reference = inspectAlignment(text(contents, `referenceAlignments.${name}`), 1);
    const nucleotide = object(bundle.alignments, "alignments")[name];
    if (nucleotide != null && inspectAlignment(text(nucleotide, `alignments.${name}`), 1).columns !== reference.columns)
      throw new Error("A stored reference row does not match its nucleotide alignment width.");
  }
  if (bundle.collapseGroups != null) for (const [sample, rawGroups] of Object.entries(object(bundle.collapseGroups, "collapseGroups"))) {
    knownSample(sample, `collapseGroups.${sample}`);
    const representatives = new Set<string>(), membersSeen = new Set<string>();
    const collapsed = inspectAlignment(text(object(bundle.alignments, "alignments")[`${sample}/nucleotide`], `alignments.${sample}/nucleotide`), 1);
    const uncollapsed = inspectAlignment(text(object(bundle.alignments, "alignments")[`${sample}/uncollapsed-nucleotide`], `alignments.${sample}/uncollapsed-nucleotide`), 1);
    const collapsedByName = new Map(collapsed.records.map((record) => [record.name, record.sequence.replaceAll("-", "")]));
    const uncollapsedByName = new Map(uncollapsed.records.map((record) => [record.name, record.sequence.replaceAll("-", "")]));
    const functionalSource = object(bundle.alignments, "alignments")[`${sample}/functional-nucleotide`];
    const functionalNames = functionalSource == null ? new Set<string>()
      : new Set(inspectAlignment(text(functionalSource, `alignments.${sample}/functional-nucleotide`), 1).records.map((record) => record.name));
    let hasCollapsedFunctionalCalls = false, collapsedFunctionalPasses = 0;
    array(rawGroups, `collapseGroups.${sample}`).forEach((rawGroup, index) => {
      const group = object(rawGroup, `collapseGroups.${sample}[${index}]`);
      if (text(group.sample, "collapse group sample") !== sample) throw new Error("A collapse group has an inconsistent sample.");
      const representative = text(group.representativeId, "collapse representative");
      if (representatives.has(representative)) throw new Error("Collapse representative identifiers must be unique.");
      representatives.add(representative);
      const members = array(group.memberIds, "collapse members").map((entry) => text(entry, "collapse member"));
      if (count(group.familyCount, "collapse family count") !== members.length) throw new Error("Collapse counts must count UMI families.");
      const representativeSequence = collapsedByName.get(representative);
      if (representativeSequence == null) throw new Error("A collapse representative is missing from the collapsed alignment.");
      const agreements: number[] = [];
      for (const member of members) {
        if (membersSeen.has(member)) throw new Error("A retained UMI family occurs in more than one collapse group.");
        membersSeen.add(member);
        if (uncollapsedByName.get(member) !== representativeSequence) throw new Error("A collapse group contains different nucleotide haplotypes.");
        const metadata = recordMetadata.get(member);
        if (!metadata || metadata.sample !== sample || metadata.alignedNt == null) throw new Error("A collapse member is not a retained UMI-family consensus.");
        agreements.push(metadata.minimumAgreement);
      }
      // Results through 0.3.2 stored a conservative, but ambiguously named,
      // minimum across the member-family minima. Keep accepting and checking
      // that legacy field without treating it as a haplotype property.
      if (group.minimumAgreement != null && numeric(group.minimumAgreement, "legacy collapse minimum agreement") !== Math.min(...agreements))
        throw new Error("A legacy collapse group has inconsistent family-agreement metadata.");
      optionalText(group.trimmedNt, "collapsed functional nucleotide"); optionalText(group.trimmedAa, "collapsed functional protein");
      if (group.functionalRejectionReasons != null) array(group.functionalRejectionReasons, "collapsed functional reasons")
        .forEach((reason) => text(reason, "collapsed functional reason"));
      if (group.functionalPass != null) {
        hasCollapsedFunctionalCalls = true; const functionalPass = bool(group.functionalPass, "collapsed functionalPass");
        if (functionalPass) collapsedFunctionalPasses++;
        if (functionalNames.has(representative) !== functionalPass)
          throw new Error("A collapsed functional decision is inconsistent with the functional alignment.");
      }
    });
    if (representatives.size !== collapsed.records.length || membersSeen.size !== uncollapsed.records.length)
      throw new Error("Collapse membership does not cover the stored nucleotide alignments.");
    const summary = array(bundle.summaries, "summaries").map((entry) => object(entry, "summary")).find((entry) => entry.sample === sample);
    if (summary?.collapsedSequences != null && count(summary.collapsedSequences, "summary collapsed count") !== representatives.size)
      throw new Error("A summary has an inconsistent collapsed haplotype count.");
    if (hasCollapsedFunctionalCalls && summary?.functionalPassed != null
      && count(summary.functionalPassed, "summary functional count") !== collapsedFunctionalPasses)
      throw new Error("A summary has an inconsistent collapsed functional-pass count.");
  }
  if (bundle.inputMappings != null) array(bundle.inputMappings, "inputMappings").forEach((rawMapping, index) => {
    const mapping = object(rawMapping, `inputMappings[${index}]`);
    text(mapping.slot, "input slot");
    if (!["reads", "configuration", "panel", "functional-reference", "contamination-panel"].includes(text(mapping.role, "input role")))
      throw new Error("An input mapping has an unknown role.");
    optionalText(mapping.expectedName, "expected filename");
    text(mapping.uploadedName, "uploaded filename"); count(mapping.uploadedSize, "uploaded size");
  });
  if (bundle.runOptions != null) {
    const options = object(bundle.runOptions, "runOptions"); bool(options.deferPhylogeny, "runOptions.deferPhylogeny");
    for (const key of ["deferContamination", "deferPostprocessing", "deferCollapse", "interactiveFiltering"])
      if (options[key] != null) bool(options[key], `runOptions.${key}`);
    if (options.spoolStorage != null && !["automatic", "external-directory"].includes(text(options.spoolStorage, "runOptions.spoolStorage")))
      throw new Error("runOptions.spoolStorage is not recognized.");
  }
  if (rawStageStatuses != null) {
    for (const stage of ["contamination", "postprocessing", "collapse", "tree"] as const) {
      const status = object(rawStageStatuses[stage], `optionalStages.${stage}`), state = text(status.state, `optionalStages.${stage}.state`);
      if (!["completed", "deferred", "skipped"].includes(state)) throw new Error(`optionalStages.${stage}.state is not recognized.`);
      text(status.detail, `optionalStages.${stage}.detail`); text(status.updatedUtc, `optionalStages.${stage}.updatedUtc`);
    }
    const state = (stage: "postprocessing" | "collapse" | "tree") => object(rawStageStatuses[stage], `optionalStages.${stage}`).state;
    if (state("postprocessing") !== "completed" && state("collapse") === "completed")
      throw new Error("optionalStages.collapse cannot be completed before post-processing.");
    if (state("collapse") !== "completed" && state("tree") === "completed")
      throw new Error("optionalStages.tree cannot be completed before haplotype collapse.");
  }
  if (bundle.alignmentEdits != null) for (const [name, rawEdit] of Object.entries(object(bundle.alignmentEdits, "alignmentEdits"))) {
    const sample = name.split("/", 1)[0]; knownSample(sample, `alignmentEdits.${name}`);
    if (name !== `${sample}/nucleotide` && name !== `${sample}/uncollapsed-nucleotide` && name !== `${sample}/functional-nucleotide`)
      throw new Error("Edited alignment keys must identify a stored nucleotide view.");
    const edit = object(rawEdit, `alignmentEdits.${name}`), fasta = text(edit.fasta, `alignmentEdits.${name}.fasta`);
    const frameOffset = count(edit.frameOffset, `alignmentEdits.${name}.frameOffset`);
    if (frameOffset > 2) throw new Error("Edited alignment frame offsets must be 0, 1, or 2.");
    const baselineFingerprint = text(edit.baselineFingerprint, `alignmentEdits.${name}.baselineFingerprint`);
    const editedFingerprint = text(edit.editedFingerprint, `alignmentEdits.${name}.editedFingerprint`);
    text(edit.source, `alignmentEdits.${name}.source`); text(edit.savedUtc, `alignmentEdits.${name}.savedUtc`);
    optionalText(edit.treeNewick, `alignmentEdits.${name}.treeNewick`);
    optionalText(edit.treeFingerprint, `alignmentEdits.${name}.treeFingerprint`);
    optionalBool(edit.treeStale, `alignmentEdits.${name}.treeStale`);
    if (edit.warnings != null) array(edit.warnings, `alignmentEdits.${name}.warnings`).forEach((warning) => text(warning, "alignment edit warning"));
    if (edit.changes != null) {
      const changes = object(edit.changes, `alignmentEdits.${name}.changes`);
      for (const field of ["rowsBefore", "rowsAfter", "columnsBefore", "columnsAfter", "removedNucleotides", "insertedNucleotides", "substitutedNucleotides"])
        count(changes[field], `alignmentEdits.${name}.changes.${field}`);
      bool(changes.rowOrderChanged, `alignmentEdits.${name}.changes.rowOrderChanged`);
      for (const field of ["rowOrderBefore", "rowOrderAfter", "removedRows", "addedRows", "changedRows"])
        array(changes[field], `alignmentEdits.${name}.changes.${field}`).forEach((value) => text(value, `alignmentEdits.${name}.changes.${field} entry`));
      array(changes.rowChanges, `alignmentEdits.${name}.changes.rowChanges`).forEach((rawRow, index) => {
        const row = object(rawRow, `alignmentEdits.${name}.changes.rowChanges[${index}]`); text(row.name, "alignment row change name");
        for (const field of ["removedNucleotides", "insertedNucleotides", "substitutedNucleotides"]) count(row[field], `alignment row change ${field}`);
        bool(row.gapPlacementChanged, "alignment row gap-placement flag");
      });
    }
    const storedAlignments = object(bundle.alignments, "alignments");
    const originalValue = storedAlignments[name] ?? (name === `${sample}/uncollapsed-nucleotide` ? storedAlignments[`${sample}/nucleotide`] : undefined);
    const original = text(originalValue, `alignments.${name}`);
    if (inspectAlignment(original, 1).fingerprint !== baselineFingerprint) throw new Error(`alignmentEdits.${name} has an inconsistent baseline fingerprint.`);
    if (inspectAlignment(fasta, 1).fingerprint !== editedFingerprint) throw new Error(`alignmentEdits.${name} has an inconsistent fingerprint.`);
    validateCorrectedAlignment(original, fasta);
    if (edit.changes != null) {
      const expected = summarizeAlignmentChanges(original, fasta), stored = edit.changes as Record<string, unknown>;
      for (const field of ["rowsBefore", "rowsAfter", "columnsBefore", "columnsAfter", "removedNucleotides", "insertedNucleotides", "substitutedNucleotides"] as const)
        if (stored[field] !== expected[field]) throw new Error(`alignmentEdits.${name}.changes.${field} is inconsistent with the stored alignments.`);
      if (stored.rowOrderChanged !== expected.rowOrderChanged) throw new Error(`alignmentEdits.${name}.changes.rowOrderChanged is inconsistent with the stored alignments.`);
      for (const field of ["rowOrderBefore", "rowOrderAfter", "removedRows", "addedRows", "changedRows"] as const)
        if (JSON.stringify(stored[field]) !== JSON.stringify(expected[field])) throw new Error(`alignmentEdits.${name}.changes.${field} is inconsistent with the stored alignments.`);
      if (JSON.stringify(stored.rowChanges) !== JSON.stringify(expected.rowChanges)) throw new Error(`alignmentEdits.${name}.changes.rowChanges is inconsistent with the stored alignments.`);
    }
  }
  if (bundle.alignmentEditHistory != null) array(bundle.alignmentEditHistory, "alignmentEditHistory").forEach((rawEntry, index) => {
    const entry = object(rawEntry, `alignmentEditHistory[${index}]`), key = text(entry.alignmentKey, `alignmentEditHistory[${index}].alignmentKey`);
    const sample = key.split("/", 1)[0]; knownSample(sample, `alignmentEditHistory[${index}]`);
    if (key !== `${sample}/nucleotide` && key !== `${sample}/uncollapsed-nucleotide` && key !== `${sample}/functional-nucleotide`)
      throw new Error("Alignment audit keys must identify a nucleotide view.");
    if (!["alignment-edit", "frame-change", "tree-recalculation", "edit-reset"].includes(text(entry.action, `alignmentEditHistory[${index}].action`)))
      throw new Error("Alignment audit action is not recognized.");
    text(entry.timestamp, `alignmentEditHistory[${index}].timestamp`); text(entry.source, `alignmentEditHistory[${index}].source`);
    array(entry.details, `alignmentEditHistory[${index}].details`).forEach((value) => text(value, "alignment audit detail"));
    optionalText(entry.beforeFingerprint, `alignmentEditHistory[${index}].beforeFingerprint`);
    optionalText(entry.afterFingerprint, `alignmentEditHistory[${index}].afterFingerprint`);
  });
  if (bundle.timings != null) array(bundle.timings, "timings").forEach((entry, index) => {
    const timing = object(entry, `timings[${index}]`);
    const stage = text(timing.stage, `timings[${index}].stage`);
    if (!["setup", "preprocessing", "umi", "consensus", "contamination", "postprocessing", "collapse", "tree", "analysis-total"].includes(stage))
      throw new Error(`timings[${index}] has an unknown stage.`);
    const seconds = numeric(timing.seconds, `timings[${index}].seconds`);
    if (seconds < 0) throw new Error(`timings[${index}].seconds must be non-negative.`);
    if (timing.workItems != null) count(timing.workItems, `timings[${index}].workItems`);
  });
  if (bundle.thresholdSelections != null) {
    const ids = new Set<string>();
    array(bundle.thresholdSelections, "thresholdSelections").forEach((rawEntry, index) => {
      const entry = object(rawEntry, `thresholdSelections[${index}]`), id = text(entry.id, `thresholdSelections[${index}].id`);
      if (ids.has(id)) throw new Error("Interactive threshold checkpoint identifiers must be unique."); ids.add(id);
      const phase = text(entry.phase, `thresholdSelections[${index}].phase`);
      if (!['umi', 'consensus-filters'].includes(phase)) throw new Error("Interactive threshold phase is not recognized.");
      text(entry.acceptedUtc, `thresholdSelections[${index}].acceptedUtc`);
      array(entry.changes, `thresholdSelections[${index}].changes`).forEach((change) => text(change, "threshold change"));
      const selectedParameters = object(entry.parameters, `thresholdSelections[${index}].parameters`);
      for (const key of ["ldaThreshold", "familySizeThreshold", "artefactFraction", "outlierQuantile", "agreementThreshold"])
        if (selectedParameters[key] != null) numeric(selectedParameters[key], `thresholdSelections[${index}].parameters.${key}`);
      const selectedSamples = new Set<string>();
      array(entry.samples, `thresholdSelections[${index}].samples`).forEach((rawSample, sampleIndex) => {
        const selected = object(rawSample, `thresholdSelections[${index}].samples[${sampleIndex}]`);
        const name = text(selected.sample, "threshold sample"); knownSample(name, "threshold sample");
        if (selectedSamples.has(name)) throw new Error("An interactive threshold checkpoint contains duplicate sample rows."); selectedSamples.add(name);
        for (const key of ["familySizeOverride", "artefactFractionOverride", "outlierQuantileOverride", "agreementOverride"])
          if (selected[key] != null) numeric(selected[key], `threshold sample ${key}`);
      });
    });
  }
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
  | "collapse-csv" | "nucleotide-alignment" | "protein-alignment" | "newick"
  | "uncollapsed-nucleotide-alignment" | "uncollapsed-protein-alignment" | "uncollapsed-newick"
  | "functional-nucleotide-alignment" | "functional-protein-alignment" | "functional-newick" | "log";

function alignmentSample(bundle: ResultBundle, sample?: string) {
  if (sample) return sample;
  if (bundle.config.samples.length === 1) return bundle.config.samples[0].name;
  throw new Error("Choose a sample when exporting a sample-specific alignment or tree.");
}

export function exportComponent(bundle: ResultBundle, kind: ExportKind, sample?: string): { extension: string; mime: string; text: string } {
  const consensuses = bundle.consensuses.filter((record) => !sample || record.sample === sample);
  const records = bundle.records.filter((record) => !sample || record.sample === sample);
  const collapsedFunctional = Object.entries(bundle.collapseGroups ?? {}).flatMap(([groupSample, groups]) =>
    (!sample || groupSample === sample) ? groups.filter((group) => group.functionalPass) : []);
  switch (kind) {
    case "consensus-fasta": return { extension: "consensus.fasta", mime: "text/x-fasta", text: fasta(consensuses.map((record) => ({ id: record.id, sequence: record.sequence }))) };
    case "passed-consensus-fasta": return { extension: "passed-consensus.fasta", mime: "text/x-fasta", text: fasta(records.filter(passed).map((record) => ({ id: record.id, sequence: record.consensusNt }))) };
    case "rejected-consensus-fasta": return { extension: "rejected-consensus.fasta", mime: "text/x-fasta", text: fasta(records.filter((record) => !passed(record)).map((record) => ({ id: record.id, sequence: record.consensusNt }))) };
    case "trimmed-nt-fasta": return { extension: "trimmed-nt.fasta", mime: "text/x-fasta", text: fasta(
      (collapsedFunctional.length ? collapsedFunctional.map((group) => ({ id: group.representativeId, sequence: group.trimmedNt ?? "" }))
        : records.filter((record) => record.functionalPass && record.trimmedNt).map((record) => ({ id: record.id, sequence: record.trimmedNt! })))
        .filter((row) => row.sequence)) };
    case "trimmed-aa-fasta": return { extension: "trimmed-aa.fasta", mime: "text/x-fasta", text: fasta(
      (collapsedFunctional.length ? collapsedFunctional.map((group) => ({ id: group.representativeId, sequence: group.trimmedAa ?? "" }))
        : records.filter((record) => record.functionalPass && record.trimmedAa).map((record) => ({ id: record.id, sequence: record.trimmedAa! })))
        .filter((row) => row.sequence)) };
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
      deduplicateContaminationCalls(bundle.contamination).filter((row) => !sample || row.sample === sample).map((row) => [row.sample, row.sequenceId, row.nearestNonselfVariant, row.nearestNonselfDistance, row.flagged, row.discarded, row.suspectOnly]),
    ) };
    case "postproc-csv": return { extension: "postproc.csv", mime: "text/csv", text: csv(
      ["sample", "id", "UMI", "fs", "minag", "panel_score", "artefact_pass", "agreement_pass", "contamination_pass", "panel_pass", "rejection_reasons"],
      records.map((row) => [row.sample, row.id, row.umi, row.familySize, row.minimumAgreement, row.panelScore, row.artefactPass, row.agreementPass, row.contaminationPass, row.panelPass, row.rejectionReasons.join(";")]),
    ) };
    case "apobec-csv": return { extension: "apobec.csv", mime: "text/csv", text: csv(
      ["sample", "id", "posterior_mean_GA_multiplier", "posterior_probability_GA_inflated", "posterior_mean_mutation_rate", "GA_mutations", "total_mutations"],
      records.filter((row) => row.apobec).map((row) => [row.sample, row.id, row.apobec!.posteriorMeanGaMultiplier, row.apobec!.posteriorGaInflated, row.apobec!.posteriorMeanMutationRate, row.apobec!.gaMutations, row.apobec!.totalMutations]),
    ) };
    case "collapse-csv": {
      const selected = alignmentSample(bundle, sample);
      return { extension: "collapsed-families.csv", mime: "text/csv", text: csv(
        ["sample", "representative_id", "family_count", "functional_pass", "functional_rejection_reasons", "member_ids"],
        (bundle.collapseGroups?.[selected] ?? []).map((group) => [selected, group.representativeId, group.familyCount,
          group.functionalPass, group.functionalRejectionReasons?.join(";"), group.memberIds.join(";")]),
      ) };
    }
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
    case "uncollapsed-nucleotide-alignment": {
      const selected = alignmentSample(bundle, sample), key = `${selected}/uncollapsed-nucleotide`, edit = bundle.alignmentEdits?.[key];
      return { extension: "uncollapsed-nucleotide-alignment.fasta", mime: "text/x-fasta", text: edit?.fasta ?? bundle.alignments[key] ?? bundle.alignments[`${selected}/nucleotide`] ?? "" };
    }
    case "uncollapsed-protein-alignment": {
      const selected = alignmentSample(bundle, sample), key = `${selected}/uncollapsed-nucleotide`, edit = bundle.alignmentEdits?.[key];
      const nucleotide = edit?.fasta ?? bundle.alignments[key] ?? bundle.alignments[`${selected}/nucleotide`];
      return { extension: "uncollapsed-protein-alignment.fasta", mime: "text/x-fasta", text: nucleotide ? translateAlignmentFasta(nucleotide, edit?.frameOffset ?? 0) : "" };
    }
    case "uncollapsed-newick": {
      const selected = alignmentSample(bundle, sample), key = `${selected}/uncollapsed-nucleotide`, edit = bundle.alignmentEdits?.[key];
      return { extension: "uncollapsed-tree.newick", mime: "text/plain", text: edit?.treeNewick ?? bundle.trees[key] ?? "" };
    }
    case "functional-nucleotide-alignment": {
      const selected = alignmentSample(bundle, sample), key = `${selected}/functional-nucleotide`, edit = bundle.alignmentEdits?.[key];
      return { extension: "functional-nucleotide-alignment.fasta", mime: "text/x-fasta", text: edit?.fasta ?? bundle.alignments[key] ?? "" };
    }
    case "functional-protein-alignment": {
      const selected = alignmentSample(bundle, sample), key = `${selected}/functional-nucleotide`, edit = bundle.alignmentEdits?.[key];
      const nucleotide = edit?.fasta ?? bundle.alignments[key];
      return { extension: "functional-protein-alignment.fasta", mime: "text/x-fasta",
        text: nucleotide ? translateAlignmentFasta(nucleotide, edit?.frameOffset ?? 0) : bundle.alignments[`${selected}/functional-protein`] ?? "" };
    }
    case "functional-newick": {
      const selected = alignmentSample(bundle, sample), key = `${selected}/functional-nucleotide`, edit = bundle.alignmentEdits?.[key];
      return { extension: "functional-tree.newick", mime: "text/plain", text: edit?.treeNewick ?? bundle.trees[key] ?? "" };
    }
    case "log": return { extension: "log.txt", mime: "text/plain", text: bundle.log.join("\n") + "\n" };
  }
}

export function safeDatasetName(value: string) { return value.replace(/[^A-Za-z0-9_.-]+/g, "_") || "webporpid"; }
