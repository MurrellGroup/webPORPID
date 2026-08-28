import assert from "node:assert/strict";
import { mkdir, mkdtemp, readFile, rm } from "node:fs/promises";
import { join, resolve } from "node:path";
import test from "node:test";
import { pathToFileURL } from "node:url";
import { createElement } from "react";
import { renderToStaticMarkup } from "react-dom/server";
import { build } from "rolldown";

import { inspectAlignment, translateAlignmentFasta, validateCorrectedAlignment } from "../src/alignment-utils.ts";
import { referenceDisplayColumns } from "../src/alignment-regions.ts";
import {
  ALIVIBE_BRIDGE_VERSION, ALIVIBE_SOURCE_REVISION, assertAlivibeInitialLoad, assertAlivibeRoundTripTarget,
  loadAlivibeNucleotideFasta, readAlivibeNucleotideFasta,
} from "../src/alivibe-roundtrip.ts";
import { decodeResult, encodeResult, exportComponent } from "../src/result-file.ts";
import { collapseAlignment } from "../src/collapse.ts";
import { nameMatchingSlot, referenceMappingRecords, referenceSlots } from "../src/input-mapping.ts";
import { mapParsimonyMutations } from "../src/phylo-mutations.ts";
import { adaptiveSpoolCutoff, PartitionStore, writeAllSync } from "../src/partition-store.ts";
import { runScalableMsa } from "../src/scalable-msa.ts";
import { parseSpoolRecordHeader, selectSpoolChunks } from "../src/spool-record.ts";
import { layoutTree, parseNewick } from "../src/tree.ts";
import { treeTipNames } from "../src/tree-names.ts";

function spoolRecord(sample, hash, umi = "AACCGGTT", name = "read", sequence = "ACGT") {
  const encoder = new TextEncoder(), umiBytes = encoder.encode(umi), nameBytes = encoder.encode(name), sequenceBytes = encoder.encode(sequence);
  const body = 20 + umiBytes.length + nameBytes.length + sequenceBytes.length * 2, output = new Uint8Array(body + 4);
  const view = new DataView(output.buffer); view.setUint32(0, body, true); view.setUint16(4, sample, true); view.setUint16(6, umiBytes.length, true);
  view.setUint32(8, nameBytes.length, true); view.setUint32(12, sequenceBytes.length, true); view.setBigUint64(16, hash, true);
  let offset = 24; output.set(umiBytes, offset); offset += umiBytes.length; output.set(nameBytes, offset); offset += nameBytes.length;
  output.set(sequenceBytes, offset); offset += sequenceBytes.length; output.set(new Uint8Array(sequenceBytes.length).fill(73), offset); return output;
}

function routedSpoolRecord(partition, record) {
  const output = new Uint8Array(5 + record.byteLength), view = new DataView(output.buffer);
  output[0] = partition; view.setUint32(1, record.byteLength, true); output.set(record, 5); return output;
}

function countSpoolRecords(bytes) {
  let count = 0, offset = 0;
  while (offset < bytes.byteLength) { const header = parseSpoolRecordHeader(bytes.subarray(offset)); offset += header.recordLength; count++; }
  assert.equal(offset, bytes.byteLength); return count;
}

test("spool cutoff filtering materializes only selected records", () => {
  const kept = spoolRecord(0, 50n), rejected = spoolRecord(0, 150n), anotherSample = spoolRecord(1, 999n);
  const selected = selectSpoolChunks([kept, rejected, anotherSample], [100n, 1000n]);
  assert.equal(selected.byteLength, kept.byteLength + anotherSample.byteLength);
  assert.equal(parseSpoolRecordHeader(selected).samplingHash, 50n);
  const corrupt = kept.slice(); new DataView(corrupt.buffer).setUint32(0, 21, true);
  assert.throws(() => selectSpoolChunks([corrupt], [100n]), /inconsistent lengths/);
});

test("OPFS short writes are retried without truncating a spool frame", () => {
  const output = new Uint8Array(23), writes = [];
  const handle = { write(bytes, { at }) {
    const length = Math.min(4, bytes.byteLength); output.set(bytes.subarray(0, length), at); writes.push({ at, length }); return length;
  } };
  const input = Uint8Array.from({ length: output.length }, (_, index) => index + 1);
  writeAllSync(handle, input, 0); assert.deepEqual(output, input); assert(writes.length > 1);
  assert.throws(() => writeAllSync({ write: () => 0 }, input, 0), /could not make progress/);
});

test("adaptive spool admission retains every final deterministic selection while bounding disk candidates", async () => {
  const maximumHash = (1n << 64n) - 1n, cap = 3;
  const hashes = [maximumHash, maximumHash - 1n, maximumHash / 2n, maximumHash / 8n, maximumHash / 16n,
    maximumHash * 3n / 4n, 1n, maximumHash / 3n, maximumHash / 32n, maximumHash - 2n];
  const finalCutoff = adaptiveSpoolCutoff(BigInt(hashes.length), cap);
  const expected = hashes.filter((hash) => hash <= finalCutoff).length;
  const store = await PartitionStore.create(2, { sampleCount: 1, maximumReadsPerSample: cap, compactionIntervalBytes: 1 });
  try {
    for (let index = 0; index < hashes.length; index++) {
      const record = spoolRecord(0, hashes[index], "AACCGGTT", `read-${index}`);
      await store.appendFrames(routedSpoolRecord(index % 2, record));
    }
    await store.compact([finalCutoff]);
    const selected = await Promise.all([store.readSelected(0, [finalCutoff]), store.readSelected(1, [finalCutoff])]);
    assert.equal(selected.reduce((sum, bytes) => sum + countSpoolRecords(bytes), 0), expected);
    const statistics = store.statistics();
    assert.equal(statistics.observedRecords, hashes.length); assert(statistics.bypassedRecords > 0);
    assert(statistics.reclaimedBytes > 0); assert(statistics.currentBytes < statistics.writtenBytes);
  } finally { await store.close(); }
});

test("original single-dataset PORPID YAML remains accepted", async () => {
  const base = resolve(".build/test-tmp"); await mkdir(base, { recursive: true }); const directory = await mkdtemp(join(base, "config-test-"));
  try {
    const output = join(directory, "config.mjs");
    await build({ input: resolve("src/config.ts"), output: { file: output, format: "es" } });
    const { parseConfigYaml, resolveReferenceFiles } = await import(pathToFileURL(output).href);
    const config = parseConfigYaml(`legacy_run:\n  sample_A:\n    cDNA_primer: CCGCTacgtaaNNNNNNNNGTCA\n    sec_str_primer: TAGG\n    panel: panels/panel.fa\n    af_override: 0.4\n`);
    assert.equal(config.dataset, "legacy_run"); assert.equal(config.samples[0].name, "sample_A");
    assert.equal(config.samples[0].artefactFractionOverride, 0.4); assert.equal(config.parameters.ldaThreshold, 0.995);
    config.parameters.contaminationFilter = false;
    const resolved = await resolveReferenceFiles(config, new Map([
      ["panels/panel.fa", async () => ">correct\nACGT\n"], ["archive/panel.fa", async () => ">wrong\nTGCA\n"],
    ]));
    assert.equal(resolved.samples[0].panelSequences[0].name, "correct", "an exact YAML path must win over a same-basename file in another slot");
  } finally { await rm(directory, { recursive: true, force: true }); }
});

test("single-tip Newick trees remain renderable", () => {
  const tree = parseNewick("FastTree output\n(only_tip:0.0);\n");
  assert.equal(tree.children.length, 1);
  assert.equal(tree.children[0].name, "only_tip");
  assert.equal(tree.children[0].length, 0);
});

test("FastTree tip identifiers remain unique after Newick sanitization", () => {
  assert.deepEqual(treeTipNames(["a b", "a=b", "a_b"]), ["a_b", "a_b__2_1", "a_b__3_1"]);
});

test("large MSA batching preserves order, residues, and rectangularity", async () => {
  const sequences = ["ACGT", ...Array.from({ length: 8_000 }, (_, index) => `A${"T".repeat(index % 4)}CGT`)];
  let calls = 0;
  const align = async (rows) => {
    calls++; const widths = rows.map((row) => row.length - 4), width = Math.max(...widths);
    return rows.map((row, index) => index === 0 && row === "ACGT" ? `A${"-".repeat(width)}CGT`
      : `A${"T".repeat(widths[index])}${"-".repeat(width - widths[index])}CGT`);
  };
  const result = await runScalableMsa(sequences, align, undefined, 3, "nucleotide");
  assert(calls > 1); assert.equal(result.length, sequences.length);
  assert(result.every((row) => row.length === result[0].length));
  result.forEach((row, index) => assert.equal(row.replaceAll("-", ""), sequences[index]));
});

test("aligned translation, permissive edit validation, and Swig tree scaling stay coordinated", () => {
  const aligned = ">a\nATG---GC-NNN\n>b\nATG---GCTNNN\n";
  assert.equal(translateAlignmentFasta(aligned), ">a\nM-XX\n>b\nM-AX\n");
  const corrected = validateCorrectedAlignment(">a\nATGCCC\n>b\nATGCCC\n", ">a\nATG-CC\n>b\nATG-CC\n");
  assert.equal(corrected.removedNucleotides, 2);
  const biologicalEdit = validateCorrectedAlignment(">a\nATGCCC\n>b\nATGCCC\n", ">a\nATGACC-\n>renamed\nATGCCCC\n");
  assert.equal(biologicalEdit.substitutedNucleotides, 1);
  assert.equal(biologicalEdit.addedRows.length, 1);
  assert.equal(biologicalEdit.removedRows.length, 1);
  const layout = layoutTree(parseNewick("(a:0.001,b:0.002);"), 500, 20, "phylogram", 20, 10);
  assert.equal(layout.leaves, 2); assert.equal(layout.maximumDistance, 0.002);
  assert(Math.max(...layout.nodes.map((node) => node.x)) > 400, "ordinary fractional branch lengths must fill the tree viewport");
  assert.equal(inspectAlignment(aligned).columns, 12);
});

test("bundled Alivibe bridge round-trips the exact nucleotide state and rejects stale editors", () => {
  const fasta = ">a\nATG---GCT\n>b\nATGAAA---\n";
  const snapshot = { version: ALIVIBE_BRIDGE_VERSION, sourceRevision: ALIVIBE_SOURCE_REVISION, alphabet: "NT", mode: "NT",
    frameOffset: 1, records: [{ name: "a", sequence: "ATG---GCT" }, { name: "b", sequence: "ATGAAA---" }], fasta };
  const bridge = { version: ALIVIBE_BRIDGE_VERSION, sourceRevision: ALIVIBE_SOURCE_REVISION,
    loadNucleotideFasta: (loaded, frame) => { assert.equal(loaded, fasta); assert.equal(frame, 1); return snapshot; },
    snapshotNucleotide: () => snapshot, installMsaRunner() {}, createMsaJob() { return null; },
    installFastTreeRunner() {}, runFastTree() { return Promise.resolve("(a:0,b:0);"); } };
  const editor = { swigAlivibeBridge: bridge };
  const loaded = loadAlivibeNucleotideFasta(editor, fasta, 1);
  assertAlivibeInitialLoad(fasta, loaded);
  assert.deepEqual(readAlivibeNucleotideFasta(editor), loaded);
  assert.doesNotThrow(() => assertAlivibeRoundTripTarget({ groupKey: "sample", alignmentFingerprint: "abc" }, { groupKey: "sample", alignmentFingerprint: "abc" }));
  assert.throws(() => assertAlivibeRoundTripTarget({ groupKey: "sample", alignmentFingerprint: "abc" }, { groupKey: "other", alignmentFingerprint: "abc" }), /lineage changed/);
  assert.throws(() => readAlivibeNucleotideFasta({ swigAlivibeBridge: { ...bridge, snapshotNucleotide: () => ({ ...snapshot, fasta: ">wrong\nATG\n" }) } }), /disagree/);
});

test("bundled Alivibe has no remote script dependency and exposes local MSA/FastTree hooks", async () => {
  const html = await readFile("public/tools/alivibe.html", "utf8");
  assert.doesNotMatch(html, /<script[^>]+src=["']https?:/i);
  for (const hook of ["installMsaRunner", "createMsaJob", "installFastTreeRunner", "runFastTree", "snapshotNucleotide"])
    assert.match(html, new RegExp(`${hook}\\s*\\(`));
  assert.match(html, new RegExp(ALIVIBE_SOURCE_REVISION));
});

test("report figures and the Swig-derived viewer render from a result payload", async () => {
  const base = resolve(".build/test-tmp"); await mkdir(base, { recursive: true }); const directory = await mkdtemp(join(base, "ui-render-"));
  try {
    await build({ input: { charts: resolve("src/components/charts.tsx"), viewer: resolve("src/components/alignment-tree-viewer.tsx") },
      external: /^react(?:\/|$)/, output: { dir: directory, format: "es", entryFileNames: "[name].mjs" } });
    const charts = await import(pathToFileURL(join(directory, "charts.mjs")).href);
    const { AlignmentTreeViewer } = await import(pathToFileURL(join(directory, "viewer.mjs")).href);
    const bundle = resultBundle(), family = bundle.umiFamilies[0], consensus = bundle.consensuses[0], record = bundle.records[0];
    const markup = [
      createElement(charts.UmiDecisionPlot, { families: [family], artefactCutoff: 1, agreementThreshold: .6, outlierQuantile: .99 }),
      createElement(charts.ArtefactDecisionPlot, { families: [family], artefactCutoff: 1, agreementThreshold: .6, outlierQuantile: .99, artefactFraction: .25 }),
      createElement(charts.AgreementPositionPlot, { consensuses: [consensus], threshold: .6 }),
      createElement(charts.MdsApobecPlot, { records: [record] }),
      createElement(charts.DinucleotideHeatmaps, { families: [family] }),
      createElement(AlignmentTreeViewer, { fasta: ">a\nATGTAA\n>b\nATGTAG\n", newick: "(a:0.001,b:0.002);", alphabet: "nt", collapsed: true,
        leafMetadata: { a: { familyCount: 4, minimumAgreement: .7 }, b: { familyCount: 1, minimumAgreement: .9 } } }),
    ].map((component) => renderToStaticMarkup(component)).join("\n");
    for (const label of ["UMI family size", "artefact cutoff", "Minimum-agreement positions", "multidimensional scaling", "dinucleotide frequencies", "Phylogram coordinated"])
      assert.match(markup, new RegExp(label, "i"));
    assert.match(markup, /<canvas/); assert.match(markup, /substitutions\/site/); assert.match(markup, /Export SVG/); assert.match(markup, /Bubble size/);
    assert.match(markup, /r="5\.6"/, "bubble radius must be proportional to the square root of family count");
    const mismatched = renderToStaticMarkup(createElement(AlignmentTreeViewer, { fasta: ">a\nATGTAA\n", newick: "(a:0.0,b:0.0);", alphabet: "nt" }));
    assert.match(mismatched, /stale tree is hidden/i);
    const coordinates = charts.classicalMds(["AAAA", "AAAT", "AATA", "ATAA", "TAAA"]);
    const secondSpan = Math.max(...coordinates.map((row) => row[1])) - Math.min(...coordinates.map((row) => row[1]));
    assert(secondSpan > .1, "the second positive MDS axis must not be replaced by a negative eigenvector");
  } finally { await rm(directory, { recursive: true, force: true }); }
});

test("collapse counts retained UMI families, not reads", () => {
  const records = [
    { id: "family_a", minimumAgreement: .91, familySize: 100 },
    { id: "family_b", minimumAgreement: .72, familySize: 2 },
    { id: "family_c", minimumAgreement: .88, familySize: 40 },
  ];
  const collapsed = collapseAlignment(">family_a\nAC-GT\n>family_b\nA-CGT\n>family_c\nAT-GT\n", "sample_1", records);
  assert.equal(collapsed.groups.length, 2);
  assert.equal(collapsed.groups[0].familyCount, 2);
  assert.equal(collapsed.groups[0].minimumAgreement, .72);
  assert.deepEqual(collapsed.groups[0].memberIds, ["family_a", "family_b"]);
  assert.match(collapsed.fasta, /^>family_a\nAC-GT\n>family_c\nAT-GT\n$/);
});

test("reference regions and mutation mapping follow the active alignment", () => {
  assert.deepEqual(referenceDisplayColumns("ATG---GCTAAC", "1;3;7-9", "nt", "nt", 0, 12), [0, 2, 9, 10, 11]);
  assert.deepEqual(referenceDisplayColumns("ATG---GCTAAC", "1;3", "aa", "nt", 0, 12), [0, 1, 2, 9, 10, 11]);
  const tree = parseNewick("((a:1,b:1):1,c:1);");
  const first = mapParsimonyMutations(tree, new Map([["a", "AAAA"], ["b", "AATA"], ["c", "TTTA"]]));
  const second = mapParsimonyMutations(tree, new Map([["a", "AAAA"], ["b", "AAAA"], ["c", "TTTA"]]));
  const signature = (mapping) => [...mapping.mutationsByClade].flatMap(([clade, rows]) => rows.map((row) => `${clade}:${row.column}:${row.from}>${row.to}`)).sort();
  assert.notDeepEqual(signature(first), signature(second), "editing a base must refresh the reconstructed branch mutations");
});

test("YAML-derived slots preserve explicit renamed-file mappings", () => {
  const config = resultBundle().config, slots = referenceSlots(config);
  assert.equal(nameMatchingSlot("panel.fa", slots)?.expectedName, "panel.fa");
  const renamed = new File([">p\nACGT\n"], "renamed-upload.fasta", { lastModified: 1 });
  const records = referenceMappingRecords(slots, { [slots[0].id]: renamed });
  assert.equal(records[0].expectedName, slots[0].expectedName);
  assert.equal(records[0].uploadedName, "renamed-upload.fasta");
});

test("the packaged FastTree WebAssembly binary is double precision", async () => {
  const bytes = await readFile("public/biowasm/fasttree/fasttree.wasm");
  assert.match(bytes.toString("latin1"), /Double precision \(No SSE3\)/);
});

function resultBundle() {
  const parameters = { errorRate: 0.05, minLength: 20, maxLength: 300, primerTolerance: 0, primerWindow: 100, primerChop: 0,
    maxReadsPerSample: 100_000, familySizeThreshold: 1, ldaThreshold: 0.995, contaminationClusterThreshold: 0.015,
    contaminationProportionThreshold: 0.2, contaminationDistanceThreshold: 0.015, contaminationFilter: false,
    agreementThreshold: 0.6, artefactFraction: 0.25, outlierQuantile: 0.99, panelThreshold: 50,
    functionalMatchThreshold: 0.7, spoolPartitions: 8, deterministicSeed: "1" };
  return {
    schema: "webporpid-results/1",
    provenance: { webporpidVersion: "0.3.0", createdUtc: "2026-08-27T00:00:00.000Z", engine: "test", workers: 2,
      inputName: "reads.fastq.gz", inputSha256: "a".repeat(64), configSha256: "b".repeat(64), deterministicSeed: "1",
      upstreamBranch: "nanopore", upstreamCommit: "c".repeat(40) },
    config: { dataset: "test", samples: [{ name: "sample_1", cdnaPrimer: "AAAaaaNNNNNNNNTTT", secondStrandPrimer: "GGG", panel: "panel.fa" }],
      contaminationPanel: "contam.fa", parameters },
    quality: { totalReads: 3, qualityReads: 3, badReads: 0, shortReads: 0, longReads: 0, primerRejects: 0, idRejects: 0,
      demultiplexedReads: 3, bpbRejects: 0, malformedRecords: 0, downsampledReads: 0, perSample: [3] },
    summaries: [{ sample: "sample_1", demultiplexedReads: 3, observedUmis: 1, likelyRealUmis: 1, consensusSequences: 1,
      contaminationPassed: 1, postprocPassed: 1, collapsedSequences: 1, functionalPassed: 1, artefactCutoff: 1 }],
    umiFamilies: [{ sample: "sample_1", sampleIndex: 0, umi: "AACCGGTT", familySize: 3, mostLikelyParent: "AACCGGTT",
      posteriorProbability: 1, logOffspringProbability: Number.NEGATIVE_INFINITY, disposition: "likely_real", minimumAgreement: 0.67 }],
    consensuses: [{ id: "c1", sample: "sample_1", sampleIndex: 0, umi: "AACCGGTT", familySize: 3, minimumAgreement: 0.67,
      sequence: "ATGTAA", lowAgreementSites: [{ position: 2, agreement: 0.67, modalReadBase: "A", modalRunLength: 1 }] }],
    contamination: [],
    records: [{ id: "c1", sample: "sample_1", umi: "AACCGGTT", familySize: 3, minimumAgreement: 0.67, consensusNt: "ATGTAA",
      alignedNt: "ATGTAA", trimmedNt: "ATGTAA", trimmedAa: "M*", panelScore: 0, artefactPass: true, agreementPass: true,
      contaminationPass: true, panelPass: true, functionalPass: true, rejectionReasons: [] }],
    alignments: { "sample_1/nucleotide": ">c1\nATGTAA\n", "sample_1/uncollapsed-nucleotide": ">c1\nATGTAA\n",
      "sample_1/protein": ">c1\nM*\n", "sample_1/uncollapsed-protein": ">c1\nM*\n" },
    referenceAlignments: { "sample_1/nucleotide": ">reference\nATGTAA\n" },
    collapseGroups: { sample_1: [{ sample: "sample_1", representativeId: "c1", memberIds: ["c1"], familyCount: 1, minimumAgreement: .67 }] },
    inputMappings: [{ slot: "panel.fa", role: "panel", expectedName: "panel.fa", uploadedName: "renamed.fasta", uploadedSize: 12 }],
    runOptions: { deferPhylogeny: false }, trees: { "sample_1/nucleotide": "(c1:0.0);" }, log: ["complete"],
  };
}

test("result bundles round-trip, export, and reject structural corruption", () => {
  const bundle = resultBundle(), encoded = encodeResult(bundle), decoded = decodeResult(encoded);
  assert.deepEqual(decoded, bundle);
  assert.equal(exportComponent(decoded, "trimmed-aa-fasta", "sample_1").text, ">c1\nM*\n");
  assert.match(exportComponent(decoded, "collapse-csv", "sample_1").text, /c1,1,0\.67,c1/);
  assert.equal(exportComponent(decoded, "uncollapsed-nucleotide-alignment", "sample_1").text, ">c1\nATGTAA\n");
  decoded.alignments["sample_1/protein"] = ">c1\nWRONG\n";
  assert.equal(exportComponent(decoded, "protein-alignment", "sample_1").text, ">c1\nM*\n", "protein export must translate the nucleotide alignment directly");
  const baselineFingerprint = inspectAlignment(decoded.alignments["sample_1/nucleotide"], 1).fingerprint;
  const editedFasta = ">c1\nATG---\n", editedFingerprint = inspectAlignment(editedFasta, 1).fingerprint;
  decoded.alignmentEdits = { "sample_1/nucleotide": { fasta: editedFasta, frameOffset: 0, baselineFingerprint, editedFingerprint,
    source: "Alivibe test", savedUtc: "2026-08-27T00:00:01.000Z", treeNewick: "(c1:0.0);" } };
  const editedRoundTrip = decodeResult(encodeResult(decoded));
  assert.equal(exportComponent(editedRoundTrip, "nucleotide-alignment", "sample_1").text, editedFasta);
  assert.equal(exportComponent(editedRoundTrip, "protein-alignment", "sample_1").text, ">c1\nM-\n");
  assert.equal(exportComponent(editedRoundTrip, "newick", "sample_1").text, "(c1:0.0);");
  const wrongMagic = encoded.slice(); wrongMagic[0] ^= 0xff;
  assert.throws(() => decodeResult(wrongMagic), /not a webPORPID results file/);
  assert.throws(() => decodeResult(encoded.subarray(0, encoded.length - 5)), /(corrupt|truncated|too large)/);
  const inconsistent = structuredClone(bundle); inconsistent.records[0].consensusNt = "ATGTAG";
  assert.throws(() => encodeResult(inconsistent), /inconsistent consensus sequence/);
});
