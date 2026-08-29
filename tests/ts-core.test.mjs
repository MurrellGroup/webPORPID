import assert from "node:assert/strict";
import { mkdir, mkdtemp, readFile, rm } from "node:fs/promises";
import { join, resolve } from "node:path";
import test from "node:test";
import { pathToFileURL } from "node:url";
import { createElement } from "react";
import { renderToStaticMarkup } from "react-dom/server";
import { build } from "rolldown";
import { gunzipSync } from "fflate";

import { inspectAlignment, summarizeAlignmentChanges, translateAlignmentFasta, validateCorrectedAlignment } from "../src/alignment-utils.ts";
import { referenceDisplayColumns } from "../src/alignment-regions.ts";
import {
  ALIVIBE_BRIDGE_VERSION, ALIVIBE_SOURCE_REVISION, assertAlivibeInitialLoad, assertAlivibeRoundTripTarget,
  loadAlivibeNucleotideFasta, readAlivibeNucleotideFasta,
} from "../src/alivibe-roundtrip.ts";
import { decodeResult, encodeResult, exportComponent } from "../src/result-file.ts";
import { classifyContamination, classifyContaminationAsync, deduplicateContaminationCalls } from "../src/contamination.ts";
import { collapseAlignment } from "../src/collapse.ts";
import { buildExportArchive, SAMPLE_EXPORT_KINDS } from "../src/export-archive.ts";
import { nameMatchingSlot, referenceMappingRecords, referenceSlots } from "../src/input-mapping.ts";
import { mapParsimonyMutations } from "../src/phylo-mutations.ts";
import { postprocess } from "../src/postprocess.ts";
import { functionalFilterStats, porpidCallStats, sampleOverviewStats } from "../src/report-stats.ts";
import { adaptiveSpoolCutoff, PartitionStore, writeAllSync } from "../src/partition-store.ts";
import { runScalableMsa } from "../src/scalable-msa.ts";
import { parseSpoolRecordHeader, selectSpoolChunks } from "../src/spool-record.ts";
import { layoutTree, midpointRoot, parseNewick, rootOnOutgroup } from "../src/tree.ts";
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

function tarEntries(archive) {
  const bytes = gunzipSync(archive), decoder = new TextDecoder(), entries = new Map(); let offset = 0;
  while (offset + 512 <= bytes.byteLength && bytes.subarray(offset, offset + 512).some(Boolean)) {
    const header = bytes.subarray(offset, offset + 512), name = decoder.decode(header.subarray(0, 100)).replace(/\0.*$/, "");
    const sizeText = decoder.decode(header.subarray(124, 136)).replace(/\0.*$/, "").trim(), size = Number.parseInt(sizeText || "0", 8);
    offset += 512; entries.set(name, bytes.slice(offset, offset + size)); offset += Math.ceil(size / 512) * 512;
  }
  return entries;
}

function pairwiseLeafDistances(root) {
  const adjacency = new Map(), leaves = [];
  const visit = (node) => {
    if (!adjacency.has(node)) adjacency.set(node, []);
    if (!node.children.length) leaves.push(node);
    for (const child of node.children) {
      adjacency.get(node).push({ node: child, length: child.length });
      if (!adjacency.has(child)) adjacency.set(child, []);
      adjacency.get(child).push({ node, length: child.length }); visit(child);
    }
  };
  visit(root); const result = new Map();
  for (let left = 0; left < leaves.length; left += 1) {
    const distances = new Map([[leaves[left], 0]]), stack = [{ node: leaves[left], parent: undefined }];
    while (stack.length) {
      const { node, parent } = stack.pop();
      for (const edge of adjacency.get(node) ?? []) if (edge.node !== parent) {
        distances.set(edge.node, distances.get(node) + edge.length); stack.push({ node: edge.node, parent: node });
      }
    }
    for (let right = left + 1; right < leaves.length; right += 1)
      result.set([leaves[left].name, leaves[right].name].sort().join("|"), distances.get(leaves[right]));
  }
  return result;
}

function rootToLeafDistances(root) {
  const output = new Map(), stack = [{ node: root, distance: 0 }];
  while (stack.length) {
    const { node, distance } = stack.pop();
    if (!node.children.length) output.set(node.name, distance);
    else node.children.forEach((child) => stack.push({ node: child, distance: distance + child.length }));
  }
  return output;
}

class MemoryScratchWritable {
  constructor(file) { this.file = file; this.chunks = []; this.finished = false; }
  async write(bytes) { assert(!this.finished); this.file.writeCalls++; this.chunks.push(bytes.slice()); }
  async close() {
    assert(!this.finished); this.finished = true;
    const length = this.chunks.reduce((sum, chunk) => sum + chunk.byteLength, 0), output = new Uint8Array(length); let offset = 0;
    for (const chunk of this.chunks) { output.set(chunk, offset); offset += chunk.byteLength; } this.file.bytes = output;
  }
  async abort() { this.finished = true; this.chunks = []; }
}

class MemoryScratchFile {
  bytes = new Uint8Array();
  writeCalls = 0;
  async createWritable() { return new MemoryScratchWritable(this); }
  async getFile() { return new Blob([this.bytes]); }
}

class MemoryScratchDirectory {
  constructor(name) { this.name = name; this.directories = new Map(); this.files = new Map(); this.removed = []; }
  async getDirectoryHandle(name, options = {}) {
    let directory = this.directories.get(name);
    if (!directory && options.create) { directory = new MemoryScratchDirectory(name); this.directories.set(name, directory); }
    if (!directory) throw new Error("directory missing"); return directory;
  }
  async getFileHandle(name, options = {}) {
    let file = this.files.get(name);
    if (!file && options.create) { file = new MemoryScratchFile(); this.files.set(name, file); }
    if (!file) throw new Error("file missing"); return file;
  }
  async removeEntry(name) { if (!this.directories.delete(name) && !this.files.delete(name)) throw new Error("entry missing"); this.removed.push(name); }
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

test("external scratch mode streams complete no-downsampling partitions and removes its temporary directory", async () => {
  const root = new MemoryScratchDirectory("chosen-disk"), maximumHash = (1n << 64n) - 1n;
  const store = await PartitionStore.create(4, { sampleCount: 1, maximumReadsPerSample: 0, externalDirectory: root });
  try {
    assert.equal(store.mode, "external-directory");
    for (let index = 0; index < 40; index++) {
      const record = spoolRecord(0, (BigInt(index + 1) * 11400714819323198485n) & maximumHash, "AACCGGTT", `external-${index}`);
      await store.appendFrames(routedSpoolRecord(index % 4, record));
    }
    await store.seal();
    const temporary = [...root.directories.values()][0];
    assert.equal([...temporary.files.values()].reduce((sum, file) => sum + file.writeCalls, 0), 4, "small frames should be coalesced into one sequential write per partition");
    const selected = await Promise.all(Array.from({ length: 4 }, (_, partition) => store.readSelected(partition, [maximumHash])));
    assert.equal(selected.reduce((sum, bytes) => sum + countSpoolRecords(bytes), 0), 40);
    assert.equal(store.statistics().bypassedRecords, 0); assert(store.statistics().currentBytes > 0);
  } finally { await store.close(); }
  assert.equal(root.directories.size, 0); assert.equal(root.removed.length, 1);
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

test("modal-tip and midpoint rooting preserve topology, tips, and every patristic distance", () => {
  const original = parseNewick("(((a:1,b:2):3,c:4):2,(d:5,e:1):6);");
  const expected = pairwiseLeafDistances(original);
  const modalRooted = rootOnOutgroup(original, "a"), midpointRooted = midpointRoot(original);
  for (const rerooted of [modalRooted, midpointRooted]) {
    const observed = pairwiseLeafDistances(rerooted);
    assert.deepEqual([...observed.keys()].sort(), [...expected.keys()].sort());
    for (const [pair, distance] of expected)
      assert(Math.abs(observed.get(pair) - distance) < 1e-12, `${pair} changed after rerooting`);
  }
  assert.equal(modalRooted.children.find((child) => child.name === "a")?.length, 0, "the modal-tip edge must be exactly zero");
  assert.equal(modalRooted.children.find((child) => child.name !== "a")?.length, 1, "the opposite root edge retains the complete pendant length");
  const midpointDepths = rootToLeafDistances(midpointRooted);
  assert.equal(midpointDepths.get("b"), 9); assert.equal(midpointDepths.get("d"), 9);
  assert.deepEqual([...rootToLeafDistances(rootOnOutgroup(parseNewick("(only:0);"), "only")).keys()], ["only"], "single-tip rerooting must not create a phantom leaf");
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
  assert.deepEqual(biologicalEdit.rowChanges.map((row) => [row.name, row.substitutedNucleotides]), [["a", 1]]);
  assert.equal(biologicalEdit.addedRows.length, 1);
  assert.equal(biologicalEdit.removedRows.length, 1);
  const shiftedInsertion = validateCorrectedAlignment(">a\nAACCGG\n", ">a\nAACTCGG\n");
  assert.equal(shiftedInsertion.insertedNucleotides, 1); assert.equal(shiftedInsertion.substitutedNucleotides, 0);
  const removedWithoutReordering = summarizeAlignmentChanges(">a\nAAAA\n>b\nAAAA\n>c\nAAAA\n", ">a\nAAAA\n>c\nAAAA\n");
  assert.equal(removedWithoutReordering.rowOrderChanged, false);
  const reordered = summarizeAlignmentChanges(">a\nAAAA\n>b\nAAAA\n>c\nAAAA\n", ">c\nAAAA\n>a\nAAAA\n>b\nAAAA\n");
  assert.equal(reordered.rowOrderChanged, true);
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
    const { AlignmentTreeViewer, abundanceBubbleRadius, chooseModalRootTip, treeAlignmentSvg } = await import(pathToFileURL(join(directory, "viewer.mjs")).href);
    const bundle = resultBundle(), family = bundle.umiFamilies[0], consensus = bundle.consensuses[0], record = bundle.records[0];
    const markup = [
      createElement(charts.UmiDecisionPlot, { families: [family], artefactCutoff: 1, agreementThreshold: .6, outlierQuantile: .99 }),
      createElement(charts.ArtefactDecisionPlot, { families: [family], artefactCutoff: 1, agreementThreshold: .6, outlierQuantile: .99, artefactFraction: .25 }),
      createElement(charts.AgreementPositionPlot, { consensuses: [consensus], threshold: .6 }),
      createElement(charts.MdsApobecPlot, { records: [record] }),
      createElement(charts.DinucleotideHeatmaps, { families: [family] }),
      createElement(AlignmentTreeViewer, { fasta: ">a\nATGTAA\n>b\nATGTAG\n", newick: "(a:0.001,b:0.002);", alphabet: "nt", collapsed: true,
        leafMetadata: { a: { familyCount: 4 }, b: { familyCount: 1 } } }),
    ].map((component) => renderToStaticMarkup(component)).join("\n");
    for (const label of ["UMI family size", "artefact cutoff", "Minimum-agreement positions", "multidimensional scaling", "dinucleotide frequencies", "Phylogram coordinated"])
      assert.match(markup, new RegExp(label, "i"));
    assert.match(markup, /<canvas/); assert.match(markup, /substitutions\/site/); assert.match(markup, /Export SVG/); assert.match(markup, /Bubble area/);
    assert.match(markup, /Tree \+ alignment SVG/); assert.match(markup, /Apply regions/); assert.match(markup, /Amino acid/);
    assert.match(markup, /Modal node/); assert.match(markup, /Midpoint root/); assert.match(markup, /exactly zero-length edge/);
    const weightedModal = chooseModalRootTip([{ name: "a", sequence: "AAA" }, { name: "b", sequence: "CCC" }], ["a", "b"], { a: { familyCount: 1 }, b: { familyCount: 4 } });
    assert.deepEqual(weightedModal, { treeName: "b", sequenceName: "b", modalSequence: "CCC", representedFamilies: 4, matchingTips: 1 });
    const uncollapsedModal = chooseModalRootTip([{ name: "first", sequence: "AAA" }, { name: "second", sequence: "CCC" }, { name: "third", sequence: "AAA" }], ["first", "second", "third"]);
    assert.deepEqual(uncollapsedModal, { treeName: "first", sequenceName: "first", modalSequence: "AAA", representedFamilies: 2, matchingTips: 2 });
    const radiusOne = abundanceBubbleRadius(1, 25), radiusFour = abundanceBubbleRadius(4, 25), radiusMillion = abundanceBubbleRadius(1_000_000, 25);
    assert(Math.abs((Math.PI * radiusFour ** 2) / (Math.PI * radiusOne ** 2) - 4) < 1e-12, "bubble area must be exactly linear in UMI-family count");
    assert(radiusMillion > 2_000, "large family counts must not encounter a display cap");
    const combined = treeAlignmentSvg({ root: parseNewick("(a:0.001,b:0.002);"), treeWidth: 500, rowHeight: 20, cellWidth: 11,
      layoutMode: "phylogram", showNames: true, leafMetadata: { a: { familyCount: 4 }, b: { familyCount: 1 } }, leafLabels: { a: "a", b: "b" },
      sequencesByTip: { a: "ATGTAA", b: "ATGTAG" }, columns: [0, 1, 2, 3, 4, 5], labels: ["1", "2", "3", "4", "5", "6"], modal: "ATGTAA",
      alphabet: "nt", highlighter: false, bubbleAreaPerFamily: 25, showAbundanceScale: true, colorByAgreement: false, mutations: new Map(), mutationLimit: 2, tipLegend: [] });
    assert.match(combined, /Tree coordinated with aligned leaf sequences/); assert.match(combined, /<rect[^>]+fill="#78c679"/); assert.match(combined, />A<\/text>/); assert.match(combined, /no radius floor, cap, or dataset normalization/);
    const collapsedViewer = renderToStaticMarkup(createElement(AlignmentTreeViewer, { fasta: ">a\nATGTAA\n", newick: "(a:0.0);", collapsed: true,
      leafMetadata: { a: { familyCount: 3, minimumAgreement: .7 } } }));
    assert.doesNotMatch(collapsedViewer, /Tip color|family minimum agreement/i, "legacy agreement metadata must not appear as a collapsed-haplotype property");
    const familyViewer = renderToStaticMarkup(createElement(AlignmentTreeViewer, { fasta: ">a\nATGTAA\n", newick: "(a:0.0);",
      leafMetadata: { a: { familyCount: 1, minimumAgreement: .7 } } }));
    assert.match(familyViewer, /Family minimum agreement/i);
    const mismatched = renderToStaticMarkup(createElement(AlignmentTreeViewer, { fasta: ">a\nATGTAA\n", newick: "(a:0.0,b:0.0);", alphabet: "nt" }));
    assert.match(mismatched, /stale tree is hidden/i);
    const coordinates = charts.classicalMds(["AAAA", "AAAT", "AATA", "ATAA", "TAAA"]);
    const secondSpan = Math.max(...coordinates.map((row) => row[1])) - Math.min(...coordinates.map((row) => row[1]));
    assert(secondSpan > .1, "the second positive MDS axis must not be replaced by a negative eigenvector");
  } finally { await rm(directory, { recursive: true, force: true }); }
});

test("collapse counts retained UMI families without assigning family agreement to haplotypes", () => {
  const collapsed = collapseAlignment(">family_a\nAC-GT\n>family_b\nA-CGT\n>family_c\nAT-GT\n", "sample_1");
  assert.equal(collapsed.groups.length, 2);
  assert.equal(collapsed.groups[0].familyCount, 2);
  assert.equal("minimumAgreement" in collapsed.groups[0], false);
  assert.deepEqual(collapsed.groups[0].memberIds, ["family_a", "family_b"]);
  assert.match(collapsed.fasta, /^>family_a\nAC-GT\n>family_c\nAT-GT\n$/);
});

test("overview statistics expose family- and read-level rejection proportions", () => {
  const bundle = resultBundle();
  bundle.summaries[0].selectedReads = 10; bundle.summaries[0].downsampledReads = 2; bundle.summaries[0].demultiplexedReads = 12;
  bundle.quality.perSample[0] = 12; bundle.quality.demultiplexedReads = 12; bundle.quality.downsampledReads = 2;
  bundle.umiFamilies = [
    { ...bundle.umiFamilies[0], familySize: 6, disposition: "likely_real" },
    { ...bundle.umiFamilies[0], umi: "TTTTGGGG", familySize: 4, disposition: "heteroduplex" },
  ];
  const calls = porpidCallStats(bundle, "sample_1"), heteroduplex = calls.find((row) => row.key === "heteroduplex");
  assert.equal(heteroduplex.families, 1); assert.equal(heteroduplex.familyPercent, 50); assert.equal(heteroduplex.reads, 4); assert.equal(heteroduplex.readPercent, 40);
  const overview = sampleOverviewStats(bundle)[0]; assert(Math.abs(overview.downsampledPercent - 100 / 6) < 1e-12); assert.equal(overview.heteroduplexReadPercent, 40);
  const functional = functionalFilterStats(bundle, "sample_1"); assert(functional.some((row) => row.key === "frameshift"), "zero-count functional categories must remain visible");
  bundle.umiFamilies.push({ ...bundle.umiFamilies[0], umi: "REJECTED", familySize: 2, disposition: "BPB-rejects" });
  const bpb = porpidCallStats(bundle, "sample_1").find((row) => row.key === "BPB-rejects");
  assert.equal(bpb.families, 0, "the aggregate BPB read bucket must not be presented as one UMI family"); assert.equal(bpb.reads, 2);
});

test("contamination calls are unique per family with the primary decision taking precedence", () => {
  const calls = deduplicateContaminationCalls([
    { sample: "s", sequenceId: "x", nearestNonselfVariant: "suspect", nearestNonselfDistance: .01, flagged: true, discarded: false, suspectOnly: true },
    { sample: "s", sequenceId: "x", nearestNonselfVariant: "primary", nearestNonselfDistance: .012, flagged: true, discarded: true, suspectOnly: false },
    { sample: "s", sequenceId: "y", nearestNonselfVariant: "only", nearestNonselfDistance: .02, flagged: true, discarded: false, suspectOnly: true },
  ]);
  assert.equal(calls.length, 2); assert.equal(calls.find((row) => row.sequenceId === "x").nearestNonselfVariant, "primary");
  assert.equal(calls.find((row) => row.sequenceId === "x").discarded, true);
});

test("chunked contamination checks preserve synchronous decisions and emit live progress", async () => {
  const base = resultBundle(), parameters = { ...base.config.parameters, deterministicSeed: 1n, contaminationFilter: true,
    contaminationDistanceThreshold: 1 };
  const samples = ["s1", "s2"].map((name) => ({ name, cdnaPrimer: "AAAaaaNNNNNNNNTTT", secondStrandPrimer: "GGG",
    panel: "panel.fa", panelSequences: [] }));
  const config = { dataset: "progress", samples, contaminationPanel: "contam.fa", contaminationPanelSequences: [], parameters };
  const consensuses = [
    { ...base.consensuses[0], id: "s1-a", sample: "s1", sampleIndex: 0, sequence: "ACGTACGTACGT" },
    { ...base.consensuses[0], id: "s2-a", sample: "s2", sampleIndex: 1, sequence: "ACGTACGTACGA" },
  ];
  const updates = [], synchronous = classifyContamination(consensuses, config);
  const asynchronous = await classifyContaminationAsync(consensuses, config, undefined, (state) => updates.push(state));
  assert.deepEqual(asynchronous, synchronous); assert(updates.length >= 4); assert.equal(updates.at(-1).fraction, 1);
  assert(updates.some((state) => state.phase === "clustering")); assert(updates.some((state) => state.phase === "classification"));
});

test("functional alignment is codon-aware and clipped to the reference endpoints", async () => {
  const base = resultBundle(), reference = "ATGAAAAAATAA", query = "ATGAAAAAAGCTGCTTAA";
  const config = { dataset: "functional-trim", contaminationPanel: "contam.fa", contaminationPanelSequences: [],
    parameters: { ...base.config.parameters, deterministicSeed: 1n, functionalMatchThreshold: 0 },
    samples: [{ name: "sample_1", cdnaPrimer: "AAAaaaNNNNNNNNTTT", secondStrandPrimer: "GGG", panel: "panel.fa",
      panelSequences: [], functionalReference: "functional.fa", functionalReferenceSequence: { name: "ref", sequence: reference } }] };
  const consensus = [{ ...base.consensuses[0], sequence: query }];
  const runner = async (sequences, _signal, _iterations, mode) => mode === "amino-acid"
    ? ["MKK*--", "MKKAA*"] : [`${sequences[0]}------`, sequences[1]];
  const output = await postprocess(consensus, [], config, undefined, runner, 1, undefined, { collapse: false });
  const functional = inspectAlignment(output.alignments["sample_1/functional-nucleotide"], 1);
  const functionalReference = inspectAlignment(output.referenceAlignments["sample_1/functional-nucleotide"], 1);
  assert.equal(functional.columns, reference.length); assert.equal(functionalReference.columns, reference.length);
  assert.equal(functional.records[0].sequence, query.slice(0, reference.length));
  assert.equal(output.records[0].functionalPass, true);
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
    provenance: { webporpidVersion: "0.3.5", createdUtc: "2026-08-27T00:00:00.000Z", engine: "test", workers: 2,
      inputName: "reads.fastq.gz", inputSha256: "a".repeat(64), configSha256: "b".repeat(64), deterministicSeed: "1",
      upstreamBranch: "nanopore", upstreamCommit: "c".repeat(40) },
    config: { dataset: "test", samples: [{ name: "sample_1", cdnaPrimer: "AAAaaaNNNNNNNNTTT", secondStrandPrimer: "GGG", panel: "panel.fa", functionalReference: "functional.fa" }],
      contaminationPanel: "contam.fa", parameters },
    quality: { totalReads: 3, qualityReads: 3, badReads: 0, shortReads: 0, longReads: 0, primerRejects: 0, idRejects: 0,
      demultiplexedReads: 3, bpbRejects: 0, malformedRecords: 0, downsampledReads: 0, perSample: [3] },
    summaries: [{ sample: "sample_1", demultiplexedReads: 3, selectedReads: 3, downsampledReads: 0, observedUmis: 1, likelyRealUmis: 1, consensusSequences: 1,
      contaminationPassed: 1, postprocPassed: 1, collapsedSequences: 1, functionalPassed: 1, artefactCutoff: 1 }],
    umiFamilies: [{ sample: "sample_1", sampleIndex: 0, umi: "AACCGGTT", familySize: 3, mostLikelyParent: "AACCGGTT",
      posteriorProbability: 1, logOffspringProbability: Number.NEGATIVE_INFINITY, disposition: "likely_real", minimumAgreement: 0.67 }],
    consensuses: [{ id: "c1", sample: "sample_1", sampleIndex: 0, umi: "AACCGGTT", familySize: 3, minimumAgreement: 0.67,
      sequence: "ATGTAA", lowAgreementSites: [{ position: 2, agreement: 0.67, modalReadBase: "A", modalRunLength: 1 }] }],
    contamination: [], contaminationReferences: [{ name: "contam_ref", sequence: "ATGTAA" }],
    downstreamResources: { samples: [{ name: "sample_1", panelSequences: [{ name: "panel_ref", sequence: "ATGTAA" }], functionalReferenceSequence: { name: "functional_ref", sequence: "ATGTAA" } }] },
    records: [{ id: "c1", sample: "sample_1", umi: "AACCGGTT", familySize: 3, minimumAgreement: 0.67, consensusNt: "ATGTAA",
      alignedNt: "ATGTAA", trimmedNt: "ATGTAA", trimmedAa: "M*", panelScore: 0, artefactPass: true, agreementPass: true,
      contaminationPass: true, panelPass: true, functionalPass: true, rejectionReasons: [] }],
    alignments: { "sample_1/nucleotide": ">c1\nATGTAA\n", "sample_1/uncollapsed-nucleotide": ">c1\nATGTAA\n",
      "sample_1/protein": ">c1\nM*\n", "sample_1/uncollapsed-protein": ">c1\nM*\n",
      "sample_1/functional-nucleotide": ">c1\nATGTAA\n", "sample_1/functional-protein": ">c1\nM*\n" },
    referenceAlignments: { "sample_1/nucleotide": ">reference\nATGTAA\n", "sample_1/uncollapsed-nucleotide": ">reference\nATGTAA\n",
      "sample_1/functional-nucleotide": ">functional_reference\nATGTAA\n" },
    collapseGroups: { sample_1: [{ sample: "sample_1", representativeId: "c1", memberIds: ["c1"], familyCount: 1 }] },
    inputMappings: [{ slot: "panel.fa", role: "panel", expectedName: "panel.fa", uploadedName: "renamed.fasta", uploadedSize: 12 }],
    runOptions: { deferPhylogeny: false, deferContamination: false, deferPostprocessing: false, deferCollapse: false, spoolStorage: "external-directory" },
    optionalStages: Object.fromEntries(["contamination", "postprocessing", "collapse", "tree"].map((stage) => [stage, { state: "completed", detail: "test complete", updatedUtc: "2026-08-27T00:00:00.000Z" }])),
    trees: { "sample_1/nucleotide": "(c1:0.0);" }, log: ["complete"],
  };
}

test("result bundles round-trip, export, and reject structural corruption", () => {
  const bundle = resultBundle(), encoded = encodeResult(bundle), decoded = decodeResult(encoded);
  assert.deepEqual(decoded, bundle);
  assert.equal(exportComponent(decoded, "trimmed-aa-fasta", "sample_1").text, ">c1\nM*\n");
  assert.match(exportComponent(decoded, "collapse-csv", "sample_1").text, /representative_id,family_count,member_ids\nsample_1,c1,1,c1/);
  assert.equal(exportComponent(decoded, "uncollapsed-nucleotide-alignment", "sample_1").text, ">c1\nATGTAA\n");
  assert.equal(exportComponent(decoded, "functional-nucleotide-alignment", "sample_1").text, ">c1\nATGTAA\n");
  assert.equal(exportComponent(decoded, "functional-protein-alignment", "sample_1").text, ">c1\nM*\n");
  decoded.alignments["sample_1/protein"] = ">c1\nWRONG\n";
  assert.equal(exportComponent(decoded, "protein-alignment", "sample_1").text, ">c1\nM*\n", "protein export must translate the nucleotide alignment directly");
  const baselineFingerprint = inspectAlignment(decoded.alignments["sample_1/nucleotide"], 1).fingerprint;
  const editedFasta = ">c1\nATG---\n", editedFingerprint = inspectAlignment(editedFasta, 1).fingerprint;
  decoded.alignmentEdits = { "sample_1/nucleotide": { fasta: editedFasta, frameOffset: 0, baselineFingerprint, editedFingerprint,
    source: "Alivibe test", savedUtc: "2026-08-27T00:00:01.000Z", treeNewick: "(c1:0.0);",
    changes: summarizeAlignmentChanges(decoded.alignments["sample_1/nucleotide"], editedFasta) } };
  decoded.alignmentEditHistory = [{ alignmentKey: "sample_1/nucleotide", action: "alignment-edit", timestamp: "2026-08-27T00:00:01.000Z",
    source: "Alivibe test", details: ["3 deleted bases"], beforeFingerprint: baselineFingerprint, afterFingerprint: editedFingerprint }];
  const editedRoundTrip = decodeResult(encodeResult(decoded));
  assert.equal(exportComponent(editedRoundTrip, "nucleotide-alignment", "sample_1").text, editedFasta);
  assert.equal(exportComponent(editedRoundTrip, "protein-alignment", "sample_1").text, ">c1\nM-\n");
  assert.equal(exportComponent(editedRoundTrip, "newick", "sample_1").text, "(c1:0.0);");
  const wrongMagic = encoded.slice(); wrongMagic[0] ^= 0xff;
  assert.throws(() => decodeResult(wrongMagic), /not a webPORPID results file/);
  assert.throws(() => decodeResult(encoded.subarray(0, encoded.length - 5)), /(corrupt|truncated|too large)/);
  const inconsistent = structuredClone(bundle); inconsistent.records[0].consensusNt = "ATGTAG";
  assert.throws(() => encodeResult(inconsistent), /inconsistent consensus sequence/);
  const unknownStorage = structuredClone(bundle); unknownStorage.runOptions.spoolStorage = "cloud";
  assert.throws(() => encodeResult(unknownStorage), /spoolStorage is not recognized/);
  const wrongSelection = structuredClone(bundle); wrongSelection.summaries[0].selectedReads = 2; wrongSelection.summaries[0].downsampledReads = 1;
  assert.throws(() => encodeResult(wrongSelection), /selected-read count does not match/);
  const legacyCollapse = structuredClone(bundle); legacyCollapse.collapseGroups.sample_1[0].minimumAgreement = .67;
  assert.doesNotThrow(() => encodeResult(legacyCollapse), "results through 0.3.2 must remain loadable");
  legacyCollapse.collapseGroups.sample_1[0].minimumAgreement = .5;
  assert.throws(() => encodeResult(legacyCollapse), /legacy collapse group has inconsistent family-agreement metadata/);
  const partial = structuredClone(bundle); partial.contamination = []; partial.records = []; partial.alignments = {}; partial.referenceAlignments = {};
  partial.collapseGroups = {}; partial.trees = {}; delete partial.summaries[0].contaminationPassed; delete partial.summaries[0].postprocPassed;
  delete partial.summaries[0].collapsedSequences; delete partial.summaries[0].functionalPassed; delete partial.summaries[0].artefactCutoff;
  partial.optionalStages = Object.fromEntries(["contamination", "postprocessing", "collapse", "tree"].map((stage) => [stage,
    { state: "deferred", detail: "not computed", updatedUtc: "2026-08-27T00:00:00.000Z" }]));
  assert.deepEqual(decodeResult(encodeResult(partial)), partial, "consensus-only projects must remain valid and resumable");
  const impossibleStages = structuredClone(bundle); impossibleStages.optionalStages.contamination.state = "deferred";
  assert.throws(() => encodeResult(impossibleStages), /cannot be completed before its prerequisite/);
  const pre035 = structuredClone(bundle); delete pre035.summaries[0].selectedReads; delete pre035.summaries[0].downsampledReads;
  delete pre035.contaminationReferences; delete pre035.downstreamResources; delete pre035.optionalStages;
  assert.doesNotThrow(() => decodeResult(encodeResult(pre035)), "pre-0.3.5 results must remain loadable without new optional statistics and references");
});

test("export all is one gzip-compressed tar with every sample output in its sample directory", () => {
  const bundle = resultBundle(), entries = tarEntries(buildExportArchive(bundle)), decoder = new TextDecoder();
  assert.equal(entries.size, SAMPLE_EXPORT_KINDS.length + 11);
  assert(entries.has("README.txt")); assert(entries.has("test.webporpid")); assert(entries.has("run.log.txt"));
  assert.equal(decoder.decode(entries.get("sample_1/trimmed-aa.fasta")), exportComponent(bundle, "trimmed-aa-fasta", "sample_1").text);
  assert.equal(decoder.decode(entries.get("sample_1/families.csv")), exportComponent(bundle, "family-csv", "sample_1").text);
  assert.match(decoder.decode(entries.get("cross-sample-overview/parameters.csv")), /maxReadsPerSample/);
  assert.match(decoder.decode(entries.get("cross-sample-overview/sample-summary.csv")), /sample_1/);
  assert.deepEqual(decodeResult(entries.get("test.webporpid")), bundle);
});

test("live demultiplexing and Swig-style navigation-loss guards stay wired into the app", async () => {
  const [app, pipeline, styles] = await Promise.all([
    readFile("src/App.tsx", "utf8"), readFile("src/pipeline-worker.ts", "utf8"), readFile("src/styles.css", "utf8"),
  ]);
  assert.match(pipeline, /sampleAssignments:\s*sampleAssignments\(\)/);
  assert.match(app, /Live sample assignments/); assert.match(app, /last pipeline update/);
  assert.match(app, /addEventListener\("beforeunload", warnBeforeLeaving\)/);
  assert.match(app, /addEventListener\("popstate", interceptHistoryDeparture\)/);
  assert.match(app, /Leave anyway/); assert.match(app, /history\.go\(-2\)/);
  assert.match(app, /showDirectoryPicker === "function" \? "external-directory" : "automatic"/);
  assert.match(app, /Defer contamination checks/); assert.match(app, /Skip this step/); assert.match(pipeline, /type === "skip-stage"/);
  assert.match(styles, /overscroll-behavior-x:\s*none/);
  assert.match(app, /className="app-version"/); assert.match(app, /packageInformation\.version/);
  assert.match(styles, /\.app-version/);
});
