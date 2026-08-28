/// <reference lib="webworker" />

import { bytesToHex } from "@noble/hashes/utils.js";
import coreWasmUrl from "/webporpid.wasm?url";
import { classifyContamination } from "./contamination";
import { runAlivibeMsa } from "./alivibe-msa-runtime";
import { compileConfig, resultConfig } from "./config";
import { finishStreamingHash, createStreamingHash, streamFastq } from "./fastq-stream";
import { PartitionStore, type ExternalScratchDirectoryHandle } from "./partition-store";
import { postprocess } from "./postprocess";
import { runFastTreeIsolated } from "./biowasm";
import { treeTipNames } from "./tree-names";
import type { InputFileMapping, PipelineConfig, PipelineProgress, QualityStats, ResultBundle } from "./types";
import { CoreWorkerPool } from "./worker-pool";
import {
  decodeConsensusOutput, decodeFamilyCounts, decodeFamilyModel, makeCutoffs, makeCutoffValues, mergeFamilyCounts,
  mergeStats,
} from "./wasm-runtime";

type RunRequest = { type: "run"; file: File; config: PipelineConfig; workers: number; deferPhylogeny?: boolean; inputMappings?: InputFileMapping[];
  spoolStorage?: "automatic" | "external-directory"; scratchDirectory?: ExternalScratchDirectoryHandle };
type CancelRequest = { type: "cancel" };
let cancellation: AbortController | undefined;

function progress(value: PipelineProgress) { self.postMessage({ type: "progress", progress: value }); }
const now = () => new Date().toISOString();

function formatBytes(bytes: number) {
  const units = ["B", "KiB", "MiB", "GiB", "TiB"]; let value = Math.max(0, bytes), unit = 0;
  while (value >= 1024 && unit < units.length - 1) { value /= 1024; unit++; }
  return `${value.toFixed(unit === 0 ? 0 : value >= 10 ? 1 : 2)} ${units[unit]}`;
}

function transferableBuffer(bytes: Uint8Array) {
  return bytes.byteOffset === 0 && bytes.buffer instanceof ArrayBuffer && bytes.byteLength === bytes.buffer.byteLength
    ? bytes.buffer : bytes.slice().buffer;
}

async function compileCore() {
  const url = coreWasmUrl;
  try { return await WebAssembly.compileStreaming(fetch(url)); }
  catch { const response = await fetch(url); if (!response.ok) throw new Error("The webPORPID WebAssembly core could not be loaded."); return WebAssembly.compile(await response.arrayBuffer()); }
}

function starTree(fasta: string) {
  const names = treeTipNames([...fasta.matchAll(/^>(.*)$/gm)].map((match) => match[1]));
  return `(${names.map((name) => `${name}:0.0`).join(",")});`;
}

async function run(request: RunRequest, signal: AbortSignal): Promise<ResultBundle> {
  const runStarted = performance.now(), timings: NonNullable<ResultBundle["timings"]> = [];
  const workers = Math.max(1, Math.floor(request.workers)), compiledConfig = compileConfig(request.config);
  if (request.spoolStorage === "external-directory" && !request.scratchDirectory)
    throw new Error("External scratch storage was selected, but no writable directory handle was provided.");
  const configForHash = compiledConfig.slice().buffer as ArrayBuffer;
  const [module, configHashBytes] = await Promise.all([compileCore(), crypto.subtle.digest("SHA-256", configForHash)]);
  const pool = await CoreWorkerPool.create(workers, module, compiledConfig);
  let store: PartitionStore;
  try {
    store = await PartitionStore.create(request.config.parameters.spoolPartitions, {
      sampleCount: request.config.samples.length,
      maximumReadsPerSample: request.config.parameters.maxReadsPerSample,
      requestPersistence: request.file.size >= 256 * 1024 * 1024,
      externalDirectory: request.spoolStorage === "external-directory" ? request.scratchDirectory : undefined,
    });
  } catch (cause) { pool.close(); throw cause; }
  const storageLabel = store.mode === "external-directory" ? "user-selected external scratch directory"
    : store.mode === "opfs" ? "browser OPFS" : "bounded memory fallback";
  const inputHash = createStreamingHash(), log = [`${now()} webPORPID 0.3.5 started`,
    `${now()} execution: ${workers} WASM workers; ${storageLabel} ${request.config.parameters.maxReadsPerSample > 0 ? "adaptive selected-read" : "all-read"} partition spool`,
    `${now()} parameters: error_rate=${request.config.parameters.errorRate}, lengths=(${request.config.parameters.minLength},${request.config.parameters.maxLength}), lda=${request.config.parameters.ldaThreshold}`];
  if (store.storage.quotaBytes != null) log.push(`${now()} browser storage: ${formatBytes(store.storage.usageBytes ?? 0)} used of ${formatBytes(store.storage.quotaBytes)} quota; persistence=${store.storage.persisted == null ? "unknown" : store.storage.persisted ? "granted" : "not granted"}`);
  for (const mapping of request.inputMappings ?? []) {
    log.push(`${now()} input mapping: ${mapping.role} slot ${mapping.slot}${mapping.expectedName ? ` (${mapping.expectedName})` : ""} <- ${mapping.uploadedName} (${mapping.uploadedSize} bytes)`);
  }
  const storeTiming = (stage: Exclude<(typeof timings)[number]["stage"], "analysis-total">, seconds: number, workItems?: number) => {
    const entry: (typeof timings)[number] = { stage, seconds };
    if (workItems != null) entry.workItems = workItems;
    timings.push(entry); log.push(`${now()} timing ${stage}: ${entry.seconds.toFixed(6)} s${workItems == null ? "" : `; ${workItems} work items`}`);
  };
  const recordTiming = (stage: Exclude<(typeof timings)[number]["stage"], "analysis-total">, started: number, workItems?: number) => {
    storeTiming(stage, (performance.now() - started) / 1000, workItems);
    return performance.now();
  };
  recordTiming("setup", runStarted);
  try {
    let stageStarted = performance.now();
    const pending = new Set<Promise<void>>(); let batches = 0, streamedRecords = 0, streamedFraction = 0, lastDemuxUpdate = 0;
    const sampleAssignments = () => {
      const observed = store.statistics().observedPerSample;
      return request.config.samples.map((sample, index) => ({ sample: sample.name, reads: Number(observed[index] ?? 0n) }));
    };
    const emitDemuxProgress = (force = false, detail?: string) => {
      const timestamp = performance.now(); if (!force && timestamp - lastDemuxUpdate < 100) return;
      lastDemuxUpdate = timestamp; const assigned = store.statistics().observedRecords;
      progress({ stage: "preprocessing", fraction: streamedFraction,
        detail: detail ?? `${streamedRecords.toLocaleString()} reads read; ${assigned.toLocaleString()} assigned to configured samples`,
        reads: streamedRecords, sampleAssignments: sampleAssignments() });
    };
    emitDemuxProgress(true, "Opening the read stream and preparing sample assignment");
    for await (const batch of streamFastq(request.file, inputHash, { signal, onProgress: (state) => {
      streamedRecords = state.records; streamedFraction = state.totalBytes ? state.compressedBytes / state.totalBytes : 0; emitDemuxProgress();
    } })) {
      const task = pool.any<ArrayBuffer>({ type: "preprocess", text: batch.text, firstOrdinal: batch.firstOrdinal })
        .then(async (buffer) => { await store.appendFrames(new Uint8Array(buffer)); emitDemuxProgress(); });
      pending.add(task); void task.then(() => pending.delete(task), () => pending.delete(task)); batches++;
      if (pending.size >= workers * 2) await Promise.race(pending);
    }
    await Promise.all(pending); if (signal.aborted) throw new DOMException("Analysis cancelled.", "AbortError");
    streamedFraction = 1; emitDemuxProgress(true, "Read streaming and sample assignment complete; calculating final filter totals");
    const statsParts = await Promise.all(pool.clients.map((_, index) => pool.at<QualityStats>(index, { type: "stats" })));
    const quality = mergeStats(statsParts, request.config.samples.length);
    const sampleCounts = quality.perSample.map(BigInt), cutoffValues = makeCutoffValues(sampleCounts, request.config.parameters.maxReadsPerSample);
    const cutoffs = makeCutoffs(sampleCounts, request.config.parameters.maxReadsPerSample);
    if (request.config.parameters.maxReadsPerSample > 0) {
      progress({ stage: "preprocessing", fraction: 1, detail: "Keeping the deterministic per-sample read selection and releasing unselected temporary records", sampleAssignments: sampleAssignments() });
      await store.compact(cutoffValues);
    }
    progress({ stage: "preprocessing", fraction: 1, detail: store.mode === "external-directory" ? "Finishing the external scratch files before UMI grouping" : "Finishing temporary read partitions before UMI grouping", sampleAssignments: sampleAssignments() });
    await store.seal();
    const spool = store.statistics();
    log.push(`${now()} preprocessing: ${quality.totalReads} raw; ${quality.qualityReads} quality; ${quality.demultiplexedReads} demultiplexed; ${batches} bounded batches`);
    log.push(`${now()} selected-read spool: ${formatBytes(spool.currentBytes)} retained; ${spool.bypassedRecords} records bypassed before disk; ${formatBytes(spool.reclaimedBytes)} reclaimed by ${spool.compactions} compaction pass${spool.compactions === 1 ? "" : "es"}`);
    if (spool.observedRecords !== quality.demultiplexedReads || spool.observedPerSample.some((count, index) => count !== sampleCounts[index]))
      throw new Error("The selected-read spool and preprocessing counters diverged; analysis stopped before UMI inference to avoid dropping a valid read.");
    stageStarted = recordTiming("preprocessing", stageStarted, quality.totalReads);

    const largestPartition = Math.max(0, ...store.sizes());
    const externalMemoryBudget = Math.max(384 * 1024 * 1024, Math.min(1536 * 1024 * 1024,
      ((navigator as Navigator & { deviceMemory?: number }).deviceMemory ?? 4) * 0.2 * 1024 * 1024 * 1024));
    const partitionWorkers = store.mode === "external-directory" && largestPartition > 0
      ? Math.max(1, Math.min(workers, Math.floor(externalMemoryBudget / (largestPartition * 3)))) : workers;
    if (partitionWorkers < workers) log.push(`${now()} partition concurrency: limited to ${partitionWorkers}/${workers} workers because the largest external-scratch partition is ${formatBytes(largestPartition)}`);

    progress({ stage: "umi", fraction: 0, detail: "Counting reads belonging to each observed UMI family" });
    const countParts: Uint8Array[] = Array(request.config.parameters.spoolPartitions);
    let countedPartitions = 0;
    await Promise.all(pool.clients.slice(0, partitionWorkers).map(async (_, worker) => {
      for (let partition = worker; partition < countParts.length; partition += partitionWorkers) {
        if (signal.aborted) throw new DOMException("Analysis cancelled.", "AbortError");
        progress({ stage: "umi", fraction: .6 * countedPartitions / Math.max(1, countParts.length), detail: `Loading temporary read block ${partition + 1} of ${countParts.length} for UMI counting` });
        const bytes = await store.readSelected(partition, cutoffValues), buffer = transferableBuffer(bytes);
        const cutoffCopy = cutoffs.slice().buffer;
        const response = await pool.at<ArrayBuffer>(worker, { type: "countFamilies", bytes: buffer, cutoffs: cutoffCopy }, [buffer, cutoffCopy]);
        countParts[partition] = new Uint8Array(response); countedPartitions++;
        progress({ stage: "umi", fraction: .6 * countedPartitions / countParts.length, detail: `Counted UMI families in ${countedPartitions} of ${countParts.length} temporary read blocks` });
      }
    }));
    progress({ stage: "umi", fraction: .64, detail: "Combining UMI counts across all temporary read blocks" });
    const mergedCounts = mergeFamilyCounts(countParts), decodedCounts = decodeFamilyCounts(mergedCounts);
    const selectedReadsBySample = Array(request.config.samples.length).fill(0) as number[];
    for (const entry of decodedCounts) selectedReadsBySample[entry.sample] += entry.count;
    const selectedReads = selectedReadsBySample.reduce((sum, count) => sum + count, 0);
    quality.downsampledReads = Math.max(0, quality.demultiplexedReads - selectedReads);
    const modelBuffer = mergedCounts.buffer.slice(mergedCounts.byteOffset, mergedCounts.byteOffset + mergedCounts.byteLength);
    progress({ stage: "umi", fraction: .72, detail: "Fitting the global UMI offspring-probability model and classifying families" });
    const familyModel = new Uint8Array(await pool.at<ArrayBuffer>(0, { type: "buildModel", bytes: modelBuffer }, [modelBuffer]));
    progress({ stage: "umi", fraction: .84, detail: "Preparing the inspectable UMI-family decisions" });
    const umiFamilies = decodeFamilyModel(familyModel, request.config);
    quality.bpbRejects = umiFamilies.filter((row) => row.disposition === "BPB-rejects").reduce((sum, row) => sum + row.familySize, 0);
    progress({ stage: "umi", fraction: .9, detail: "Sending the fitted UMI model to the consensus workers" });
    await Promise.all(pool.clients.map((_, index) => {
      const copy = familyModel.slice().buffer; return pool.at(index, { type: "initModel", bytes: copy }, [copy]);
    }));
    progress({ stage: "umi", fraction: 1, detail: `UMI grouping complete: ${umiFamilies.filter((row) => row.disposition === "likely_real").length.toLocaleString()} families selected for consensus calling` });
    log.push(`${now()} UMI model: ${umiFamilies.filter((row) => row.disposition !== "BPB-rejects").length} observed families; ${quality.bpbRejects} BPB rejects; ${umiFamilies.filter((row) => row.disposition === "likely_real").length} initially likely real`);
    stageStarted = recordTiming("umi", stageStarted, quality.demultiplexedReads - quality.downsampledReads);

    progress({ stage: "consensus", fraction: 0, detail: "Starting indel-aware consensus calling for selected UMI families" });
    const consensusParts: ReturnType<typeof decodeConsensusOutput>[] = Array(countParts.length);
    let consensusPartitions = 0, consensusSequences = 0;
    await Promise.all(pool.clients.slice(0, partitionWorkers).map(async (_, worker) => {
      for (let partition = worker; partition < consensusParts.length; partition += partitionWorkers) {
        if (signal.aborted) throw new DOMException("Analysis cancelled.", "AbortError");
        progress({ stage: "consensus", fraction: .9 * consensusPartitions / Math.max(1, consensusParts.length), detail: `Loading temporary read block ${partition + 1} of ${consensusParts.length} for consensus calling` });
        const bytes = await store.readSelected(partition, cutoffValues), buffer = transferableBuffer(bytes), cutoffCopy = cutoffs.slice().buffer;
        const response = await pool.at<ArrayBuffer>(worker, { type: "consensus", bytes: buffer, cutoffs: cutoffCopy }, [buffer, cutoffCopy]);
        consensusParts[partition] = decodeConsensusOutput(new Uint8Array(response), request.config);
        consensusPartitions++; consensusSequences += consensusParts[partition].consensuses.length;
        progress({ stage: "consensus", fraction: .9 * consensusPartitions / consensusParts.length,
          detail: `Finished ${consensusPartitions} of ${consensusParts.length} read blocks; ${consensusSequences.toLocaleString()} consensus sequences called` });
      }
    }));
    progress({ stage: "consensus", fraction: .93, detail: "Combining consensus calls and checking for heteroduplex UMI families" });
    const consensuses = consensusParts.flatMap((part) => part.consensuses).sort((a, b) => a.sampleIndex - b.sampleIndex || a.umi.localeCompare(b.umi));
    const heteroduplexes = new Set(consensusParts.flatMap((part) => part.heteroduplexes));
    const consensusByFamily = new Map<string, (typeof consensuses)[number]>(
      consensuses.map((record) => [`${record.sampleIndex}\0${record.umi}`, record]),
    );
    for (const family of umiFamilies) {
      const key = `${family.sampleIndex}\0${family.umi}`;
      if (heteroduplexes.has(key)) family.disposition = "heteroduplex";
      const consensus = consensusByFamily.get(key);
      if (consensus) family.minimumAgreement = consensus.minimumAgreement;
    }
    log.push(`${now()} consensus: ${consensuses.length} sequences; ${heteroduplexes.size} heteroduplex families`);
    progress({ stage: "consensus", fraction: .99, detail: "Consensus calls are complete; removing temporary read files" });
    await store.close();
    progress({ stage: "consensus", fraction: 1, detail: `Consensus calling complete: ${consensuses.length.toLocaleString()} sequences retained for downstream checks` });
    stageStarted = recordTiming("consensus", stageStarted, consensuses.length);

    progress({ stage: "contamination", fraction: 0, detail: "Comparing consensus sequences across samples for possible contamination" });
    const contamination = classifyContamination(consensuses, request.config);
    progress({ stage: "contamination", fraction: 1, detail: `Contamination comparison complete: ${contamination.filter((call) => call.discarded).length.toLocaleString()} sequences excluded` });
    log.push(`${now()} contamination: ${contamination.filter((call) => call.discarded).length} discarded; ${contamination.filter((call) => call.suspectOnly).length} suspect calls`);
    stageStarted = recordTiming("contamination", stageStarted, consensuses.length);
    progress({ stage: "postprocessing", fraction: 0, detail: "Starting panel screening, retained-sequence alignment, functional checks, and sequence annotation" });
    const downstreamStarted = performance.now();
    const downstream = await postprocess(consensuses, contamination, request.config, signal, runAlivibeMsa, workers,
      (state) => progress({ stage: "postprocessing", fraction: state.fraction, detail: state.detail }));
    downstream.summaries.forEach((summary, index) => {
      summary.demultiplexedReads = quality.perSample[index] ?? 0;
      summary.selectedReads = selectedReadsBySample[index] ?? 0;
      summary.downsampledReads = Math.max(0, summary.demultiplexedReads - summary.selectedReads);
      summary.observedUmis = umiFamilies.filter((family) => family.sampleIndex === index && family.disposition !== "BPB-rejects").length;
      summary.likelyRealUmis = umiFamilies.filter((family) => family.sampleIndex === index && family.disposition === "likely_real").length;
    });
    const downstreamFinished = performance.now(), downstreamSeconds = (downstreamFinished - downstreamStarted) / 1000;
    const collapseSeconds = Math.max(0, Math.min(downstreamSeconds, downstream.collapseSeconds));
    storeTiming("postprocessing", downstreamSeconds - collapseSeconds, downstream.records.length);
    const collapsedHaplotypes = Object.values(downstream.collapseGroups).reduce((sum, groups) => sum + groups.length, 0);
    progress({ stage: "collapse", fraction: 1, detail: `Collapsed identical retained sequences into ${collapsedHaplotypes.toLocaleString()} distinct haplotypes; counts represent UMI families` });
    storeTiming("collapse", collapseSeconds, collapsedHaplotypes); stageStarted = downstreamFinished;
    log.push(`${now()} collapse: ${collapsedHaplotypes} distinct haplotypes from ${downstream.records.filter((record) => record.alignedNt).length} retained UMI families; multiplicities count families, not reads`);

    const trees: Record<string, string> = {}, treeInputs = Object.entries(downstream.alignments).filter(([name]) => name.endsWith("/nucleotide"));
    let treeCursor = 0, treesFinished = 0;
    if (!request.deferPhylogeny) await Promise.all(Array.from({ length: Math.min(workers, Math.max(1, treeInputs.length)) }, async () => {
      while (true) {
        const index = treeCursor++; if (index >= treeInputs.length) return;
        const [name, alignment] = treeInputs[index];
        progress({ stage: "tree", fraction: treesFinished / Math.max(1, treeInputs.length), detail: `Inferring the collapsed phylogeny for ${name.split("/")[0]} (${treesFinished} of ${treeInputs.length} complete)` });
        try { trees[name] = await runFastTreeIsolated(alignment); }
        catch (cause) { trees[name] = starTree(alignment); log.push(`${now()} FastTree warning for ${name}: ${cause instanceof Error ? cause.message : String(cause)}; stored a zero-branch star fallback`); }
        treesFinished++;
        progress({ stage: "tree", fraction: treesFinished / Math.max(1, treeInputs.length), detail: `Finished ${treesFinished} of ${treeInputs.length} collapsed phylogenies` });
      }
    }));
    else {
      progress({ stage: "tree", fraction: 1, detail: "Phylogeny inference was deferred; retained alignments are ready for on-demand tree building" });
      log.push(`${now()} phylogeny: deferred by user; collapsed alignments are stored and trees can be inferred in the results explorer`);
    }
    if (!treeInputs.length) progress({ stage: "tree", fraction: 1, detail: "No retained alignments required phylogeny inference" });
    recordTiming("tree", stageStarted, Object.keys(trees).length);
    timings.push({ stage: "analysis-total", seconds: (performance.now() - runStarted) / 1000 });
    log.push(`${now()} postprocessing: ${downstream.records.filter((record) => record.artefactPass && record.agreementPass && record.contaminationPass && record.panelPass).length} sequences passed all non-functional filters`);
    progress({ stage: "complete", fraction: 1, detail: "Results ready" });
    return {
      schema: "webporpid-results/1",
      provenance: { webporpidVersion: "0.3.5", createdUtc: now(), engine: "C++20 WASM/WASI SIMD",
        workers, inputName: request.file.name, inputSha256: finishStreamingHash(inputHash),
        configSha256: bytesToHex(new Uint8Array(configHashBytes)), deterministicSeed: request.config.parameters.deterministicSeed.toString(),
        upstreamBranch: "nanopore", upstreamCommit: "201af7942029cfb7974880e41674be9f0ddfaf3b" },
      config: resultConfig(request.config), quality, summaries: downstream.summaries, umiFamilies,
      consensuses, contamination, contaminationReferences: request.config.contaminationPanelSequences,
      records: downstream.records, alignments: downstream.alignments, trees,
      referenceAlignments: downstream.referenceAlignments, collapseGroups: downstream.collapseGroups,
      inputMappings: request.inputMappings, runOptions: { deferPhylogeny: Boolean(request.deferPhylogeny),
        spoolStorage: request.spoolStorage === "external-directory" ? "external-directory" : "automatic" }, timings, log,
    };
  } finally {
    pool.close(); try { await store.close(); } catch { /* already closed or best-effort cleanup */ }
  }
}

self.addEventListener("message", (event: MessageEvent<RunRequest | CancelRequest>) => {
  if (event.data.type === "cancel") { cancellation?.abort(); return; }
  cancellation?.abort(); cancellation = new AbortController();
  void run(event.data, cancellation.signal).then((result) => self.postMessage({ type: "result", result }))
    .catch((cause) => self.postMessage({ type: "error", message: cause instanceof Error ? cause.message : String(cause) }));
});
