/// <reference lib="webworker" />

import { bytesToHex } from "@noble/hashes/utils.js";
import coreWasmUrl from "/webporpid.wasm?url";
import { classifyContaminationAsync } from "./contamination";
import { runAlivibeMsa } from "./alivibe-msa-runtime";
import { compileConfig, resultConfig } from "./config";
import { finishStreamingHash, createStreamingHash, streamFastq } from "./fastq-stream";
import { PartitionStore, type ExternalScratchDirectoryHandle } from "./partition-store";
import { collapsePostprocess, postprocess, type PostprocessOutput } from "./postprocess";
import { runFastTreeIsolated } from "./biowasm";
import { downstreamResources, statusRecord } from "./optional-stages";
import { applyThresholdSelection, buildConsensusThresholdReview, buildUmiThresholdReview } from "./threshold-review";
import { treeTipNames } from "./tree-names";
import type { InputFileMapping, OptionalStageName, PipelineConfig, PipelineProgress, QualityStats, ResultBundle, SampleSummary,
  ThresholdReview, ThresholdSelection } from "./types";
import { CoreWorkerPool } from "./worker-pool";
import {
  decodeConsensusOutput, decodeFamilyCounts, decodeFamilyModel, encodeFamilyModel, makeCutoffs, makeCutoffValues, mergeFamilyCounts,
  mergeStats,
} from "./wasm-runtime";

type RunRequest = { type: "run"; file: File; config: PipelineConfig; workers: number; deferPhylogeny?: boolean; deferContamination?: boolean;
  deferPostprocessing?: boolean; deferCollapse?: boolean; interactiveFiltering?: boolean; inputMappings?: InputFileMapping[];
  spoolStorage?: "automatic" | "external-directory"; scratchDirectory?: ExternalScratchDirectoryHandle };
type CancelRequest = { type: "cancel" };
type SkipStageRequest = { type: "skip-stage"; stage: OptionalStageName };
type ThresholdSelectionRequest = { type: "threshold-selection"; selection: ThresholdSelection };
let cancellation: AbortController | undefined;
let activeOptionalStage: OptionalStageName | undefined, optionalStageCancellation: AbortController | undefined;
let requestedSkip: OptionalStageName | undefined;
let thresholdWaiter: { id: string; accept(selection: ThresholdSelection): void; cancel(): void } | undefined;

function progress(value: PipelineProgress) { self.postMessage({ type: "progress", progress: value }); }
const now = () => new Date().toISOString();

function requestThresholdSelection(review: ThresholdReview, signal: AbortSignal): Promise<ThresholdSelection> {
  if (thresholdWaiter) throw new Error("An interactive threshold checkpoint is already open.");
  self.postMessage({ type: "threshold-review", review });
  return new Promise((resolve, reject) => {
    const cleanup = () => { signal.removeEventListener("abort", abort); thresholdWaiter = undefined; };
    const abort = () => { cleanup(); reject(new DOMException("Analysis cancelled.", "AbortError")); };
    thresholdWaiter = { id: review.id,
      accept(selection) { cleanup(); resolve(selection); },
      cancel: abort };
    signal.addEventListener("abort", abort, { once: true });
  });
}

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

async function optionalStage<T>(stage: OptionalStageName, runSignal: AbortSignal,
  operation: (signal: AbortSignal) => Promise<T>): Promise<{ skipped: false; value: T } | { skipped: true }> {
  activeOptionalStage = stage; requestedSkip = undefined; optionalStageCancellation = new AbortController();
  const controller = optionalStageCancellation;
  const cancelWithRun = () => controller.abort();
  runSignal.addEventListener("abort", cancelWithRun, { once: true });
  try {
    return { skipped: false, value: await operation(controller.signal) };
  } catch (cause) {
    if (runSignal.aborted) throw new DOMException("Analysis cancelled.", "AbortError");
    if (!runSignal.aborted && controller.signal.aborted && requestedSkip === stage) return { skipped: true };
    throw cause;
  } finally {
    runSignal.removeEventListener("abort", cancelWithRun);
    if (activeOptionalStage === stage) activeOptionalStage = undefined;
    if (optionalStageCancellation === controller) optionalStageCancellation = undefined;
  }
}

async function run(request: RunRequest, signal: AbortSignal): Promise<ResultBundle> {
  const runStarted = performance.now(), timings: NonNullable<ResultBundle["timings"]> = [];
  const thresholdSelections: NonNullable<ResultBundle["thresholdSelections"]> = [];
  let thresholdReviewPauseMs = 0;
  const workers = Math.max(1, Math.floor(request.workers)), compiledConfig = compileConfig(request.config);
  if (request.spoolStorage === "external-directory" && !request.scratchDirectory)
    throw new Error("External scratch storage was selected, but no writable directory handle was provided.");
  // Include browser/CLI analysis settings that are intentionally not part of
  // the fixed C++ configuration ABI (for example the panel-filter strategy).
  const configForHash = new TextEncoder().encode(JSON.stringify(resultConfig(request.config))).buffer as ArrayBuffer;
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
	  const inputHash = createStreamingHash(), log = [`${now()} webPORPID 0.3.10 started`,
    `${now()} execution: ${workers} WASM workers; ${storageLabel} ${request.config.parameters.maxReadsPerSample > 0 ? "adaptive selected-read" : "all-read"} partition spool`,
    `${now()} parameters: error_rate=${request.config.parameters.errorRate}, lengths=(${request.config.parameters.minLength},${request.config.parameters.maxLength}), lda=${request.config.parameters.ldaThreshold}, panel_filter=${request.config.parameters.panelFilterMode ?? "mafft-batch"}`];
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

    const countParts: Uint8Array[] = Array(request.config.parameters.spoolPartitions);
    const readBlocks: NonNullable<PipelineProgress["readBlocks"]> = Array(countParts.length).fill("waiting");
    progress({ stage: "umi", fraction: 0, detail: "Counting reads belonging to each observed UMI family", readBlocks: [...readBlocks] });
    let countedPartitions = 0;
    await Promise.all(pool.clients.slice(0, partitionWorkers).map(async (_, worker) => {
      for (let partition = worker; partition < countParts.length; partition += partitionWorkers) {
        if (signal.aborted) throw new DOMException("Analysis cancelled.", "AbortError");
        progress({ stage: "umi", fraction: .6 * countedPartitions / Math.max(1, countParts.length), detail: `Loading temporary read block ${partition + 1} of ${countParts.length} for UMI counting`, readBlocks: [...readBlocks] });
        const bytes = await store.readSelected(partition, cutoffValues), buffer = transferableBuffer(bytes);
        readBlocks[partition] = "loaded";
        progress({ stage: "umi", fraction: .6 * countedPartitions / Math.max(1, countParts.length),
          detail: `Read block ${partition + 1} of ${countParts.length} loaded; counting its UMI families`, readBlocks: [...readBlocks] });
        const cutoffCopy = cutoffs.slice().buffer;
        const response = await pool.at<ArrayBuffer>(worker, { type: "countFamilies", bytes: buffer, cutoffs: cutoffCopy }, [buffer, cutoffCopy]);
        countParts[partition] = new Uint8Array(response); countedPartitions++;
        progress({ stage: "umi", fraction: .6 * countedPartitions / countParts.length, detail: `Counted UMI families in ${countedPartitions} of ${countParts.length} temporary read blocks`, readBlocks: [...readBlocks] });
      }
    }));
    progress({ stage: "umi", fraction: .64, detail: "Combining UMI counts across all temporary read blocks", readBlocks: [...readBlocks] });
    const mergedCounts = mergeFamilyCounts(countParts), decodedCounts = decodeFamilyCounts(mergedCounts);
    const selectedReadsBySample = Array(request.config.samples.length).fill(0) as number[];
    for (const entry of decodedCounts) selectedReadsBySample[entry.sample] += entry.count;
    const selectedReads = selectedReadsBySample.reduce((sum, count) => sum + count, 0);
    quality.downsampledReads = Math.max(0, quality.demultiplexedReads - selectedReads);
    const modelBuffer = mergedCounts.buffer.slice(mergedCounts.byteOffset, mergedCounts.byteOffset + mergedCounts.byteLength);
    progress({ stage: "umi", fraction: .72, detail: "Fitting the global UMI offspring-probability model and classifying families", readBlocks: [...readBlocks] });
    let familyModel = new Uint8Array(await pool.at<ArrayBuffer>(0, { type: "buildModel", bytes: modelBuffer }, [modelBuffer]));
    progress({ stage: "umi", fraction: .84, detail: "Preparing the inspectable UMI-family decisions", readBlocks: [...readBlocks] });
    const umiFamilies = decodeFamilyModel(familyModel, request.config);
    quality.bpbRejects = umiFamilies.filter((row) => row.disposition === "BPB-rejects").reduce((sum, row) => sum + row.familySize, 0);
    if (request.interactiveFiltering) {
      progress({ stage: "umi", fraction: .86, detail: "UMI probability calculations are complete; waiting for threshold review", readBlocks: [...readBlocks] });
      const review = buildUmiThresholdReview(umiFamilies, request.config);
      const pauseStarted = performance.now();
      const selection = await requestThresholdSelection(review, signal);
      const pauseMs = performance.now() - pauseStarted;
      thresholdReviewPauseMs += pauseMs;
      stageStarted += pauseMs; // Human review time is not computational UMI time.
      if (selection.id !== review.id || selection.phase !== review.phase) throw new Error("The interactive UMI-threshold response did not match the open checkpoint.");
      const accepted = applyThresholdSelection(request.config, umiFamilies, selection); thresholdSelections.push(accepted);
      familyModel = new Uint8Array(encodeFamilyModel(umiFamilies));
      for (const change of accepted.changes) log.push(`${now()} interactive UMI threshold: ${change}`);
      progress({ stage: "umi", fraction: .88, detail: "Accepted UMI thresholds; updating family decisions for every observed UMI", readBlocks: [...readBlocks] });
    }
    progress({ stage: "umi", fraction: .9, detail: "Sending the fitted UMI model to the consensus workers", readBlocks: [...readBlocks] });
    await Promise.all(pool.clients.map((_, index) => {
      const copy = familyModel.slice().buffer; return pool.at(index, { type: "initModel", bytes: copy }, [copy]);
    }));
    progress({ stage: "umi", fraction: 1, detail: `UMI grouping complete: ${umiFamilies.filter((row) => row.disposition === "likely_real").length.toLocaleString()} families selected for consensus calling`, readBlocks: [...readBlocks] });
    log.push(`${now()} UMI model: ${umiFamilies.filter((row) => row.disposition !== "BPB-rejects").length} observed families; ${quality.bpbRejects} BPB rejects; ${umiFamilies.filter((row) => row.disposition === "likely_real").length} initially likely real`);
    stageStarted = recordTiming("umi", stageStarted, quality.demultiplexedReads - quality.downsampledReads);

    progress({ stage: "consensus", fraction: 0, detail: "Starting indel-aware consensus calling for selected UMI families", readBlocks: [...readBlocks] });
    const consensusParts: ReturnType<typeof decodeConsensusOutput>[] = Array(countParts.length);
    let consensusPartitions = 0, consensusSequences = 0;
    await Promise.all(pool.clients.slice(0, partitionWorkers).map(async (_, worker) => {
      for (let partition = worker; partition < consensusParts.length; partition += partitionWorkers) {
        if (signal.aborted) throw new DOMException("Analysis cancelled.", "AbortError");
        progress({ stage: "consensus", fraction: .9 * consensusPartitions / Math.max(1, consensusParts.length), detail: `Loading temporary read block ${partition + 1} of ${consensusParts.length} for consensus calling`, readBlocks: [...readBlocks] });
        const bytes = await store.readSelected(partition, cutoffValues), buffer = transferableBuffer(bytes), cutoffCopy = cutoffs.slice().buffer;
        const response = await pool.at<ArrayBuffer>(worker, { type: "consensus", bytes: buffer, cutoffs: cutoffCopy }, [buffer, cutoffCopy]);
        consensusParts[partition] = decodeConsensusOutput(new Uint8Array(response), request.config);
        readBlocks[partition] = "complete"; consensusPartitions++; consensusSequences += consensusParts[partition].consensuses.length;
        progress({ stage: "consensus", fraction: .9 * consensusPartitions / consensusParts.length,
          detail: `Finished ${consensusPartitions} of ${consensusParts.length} read blocks; ${consensusSequences.toLocaleString()} consensus sequences called`, readBlocks: [...readBlocks] });
      }
    }));
    progress({ stage: "consensus", fraction: .93, detail: "Combining consensus calls and checking for heteroduplex UMI families", readBlocks: [...readBlocks] });
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
    progress({ stage: "consensus", fraction: .99, detail: "Consensus calls are complete; removing temporary read files", readBlocks: [...readBlocks] });
    await store.close();
    progress({ stage: "consensus", fraction: 1, detail: `Consensus calling complete: ${consensuses.length.toLocaleString()} sequences retained for downstream checks`, readBlocks: [...readBlocks] });
    stageStarted = recordTiming("consensus", stageStarted, consensuses.length);

    // One pass avoids repeatedly scanning every family and consensus for every sample on large multiplexed runs.
    const observedBySample = new Uint32Array(request.config.samples.length), likelyBySample = new Uint32Array(request.config.samples.length);
    const consensusBySample = new Uint32Array(request.config.samples.length);
    for (const family of umiFamilies) {
      if (family.disposition !== "BPB-rejects") observedBySample[family.sampleIndex]++;
      if (family.disposition === "likely_real") likelyBySample[family.sampleIndex]++;
    }
    for (const record of consensuses) consensusBySample[record.sampleIndex]++;
    const baseSummaries: SampleSummary[] = request.config.samples.map((sample, index) => ({
      sample: sample.name, demultiplexedReads: quality.perSample[index] ?? 0,
      selectedReads: selectedReadsBySample[index] ?? 0,
      downsampledReads: Math.max(0, (quality.perSample[index] ?? 0) - (selectedReadsBySample[index] ?? 0)),
      observedUmis: observedBySample[index], likelyRealUmis: likelyBySample[index], consensusSequences: consensusBySample[index],
    }));
    const optionalStages: NonNullable<ResultBundle["optionalStages"]> = {
      contamination: statusRecord("deferred", "Waiting after consensus."),
      postprocessing: statusRecord("deferred", "Waiting after consensus."),
      collapse: statusRecord("deferred", "Waiting for downstream filtering."),
      tree: statusRecord("deferred", "Waiting for haplotype collapse."),
    };
    let contamination: ResultBundle["contamination"] = [];
    let postprocessingContaminationMode: ResultBundle["postprocessingContaminationMode"];
    let downstream: PostprocessOutput = { records: [], summaries: baseSummaries, alignments: {}, referenceAlignments: {}, collapseGroups: {}, collapseSeconds: 0 };
    const trees: Record<string, string> = {};

    if (request.deferContamination && request.config.parameters.contaminationFilter) {
      const detail = "Deferred by user after consensus; no contamination decision has been assigned. Other requested stages continue without applying this filter.";
      optionalStages.contamination = statusRecord("deferred", detail);
      progress({ stage: "contamination", fraction: 1, detail: "Contamination checks deferred; continuing without excluding any sequence as contamination" });
      log.push(`${now()} contamination: deferred by user; contamination gate bypassed for downstream work`);
    } else {
      progress({ stage: "contamination", fraction: 0, detail: "Preparing run-wide sequence signatures for contamination checks" });
      const started = performance.now(), result = await optionalStage("contamination", signal, (stageSignal) =>
        classifyContaminationAsync(consensuses, request.config, stageSignal, (state) => progress({ stage: "contamination",
          fraction: state.fraction, detail: state.detail })));
      storeTiming("contamination", (performance.now() - started) / 1000, consensuses.length);
      if (result.skipped) {
        optionalStages.contamination = statusRecord("skipped", "Skipped by user while contamination checks were running; no decisions were retained and no sequence is excluded at this gate.");
        progress({ stage: "contamination", fraction: 1, detail: "Contamination checks skipped; continuing without excluding any sequence as contamination" });
        log.push(`${now()} contamination: skipped by user; contamination gate bypassed for downstream work`);
      } else {
        contamination = result.value;
        const discardedBySample = new Uint32Array(request.config.samples.length); let discardedCount = 0, suspectCount = 0;
        const sampleIndexByName = new Map(request.config.samples.map((sample, index) => [sample.name, index]));
        for (const call of contamination) {
          if (call.discarded) { discardedCount++; const sampleIndex = sampleIndexByName.get(call.sample); if (sampleIndex !== undefined) discardedBySample[sampleIndex]++; }
          if (call.suspectOnly) suspectCount++;
        }
        optionalStages.contamination = statusRecord("completed", request.config.parameters.contaminationFilter
          ? `${discardedCount} consensus sequences excluded.`
          : "Contamination filtering was disabled in the run configuration.");
        baseSummaries.forEach((summary, index) => { summary.contaminationPassed = summary.consensusSequences - discardedBySample[index]; });
        progress({ stage: "contamination", fraction: 1, detail: `Contamination checks complete: ${discardedCount.toLocaleString()} sequences excluded` });
        log.push(`${now()} contamination: ${discardedCount} discarded; ${suspectCount} suspect calls`);
      }
    }

    const contaminationApplied = optionalStages.contamination.state === "completed" && request.config.parameters.contaminationFilter;
    if (request.interactiveFiltering) {
      progress({ stage: "postprocessing", fraction: 0,
        detail: "Consensus and contamination eligibility are complete; waiting for the downstream threshold review" });
      const discardedIds = contaminationApplied
        ? new Set(contamination.filter((call) => call.discarded).map((call) => call.sequenceId))
        : new Set<string>();
      const review = buildConsensusThresholdReview(consensuses, discardedIds, request.config, umiFamilies);
      const pauseStarted = performance.now();
      const selection = await requestThresholdSelection(review, signal);
      thresholdReviewPauseMs += performance.now() - pauseStarted;
      if (selection.id !== review.id || selection.phase !== review.phase) throw new Error("The interactive consensus-filter response did not match the open checkpoint.");
      const accepted = applyThresholdSelection(request.config, umiFamilies, selection); thresholdSelections.push(accepted);
      for (const change of accepted.changes) log.push(`${now()} interactive consensus filter: ${change}`);
      progress({ stage: "postprocessing", fraction: 0,
        detail: "Accepted consensus-family thresholds; continuing to reference-panel screening and retained-family alignment" });
    }
    if (request.deferPostprocessing) {
        optionalStages.postprocessing = statusRecord("deferred", "Deferred by user; panel screening, retained-family alignment, and annotations have not run.");
        optionalStages.collapse = statusRecord("deferred", "Waiting for downstream filtering.");
        optionalStages.tree = statusRecord("deferred", "Waiting for haplotype collapse.");
        progress({ stage: "postprocessing", fraction: 1, detail: "Alignment and downstream filtering deferred; it can be computed from the stored consensus calls" });
        log.push(`${now()} postprocessing: deferred by user`);
    } else {
        const contaminationNote = contaminationApplied ? "Computed contamination decisions will be applied." : "The contamination gate is bypassed; every consensus remains eligible.";
        progress({ stage: "postprocessing", fraction: 0, detail: `Starting panel screening, retained-family alignment, and sequence annotation. ${contaminationNote}` });
        const started = performance.now(), result = await optionalStage("postprocessing", signal, (stageSignal) =>
          postprocess(consensuses, contaminationApplied ? contamination : [], request.config, stageSignal, runAlivibeMsa, workers,
            (state) => progress({ stage: "postprocessing", fraction: state.fraction, detail: state.detail }), { collapse: false }));
        storeTiming("postprocessing", (performance.now() - started) / 1000, result.skipped ? undefined : result.value.records.length);
        if (result.skipped) {
          optionalStages.postprocessing = statusRecord("skipped", "Skipped by user while downstream filtering was running; partial decisions were discarded.");
          optionalStages.collapse = statusRecord("deferred", "Waiting for downstream filtering.");
          optionalStages.tree = statusRecord("deferred", "Waiting for haplotype collapse.");
          progress({ stage: "postprocessing", fraction: 1, detail: "Alignment and downstream filtering skipped; partial outputs were not retained" });
          log.push(`${now()} postprocessing: skipped by user; partial outputs discarded`);
        } else {
          downstream = result.value;
          postprocessingContaminationMode = contaminationApplied ? "applied" : "bypassed";
          downstream.summaries.forEach((summary, index) => {
            const base = baseSummaries[index]; summary.demultiplexedReads = base.demultiplexedReads;
            summary.selectedReads = base.selectedReads; summary.downsampledReads = base.downsampledReads;
            summary.observedUmis = base.observedUmis; summary.likelyRealUmis = base.likelyRealUmis;
            if (!contaminationApplied) delete summary.contaminationPassed;
          });
          optionalStages.postprocessing = statusRecord("completed", `${downstream.records.length} consensus-family records evaluated${contaminationApplied ? " with contamination decisions applied" : "; contamination was bypassed and excluded zero sequences"}.`);
          progress({ stage: "postprocessing", fraction: 1, detail: `Panel screening, retained-family alignment, and annotations complete${contaminationApplied ? "" : "; contamination was not applied"}` });
          log.push(`${now()} postprocessing: ${downstream.records.filter((record) => record.artefactPass && record.agreementPass && record.contaminationPass && record.panelPass).length} sequences passed all non-functional filters; contamination=${contaminationApplied ? "applied" : "bypassed"}`);
        }
    }

    if (optionalStages.postprocessing.state === "completed") {
      if (request.deferCollapse) {
        optionalStages.collapse = statusRecord("deferred", "Deferred by user; uncollapsed retained-family alignments are stored.");
        optionalStages.tree = statusRecord("deferred", "Waiting for haplotype collapse.");
        progress({ stage: "collapse", fraction: 1, detail: "Haplotype collapse deferred; uncollapsed alignments are ready for later computation" });
        log.push(`${now()} collapse: deferred by user`);
      } else {
        progress({ stage: "collapse", fraction: 0, detail: "Starting identical-haplotype collapse and collapsed-variant functional filtering; multiplicities count UMI families, not reads" });
        const started = performance.now(), result = await optionalStage("collapse", signal, (stageSignal) =>
          collapsePostprocess(downstream, request.config, stageSignal,
            (state) => progress({ stage: "collapse", fraction: state.fraction, detail: state.detail })));
        storeTiming("collapse", (performance.now() - started) / 1000,
          result.skipped ? undefined : Object.values(result.value.collapseGroups).reduce((sum, groups) => sum + groups.length, 0));
        if (result.skipped) {
          optionalStages.collapse = statusRecord("skipped", "Skipped by user; uncollapsed retained-family alignments remain available.");
          optionalStages.tree = statusRecord("deferred", "Waiting for haplotype collapse.");
          progress({ stage: "collapse", fraction: 1, detail: "Haplotype collapse skipped; uncollapsed alignments remain available" });
          log.push(`${now()} collapse: skipped by user`);
        } else {
          downstream = result.value;
          const collapsedHaplotypes = Object.values(downstream.collapseGroups).reduce((sum, groups) => sum + groups.length, 0);
          const functionalHaplotypes = downstream.summaries.reduce((sum, summary) => sum + (summary.functionalPassed ?? 0), 0);
          optionalStages.collapse = statusRecord("completed", `${collapsedHaplotypes} haplotypes; ${functionalHaplotypes} collapsed variants passed configured functional filters; multiplicities count UMI families.`);
          progress({ stage: "collapse", fraction: 1, detail: `Collapsed retained sequences into ${collapsedHaplotypes.toLocaleString()} distinct haplotypes; ${functionalHaplotypes.toLocaleString()} passed configured functional filters; counts represent UMI families` });
          log.push(`${now()} collapse: ${collapsedHaplotypes} distinct haplotypes from ${downstream.records.filter((record) => record.alignedNt).length} retained UMI families; ${functionalHaplotypes} collapsed functional passes; multiplicities count families, not reads`);
        }
      }
    }

    if (optionalStages.collapse.state === "completed") {
      const treeInputs = request.config.samples.flatMap((sample) => {
        const key = `${sample.name}/nucleotide`, alignment = downstream.alignments[key]; return alignment ? [[key, alignment] as const] : [];
      });
      if (request.deferPhylogeny) {
        optionalStages.tree = statusRecord("deferred", "Deferred by user; collapsed alignments are stored for on-demand inference.");
        progress({ stage: "tree", fraction: 1, detail: "Phylogeny inference deferred; collapsed alignments are ready for on-demand tree building" });
        log.push(`${now()} phylogeny: deferred by user; collapsed alignments are stored and trees can be inferred in the results explorer`);
      } else if (!treeInputs.length) {
        optionalStages.tree = statusRecord("completed", "No retained collapsed alignments required a tree.");
        progress({ stage: "tree", fraction: 1, detail: "No retained alignments required phylogeny inference" });
      } else {
        let treeCursor = 0, treesFinished = 0;
        const started = performance.now(), result = await optionalStage("tree", signal, async (stageSignal) => {
          await Promise.all(Array.from({ length: Math.min(workers, treeInputs.length) }, async () => {
            while (true) {
              if (stageSignal.aborted) throw new DOMException("Phylogeny inference skipped.", "AbortError");
              const index = treeCursor++; if (index >= treeInputs.length) return;
              const [name, alignment] = treeInputs[index];
              progress({ stage: "tree", fraction: treesFinished / treeInputs.length, detail: `Inferring the collapsed phylogeny for ${name.split("/")[0]} (${treesFinished} of ${treeInputs.length} complete)` });
              try { trees[name] = await runFastTreeIsolated(alignment, stageSignal); }
              catch (cause) {
                if (stageSignal.aborted) throw cause;
                trees[name] = starTree(alignment); log.push(`${now()} FastTree warning for ${name}: ${cause instanceof Error ? cause.message : String(cause)}; stored a zero-branch star fallback`);
              }
              treesFinished++;
              progress({ stage: "tree", fraction: treesFinished / treeInputs.length, detail: `Finished ${treesFinished} of ${treeInputs.length} collapsed phylogenies` });
            }
          }));
        });
        storeTiming("tree", (performance.now() - started) / 1000, Object.keys(trees).length);
        if (result.skipped) {
          optionalStages.tree = statusRecord("skipped", `Skipped by user after ${Object.keys(trees).length} of ${treeInputs.length} trees completed.`);
          progress({ stage: "tree", fraction: 1, detail: `Phylogeny inference skipped; ${Object.keys(trees).length} completed trees were retained` });
          log.push(`${now()} phylogeny: skipped by user after ${Object.keys(trees).length}/${treeInputs.length} trees`);
        } else optionalStages.tree = statusRecord("completed", `${Object.keys(trees).length} collapsed phylogenies inferred.`);
      }
    }

    timings.push({ stage: "analysis-total", seconds: Math.max(0, performance.now() - runStarted - thresholdReviewPauseMs) / 1000 });
    progress({ stage: "complete", fraction: 1, detail: "Results ready" });
    return {
      schema: "webporpid-results/1",
      provenance: { webporpidVersion: "0.3.10", createdUtc: now(), engine: "C++20 WASM/WASI SIMD",
        workers, inputName: request.file.name, inputSha256: finishStreamingHash(inputHash),
        configSha256: bytesToHex(new Uint8Array(configHashBytes)), deterministicSeed: request.config.parameters.deterministicSeed.toString(),
        upstreamBranch: "nanopore", upstreamCommit: "201af7942029cfb7974880e41674be9f0ddfaf3b" },
      config: resultConfig(request.config), quality, summaries: downstream.summaries, umiFamilies,
      consensuses, contamination, contaminationReferences: request.config.contaminationPanelSequences,
      downstreamResources: downstreamResources(request.config),
      records: downstream.records, alignments: downstream.alignments, trees,
      referenceAlignments: downstream.referenceAlignments, collapseGroups: downstream.collapseGroups,
      inputMappings: request.inputMappings, runOptions: { deferPhylogeny: Boolean(request.deferPhylogeny),
        deferContamination: Boolean(request.deferContamination), deferPostprocessing: Boolean(request.deferPostprocessing),
        deferCollapse: Boolean(request.deferCollapse),
        interactiveFiltering: Boolean(request.interactiveFiltering),
        spoolStorage: request.spoolStorage === "external-directory" ? "external-directory" : "automatic" },
      optionalStages, postprocessingContaminationMode, timings, thresholdSelections, log,
    };
  } finally {
    pool.close(); try { await store.close(); } catch { /* already closed or best-effort cleanup */ }
  }
}

self.addEventListener("message", (event: MessageEvent<RunRequest | CancelRequest | SkipStageRequest | ThresholdSelectionRequest>) => {
  if (event.data.type === "cancel") { cancellation?.abort(); optionalStageCancellation?.abort(); thresholdWaiter?.cancel(); return; }
  if (event.data.type === "skip-stage") {
    if (event.data.stage === activeOptionalStage) { requestedSkip = event.data.stage; optionalStageCancellation?.abort(); }
    return;
  }
  if (event.data.type === "threshold-selection") {
    if (event.data.selection.id === thresholdWaiter?.id) thresholdWaiter.accept(event.data.selection);
    return;
  }
  cancellation?.abort(); cancellation = new AbortController();
  void run(event.data, cancellation.signal).then((result) => self.postMessage({ type: "result", result }))
    .catch((cause) => self.postMessage({ type: "error", message: cause instanceof Error ? cause.message : String(cause) }));
});
