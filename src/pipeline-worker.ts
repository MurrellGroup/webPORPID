/// <reference lib="webworker" />

import { bytesToHex } from "@noble/hashes/utils.js";
import coreWasmUrl from "/webporpid.wasm?url";
import { classifyContamination } from "./contamination";
import { runAlivibeMsa } from "./alivibe-msa-runtime";
import { compileConfig, resultConfig } from "./config";
import { finishStreamingHash, createStreamingHash, streamFastq } from "./fastq-stream";
import { PartitionStore } from "./partition-store";
import { postprocess } from "./postprocess";
import { runFastTreeIsolated } from "./biowasm";
import { treeTipNames } from "./tree-names";
import type { InputFileMapping, PipelineConfig, PipelineProgress, QualityStats, ResultBundle } from "./types";
import { CoreWorkerPool } from "./worker-pool";
import {
  decodeConsensusOutput, decodeFamilyCounts, decodeFamilyModel, makeCutoffs, makeCutoffValues, mergeFamilyCounts,
  mergeStats,
} from "./wasm-runtime";

type RunRequest = { type: "run"; file: File; config: PipelineConfig; workers: number; deferPhylogeny?: boolean; inputMappings?: InputFileMapping[] };
type CancelRequest = { type: "cancel" };
let cancellation: AbortController | undefined;

function progress(value: PipelineProgress) { self.postMessage({ type: "progress", progress: value }); }
const now = () => new Date().toISOString();

function formatBytes(bytes: number) {
  const units = ["B", "KiB", "MiB", "GiB", "TiB"]; let value = Math.max(0, bytes), unit = 0;
  while (value >= 1024 && unit < units.length - 1) { value /= 1024; unit++; }
  return `${value.toFixed(unit === 0 ? 0 : value >= 10 ? 1 : 2)} ${units[unit]}`;
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
  const configForHash = compiledConfig.slice().buffer as ArrayBuffer;
  const [module, configHashBytes] = await Promise.all([compileCore(), crypto.subtle.digest("SHA-256", configForHash)]);
  const pool = await CoreWorkerPool.create(workers, module, compiledConfig), store = await PartitionStore.create(
    request.config.parameters.spoolPartitions,
    {
      sampleCount: request.config.samples.length,
      maximumReadsPerSample: request.config.parameters.maxReadsPerSample,
      requestPersistence: request.file.size >= 256 * 1024 * 1024,
    },
  );
  const inputHash = createStreamingHash(), log = [`${now()} webPORPID 0.3.1 started`,
    `${now()} execution: ${workers} WASM workers; ${store.persistent ? "OPFS" : "bounded memory"} adaptive selected-read partition spool`,
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
    const pending = new Set<Promise<void>>(); let batches = 0;
    for await (const batch of streamFastq(request.file, inputHash, { signal, onProgress: (state) => progress({
      stage: "preprocessing", fraction: state.totalBytes ? state.compressedBytes / state.totalBytes : 0,
      detail: `${state.records.toLocaleString()} FASTQ reads streamed`, reads: state.records,
    }) })) {
      const task = pool.any<ArrayBuffer>({ type: "preprocess", text: batch.text, firstOrdinal: batch.firstOrdinal })
        .then((buffer) => store.appendFrames(new Uint8Array(buffer)));
      pending.add(task); void task.then(() => pending.delete(task), () => pending.delete(task)); batches++;
      if (pending.size >= workers * 2) await Promise.race(pending);
    }
    await Promise.all(pending); if (signal.aborted) throw new DOMException("Analysis cancelled.", "AbortError");
    const statsParts = await Promise.all(pool.clients.map((_, index) => pool.at<QualityStats>(index, { type: "stats" })));
    const quality = mergeStats(statsParts, request.config.samples.length);
    const sampleCounts = quality.perSample.map(BigInt), cutoffValues = makeCutoffValues(sampleCounts, request.config.parameters.maxReadsPerSample);
    const cutoffs = makeCutoffs(sampleCounts, request.config.parameters.maxReadsPerSample);
    if (request.config.parameters.maxReadsPerSample > 0) {
      progress({ stage: "preprocessing", fraction: 1, detail: "Compacting the selected-read spool" });
      await store.compact(cutoffValues);
    }
    const spool = store.statistics();
    log.push(`${now()} preprocessing: ${quality.totalReads} raw; ${quality.qualityReads} quality; ${quality.demultiplexedReads} demultiplexed; ${batches} bounded batches`);
    log.push(`${now()} selected-read spool: ${formatBytes(spool.currentBytes)} retained; ${spool.bypassedRecords} records bypassed before disk; ${formatBytes(spool.reclaimedBytes)} reclaimed by ${spool.compactions} compaction pass${spool.compactions === 1 ? "" : "es"}`);
    if (spool.observedRecords !== quality.demultiplexedReads || spool.observedPerSample.some((count, index) => count !== sampleCounts[index]))
      throw new Error("The selected-read spool and preprocessing counters diverged; analysis stopped before UMI inference to avoid dropping a valid read.");
    stageStarted = recordTiming("preprocessing", stageStarted, quality.totalReads);

    progress({ stage: "umi", fraction: 0, detail: "Counting UMI families from disk-backed partitions" });
    const countParts: Uint8Array[] = Array(request.config.parameters.spoolPartitions);
    await Promise.all(pool.clients.map(async (_, worker) => {
      for (let partition = worker; partition < countParts.length; partition += workers) {
        if (signal.aborted) throw new DOMException("Analysis cancelled.", "AbortError");
        const bytes = await store.readSelected(partition, cutoffValues), buffer = bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength);
        const cutoffCopy = cutoffs.slice().buffer;
        const response = await pool.at<ArrayBuffer>(worker, { type: "countFamilies", bytes: buffer, cutoffs: cutoffCopy }, [buffer, cutoffCopy]);
        countParts[partition] = new Uint8Array(response);
        progress({ stage: "umi", fraction: (countParts.filter(Boolean).length) / countParts.length, detail: `Counted partition ${partition + 1}/${countParts.length}` });
      }
    }));
    const mergedCounts = mergeFamilyCounts(countParts), selectedReads = decodeFamilyCounts(mergedCounts).reduce((sum, entry) => sum + entry.count, 0);
    quality.downsampledReads = Math.max(0, quality.demultiplexedReads - selectedReads);
    const modelBuffer = mergedCounts.buffer.slice(mergedCounts.byteOffset, mergedCounts.byteOffset + mergedCounts.byteLength);
    const familyModel = new Uint8Array(await pool.at<ArrayBuffer>(0, { type: "buildModel", bytes: modelBuffer }, [modelBuffer]));
    const umiFamilies = decodeFamilyModel(familyModel, request.config);
    quality.bpbRejects = umiFamilies.filter((row) => row.disposition === "BPB-rejects").reduce((sum, row) => sum + row.familySize, 0);
    await Promise.all(pool.clients.map((_, index) => {
      const copy = familyModel.slice().buffer; return pool.at(index, { type: "initModel", bytes: copy }, [copy]);
    }));
    log.push(`${now()} UMI model: ${umiFamilies.filter((row) => row.disposition !== "BPB-rejects").length} observed families; ${quality.bpbRejects} BPB rejects; ${umiFamilies.filter((row) => row.disposition === "likely_real").length} initially likely real`);
    stageStarted = recordTiming("umi", stageStarted, quality.demultiplexedReads - quality.downsampledReads);

    progress({ stage: "consensus", fraction: 0, detail: "Generating indel-aware UMI consensuses" });
    const consensusParts: ReturnType<typeof decodeConsensusOutput>[] = Array(countParts.length);
    await Promise.all(pool.clients.map(async (_, worker) => {
      for (let partition = worker; partition < consensusParts.length; partition += workers) {
        if (signal.aborted) throw new DOMException("Analysis cancelled.", "AbortError");
        const bytes = await store.readSelected(partition, cutoffValues), buffer = bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength), cutoffCopy = cutoffs.slice().buffer;
        const response = await pool.at<ArrayBuffer>(worker, { type: "consensus", bytes: buffer, cutoffs: cutoffCopy }, [buffer, cutoffCopy]);
        consensusParts[partition] = decodeConsensusOutput(new Uint8Array(response), request.config);
        progress({ stage: "consensus", fraction: consensusParts.filter(Boolean).length / consensusParts.length, detail: `Consensus partition ${partition + 1}/${consensusParts.length}` });
      }
    }));
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
    await store.close();
    stageStarted = recordTiming("consensus", stageStarted, consensuses.length);

    progress({ stage: "contamination", fraction: 0, detail: "Building run-aware contamination database" });
    const contamination = classifyContamination(consensuses, request.config);
    log.push(`${now()} contamination: ${contamination.filter((call) => call.discarded).length} discarded; ${contamination.filter((call) => call.suspectOnly).length} suspect calls`);
    stageStarted = recordTiming("contamination", stageStarted, consensuses.length);
    progress({ stage: "postprocessing", fraction: 0, detail: "Panel alignment, functional filter, and APOBEC model" });
    const downstreamStarted = performance.now();
    const downstream = await postprocess(consensuses, contamination, request.config, signal, runAlivibeMsa, workers);
    downstream.summaries.forEach((summary, index) => {
      summary.demultiplexedReads = quality.perSample[index] ?? 0;
      summary.observedUmis = umiFamilies.filter((family) => family.sampleIndex === index && family.disposition !== "BPB-rejects").length;
      summary.likelyRealUmis = umiFamilies.filter((family) => family.sampleIndex === index && family.disposition === "likely_real").length;
    });
    const downstreamFinished = performance.now(), downstreamSeconds = (downstreamFinished - downstreamStarted) / 1000;
    const collapseSeconds = Math.max(0, Math.min(downstreamSeconds, downstream.collapseSeconds));
    storeTiming("postprocessing", downstreamSeconds - collapseSeconds, downstream.records.length);
    const collapsedHaplotypes = Object.values(downstream.collapseGroups).reduce((sum, groups) => sum + groups.length, 0);
    progress({ stage: "collapse", fraction: 1, detail: `${collapsedHaplotypes.toLocaleString()} distinct haplotypes from retained UMI families` });
    storeTiming("collapse", collapseSeconds, collapsedHaplotypes); stageStarted = downstreamFinished;
    log.push(`${now()} collapse: ${collapsedHaplotypes} distinct haplotypes from ${downstream.records.filter((record) => record.alignedNt).length} retained UMI families; multiplicities count families, not reads`);

    const trees: Record<string, string> = {}, treeInputs = Object.entries(downstream.alignments).filter(([name]) => name.endsWith("/nucleotide"));
    let treeCursor = 0, treesFinished = 0;
    if (!request.deferPhylogeny) await Promise.all(Array.from({ length: Math.min(workers, Math.max(1, treeInputs.length)) }, async () => {
      while (true) {
        const index = treeCursor++; if (index >= treeInputs.length) return;
        const [name, alignment] = treeInputs[index];
        progress({ stage: "tree", fraction: treesFinished / Math.max(1, treeInputs.length), detail: `FastTree: ${name.split("/")[0]}` });
        try { trees[name] = await runFastTreeIsolated(alignment); }
        catch (cause) { trees[name] = starTree(alignment); log.push(`${now()} FastTree warning for ${name}: ${cause instanceof Error ? cause.message : String(cause)}; stored a zero-branch star fallback`); }
        treesFinished++;
      }
    }));
    else log.push(`${now()} phylogeny: deferred by user; collapsed alignments are stored and trees can be inferred in the results explorer`);
    recordTiming("tree", stageStarted, Object.keys(trees).length);
    timings.push({ stage: "analysis-total", seconds: (performance.now() - runStarted) / 1000 });
    log.push(`${now()} postprocessing: ${downstream.records.filter((record) => record.artefactPass && record.agreementPass && record.contaminationPass && record.panelPass).length} sequences passed all non-functional filters`);
    progress({ stage: "complete", fraction: 1, detail: "Results ready" });
    return {
      schema: "webporpid-results/1",
      provenance: { webporpidVersion: "0.3.1", createdUtc: now(), engine: "C++20 WASM/WASI SIMD",
        workers, inputName: request.file.name, inputSha256: finishStreamingHash(inputHash),
        configSha256: bytesToHex(new Uint8Array(configHashBytes)), deterministicSeed: request.config.parameters.deterministicSeed.toString(),
        upstreamBranch: "nanopore", upstreamCommit: "201af7942029cfb7974880e41674be9f0ddfaf3b" },
      config: resultConfig(request.config), quality, summaries: downstream.summaries, umiFamilies,
      consensuses, contamination, records: downstream.records, alignments: downstream.alignments, trees,
      referenceAlignments: downstream.referenceAlignments, collapseGroups: downstream.collapseGroups,
      inputMappings: request.inputMappings, runOptions: { deferPhylogeny: Boolean(request.deferPhylogeny) }, timings, log,
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
