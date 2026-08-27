/// <reference lib="webworker" />

import { bytesToHex } from "@noble/hashes/utils.js";
import coreWasmUrl from "/webporpid.wasm?url";
import { classifyContamination } from "./contamination";
import { compileConfig, resultConfig } from "./config";
import { finishStreamingHash, createStreamingHash, streamFastq } from "./fastq-stream";
import { PartitionStore } from "./partition-store";
import { postprocess } from "./postprocess";
import { runFastTree } from "./biowasm";
import type { PipelineConfig, PipelineProgress, QualityStats, ResultBundle } from "./types";
import { CoreWorkerPool } from "./worker-pool";
import {
  decodeConsensusOutput, decodeFamilyCounts, decodeFamilyModel, makeCutoffs, makeCutoffValues, mergeFamilyCounts,
  mergeStats,
} from "./wasm-runtime";

type RunRequest = { type: "run"; file: File; config: PipelineConfig; workers: number };
type CancelRequest = { type: "cancel" };
let cancellation: AbortController | undefined;

function progress(value: PipelineProgress) { self.postMessage({ type: "progress", progress: value }); }
const now = () => new Date().toISOString();

async function compileCore() {
  const url = coreWasmUrl;
  try { return await WebAssembly.compileStreaming(fetch(url)); }
  catch { const response = await fetch(url); if (!response.ok) throw new Error("The webPORPID WebAssembly core could not be loaded."); return WebAssembly.compile(await response.arrayBuffer()); }
}

function starTree(fasta: string) {
  const names = [...fasta.matchAll(/^>(.*)$/gm)].map((match) => match[1].replace(/[^A-Za-z0-9_.|*+\-]/g, "_"));
  return `(${names.map((name) => `${name}:0.0`).join(",")});`;
}

async function run(request: RunRequest, signal: AbortSignal): Promise<ResultBundle> {
  const workers = Math.max(1, Math.floor(request.workers)), compiledConfig = compileConfig(request.config);
  const configForHash = compiledConfig.slice().buffer as ArrayBuffer;
  const [module, configHashBytes] = await Promise.all([compileCore(), crypto.subtle.digest("SHA-256", configForHash)]);
  const pool = await CoreWorkerPool.create(workers, module, compiledConfig), store = await PartitionStore.create(request.config.parameters.spoolPartitions);
  const inputHash = createStreamingHash(), log = [`${now()} webPORPID 0.1.0 started`,
    `${now()} execution: ${workers} WASM workers; ${store.persistent ? "OPFS" : "bounded memory"} partition spool`,
    `${now()} parameters: error_rate=${request.config.parameters.errorRate}, lengths=(${request.config.parameters.minLength},${request.config.parameters.maxLength}), lda=${request.config.parameters.ldaThreshold}`];
  try {
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
    log.push(`${now()} preprocessing: ${quality.totalReads} raw; ${quality.qualityReads} quality; ${quality.demultiplexedReads} demultiplexed; ${batches} bounded batches`);

    const sampleCounts = quality.perSample.map(BigInt), cutoffValues = makeCutoffValues(sampleCounts, request.config.parameters.maxReadsPerSample);
    const cutoffs = makeCutoffs(sampleCounts, request.config.parameters.maxReadsPerSample);
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

    progress({ stage: "contamination", fraction: 0, detail: "Building run-aware contamination database" });
    const contamination = classifyContamination(consensuses, request.config);
    log.push(`${now()} contamination: ${contamination.filter((call) => call.discarded).length} discarded; ${contamination.filter((call) => call.suspectOnly).length} suspect calls`);
    progress({ stage: "postprocessing", fraction: 0, detail: "Panel alignment, functional filter, and APOBEC model" });
    const downstream = await postprocess(consensuses, contamination, request.config, signal);
    downstream.summaries.forEach((summary, index) => {
      summary.demultiplexedReads = quality.perSample[index] ?? 0;
      summary.observedUmis = umiFamilies.filter((family) => family.sampleIndex === index && family.disposition !== "BPB-rejects").length;
      summary.likelyRealUmis = umiFamilies.filter((family) => family.sampleIndex === index && family.disposition === "likely_real").length;
    });

    const trees: Record<string, string> = {};
    for (const [name, alignment] of Object.entries(downstream.alignments).filter(([name]) => name.endsWith("/nucleotide"))) {
      progress({ stage: "tree", fraction: Object.keys(trees).length / Math.max(1, request.config.samples.length), detail: `FastTree: ${name.split("/")[0]}` });
      try { trees[name] = await runFastTree(alignment); }
      catch (cause) { trees[name] = starTree(alignment); log.push(`${now()} FastTree warning for ${name}: ${cause instanceof Error ? cause.message : String(cause)}; stored a zero-branch star fallback`); }
    }
    log.push(`${now()} postprocessing: ${downstream.records.filter((record) => record.artefactPass && record.agreementPass && record.contaminationPass && record.panelPass).length} sequences passed all non-functional filters`);
    progress({ stage: "complete", fraction: 1, detail: "Results ready" });
    return {
      schema: "webporpid-results/1",
      provenance: { webporpidVersion: "0.1.0", createdUtc: now(), engine: "C++20 WASM/WASI SIMD",
        workers, inputName: request.file.name, inputSha256: finishStreamingHash(inputHash),
        configSha256: bytesToHex(new Uint8Array(configHashBytes)), deterministicSeed: request.config.parameters.deterministicSeed.toString(),
        upstreamBranch: "nanopore", upstreamCommit: "201af7942029cfb7974880e41674be9f0ddfaf3b" },
      config: resultConfig(request.config), quality, summaries: downstream.summaries, umiFamilies,
      consensuses, contamination, records: downstream.records, alignments: downstream.alignments, trees, log,
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
