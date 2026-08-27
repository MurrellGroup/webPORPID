import { createHash } from "node:crypto";
import { createReadStream } from "node:fs";
import { mkdir, mkdtemp, open, readFile, rm, stat, writeFile } from "node:fs/promises";
import { availableParallelism, tmpdir } from "node:os";
import { basename, dirname, isAbsolute, join, resolve } from "node:path";
import { createInterface } from "node:readline";
import { createGunzip } from "node:zlib";
import { Worker as NodeWorker } from "node:worker_threads";

import { bytesToHex } from "@noble/hashes/utils.js";
import { classifyContamination } from "../src/contamination.ts";
import { compileConfig, parseConfigYaml, resolveReferenceFiles, resultConfig } from "../src/config.ts";
import { postprocess } from "../src/postprocess.ts";
import { decodeResult, encodeResult, exportComponent, safeDatasetName } from "../src/result-file.ts";
import { concatenateSpoolRecords, parseSpoolRecordHeader, selectedSpoolRecord, SPOOL_HEADER_BYTES } from "../src/spool-record.ts";
import {
  decodeConsensusOutput, decodeFamilyCounts, decodeFamilyModel, makeCutoffs, makeCutoffValues, mergeFamilyCounts, mergeStats,
} from "../src/wasm-runtime.ts";
import { createFastTreeRunner } from "./direct-fasttree.mjs";
import { createMsaRunner } from "./direct-msa.mjs";

const VERSION = "0.1.1";
const UPSTREAM_COMMIT = "201af7942029cfb7974880e41674be9f0ddfaf3b";
const CLI_DIRECTORY = dirname(new URL(import.meta.url).pathname);

export function defaultCliAssets() {
  const directory = join(CLI_DIRECTORY, "assets");
  return { wasmPath: join(directory, "webporpid.wasm"), msaPath: join(directory, "alivibe-msa.wasm"),
    fastTreeJavascriptPath: join(directory, "fasttree.cjs"), fastTreeWasmPath: join(directory, "fasttree.wasm") };
}

function usage() {
  return `porpid-cli ${VERSION}\n\n` +
    `Run the complete nanopore/PacBio pipeline:\n` +
    `  porpid-cli run reads.fastq.gz --config config.yaml --output results.webporpid [--workers N]\n\n` +
    `Inspect or export a saved analysis:\n` +
    `  porpid-cli inspect results.webporpid\n` +
    `  porpid-cli export results.webporpid --component consensus-fasta [--sample NAME] --output consensus.fasta\n\n` +
    `Workers default to all logical CPUs (${availableParallelism()}). Temporary read partitions are streamed to disk and removed after consensus.\n` +
    `Components: consensus-fasta, passed-consensus-fasta, rejected-consensus-fasta, trimmed-nt-fasta, trimmed-aa-fasta,\n` +
    `            family-csv, low-agreement-csv, contamination-csv, postproc-csv, apobec-csv,\n` +
    `            nucleotide-alignment, protein-alignment, newick, log`;
}

function option(args, name) {
  const index = args.indexOf(name); if (index >= 0) return args[index + 1];
  const inline = args.find((value) => value.startsWith(`${name}=`)); return inline?.slice(name.length + 1);
}

function integer(value, label, fallback) {
  if (value == null) return fallback;
  if (!/^\d+$/.test(value) || Number(value) < 1 || !Number.isSafeInteger(Number(value))) throw new Error(`${label} requires an integer of at least one.`);
  return Number(value);
}

function status(message) { process.stderr.write(`[webPORPID] ${message}\n`); }
const now = () => new Date().toISOString();

class WorkerClient {
  constructor(worker, webWorker) {
    this.worker = worker; this.webWorker = webWorker; this.pending = new Map(); this.nextId = 1; this.tail = Promise.resolve();
    const receive = (message) => { const pending = this.pending.get(message.id); if (!pending) return; this.pending.delete(message.id); message.error ? pending.reject(new Error(message.error)) : pending.resolve(message.result); };
    const fail = (cause) => { const error = cause instanceof Error ? cause : new Error(cause?.message ?? String(cause)); for (const pending of this.pending.values()) pending.reject(error); this.pending.clear(); };
    if (webWorker) { worker.addEventListener("message", (event) => receive(event.data)); worker.addEventListener("error", fail); }
    else { worker.on("message", receive); worker.on("error", fail); }
  }
  raw(message, transfer = []) {
    return new Promise((resolvePromise, reject) => {
      const id = this.nextId++; this.pending.set(id, { resolve: resolvePromise, reject });
      if (this.webWorker) this.worker.postMessage({ id, ...message }, { transfer }); else this.worker.postMessage({ id, ...message }, transfer);
    });
  }
  call(message, transfer = []) {
    const start = () => this.raw(message, transfer), result = this.tail.then(start, start); this.tail = result.catch(() => {}); return result;
  }
  terminate() { return Promise.resolve(this.worker.terminate()); }
}

class WorkerPool {
  constructor(clients) { this.clients = clients; this.cursor = 0; }
  static async create(size, wasmPath, compiledConfig) {
    const clients = [];
    for (let index = 0; index < size; index++) {
      const webWorker = Boolean(process.versions.bun && globalThis.Worker);
      const worker = webWorker ? new globalThis.Worker(new URL("./porpid-worker.mjs", import.meta.url))
        : new NodeWorker(new URL("./porpid-worker.mjs", import.meta.url));
      clients.push(new WorkerClient(worker, webWorker));
    }
    await Promise.all(clients.map((client) => { const copy = compiledConfig.slice().buffer; return client.call({ type: "init", wasmPath, config: copy }, [copy]); }));
    return new WorkerPool(clients);
  }
  any(message, transfer = []) { return this.clients[this.cursor++ % this.clients.length].call(message, transfer); }
  at(index, message, transfer = []) { return this.clients[index].call(message, transfer); }
  async close() { await Promise.all(this.clients.map((client) => client.terminate())); }
}

class DiskPartitions {
  static async create(count) {
    let base = process.env.WEBPORPID_TMPDIR || tmpdir();
    try { await mkdir(base, { recursive: true }); } catch { base = process.cwd(); }
    const directory = await mkdtemp(join(base, "webporpid-")), handles = [], paths = [];
    for (let index = 0; index < count; index++) { const path = join(directory, `partition-${index}.bin`); paths.push(path); handles.push(await open(path, "w+")); }
    return new DiskPartitions(directory, handles, paths, Array.from({ length: count }, () => Promise.resolve()), Array(count).fill(0));
  }
  constructor(directory, handles, paths, chains, lengths) { Object.assign(this, { directory, handles, paths, chains, lengths }); this.closed = false; }
  async appendFrames(bytes) {
    const grouped = new Map(); const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength); let offset = 0;
    while (offset < bytes.byteLength) {
      if (offset + 5 > bytes.byteLength) throw new Error("Truncated spool routing frame.");
      const partition = bytes[offset++], length = view.getUint32(offset, true); offset += 4;
      if (partition >= this.handles.length || offset + length > bytes.byteLength) throw new Error("Invalid spool routing frame.");
      const parts = grouped.get(partition) ?? []; parts.push(bytes.slice(offset, offset + length)); grouped.set(partition, parts); offset += length;
    }
    await Promise.all([...grouped].map(([partition, parts]) => {
      const length = parts.reduce((sum, part) => sum + part.length, 0), joined = new Uint8Array(length); let at = 0;
      for (const part of parts) { joined.set(part, at); at += part.length; }
      const operation = this.chains[partition].then(async () => {
        const result = await this.handles[partition].write(joined, 0, joined.length, this.lengths[partition]);
        if (result.bytesWritten !== joined.length) throw new Error("Temporary partition write was incomplete."); this.lengths[partition] += result.bytesWritten;
      });
      this.chains[partition] = operation; return operation;
    }));
  }
  async readFully(partition, output, position) {
    let offset = 0;
    while (offset < output.byteLength) {
      const result = await this.handles[partition].read(output, offset, output.byteLength - offset, position + offset);
      if (!result.bytesRead) throw new Error("Temporary partition read was incomplete."); offset += result.bytesRead;
    }
  }
  async readSelected(partition, cutoffs) {
    await this.chains[partition]; const records = [], headerBytes = new Uint8Array(SPOOL_HEADER_BYTES); let offset = 0;
    while (offset < this.lengths[partition]) {
      await this.readFully(partition, headerBytes, offset); const header = parseSpoolRecordHeader(headerBytes);
      if (offset + header.recordLength > this.lengths[partition]) throw new Error("Temporary partition contains a truncated spool record.");
      if (selectedSpoolRecord(header, cutoffs)) {
        const record = new Uint8Array(header.recordLength); record.set(headerBytes);
        await this.readFully(partition, record.subarray(SPOOL_HEADER_BYTES), offset + SPOOL_HEADER_BYTES); records.push(record);
      }
      offset += header.recordLength;
    }
    if (offset !== this.lengths[partition]) throw new Error("Temporary partition contains trailing spool bytes.");
    return concatenateSpoolRecords(records);
  }
  async close() {
    if (this.closed) return; this.closed = true;
    await Promise.all(this.chains); await Promise.all(this.handles.map((handle) => handle.close()));
    await rm(this.directory, { recursive: true, force: true });
  }
}

async function* fastqBatches(path, hash, batchRecords = 256, maximumBatchBytes = 4 * 1024 * 1024) {
  const raw = createReadStream(path, { highWaterMark: 1024 * 1024 }); raw.on("data", (chunk) => hash.update(chunk));
  const input = path.toLowerCase().endsWith(".gz") ? raw.pipe(createGunzip({ chunkSize: 1024 * 1024 })) : raw;
  const lines = createInterface({ input, crlfDelay: Infinity }); let record = [], batch = "", count = 0, records = 0, firstOrdinal = 1;
  for await (const line of lines) {
    record.push(line); if (record.length !== 4) continue;
    if (!record[0].startsWith("@") || !record[2].startsWith("+") || record[1].length !== record[3].length)
      throw new Error(`Malformed FASTQ record ${records + 1}.`);
    batch += `${record.join("\n")}\n`; record = []; count++; records++;
    if (count >= batchRecords || batch.length >= maximumBatchBytes) {
      yield { text: batch, count, firstOrdinal, records }; firstOrdinal += count; batch = ""; count = 0;
    }
  }
  if (record.length) throw new Error("The FASTQ stream is truncated at its final record.");
  if (count) yield { text: batch, count, firstOrdinal, records };
  if (!records) throw new Error("The input contains no FASTQ records.");
}

async function loadConfiguration(path) {
  const source = await readFile(path, "utf8"), config = parseConfigYaml(source), base = dirname(resolve(path));
  const requested = new Set([config.contaminationPanel]);
  for (const sample of config.samples) { requested.add(sample.panel); if (sample.functionalReference) requested.add(sample.functionalReference); }
  const files = new Map([...requested].map((name) => [name, () => readFile(isAbsolute(name) ? name : resolve(base, name), "utf8")]));
  return resolveReferenceFiles(config, files);
}

async function runPipeline({ inputPath, configPath, outputPath, workers, assets }) {
  const input = resolve(inputPath), configuration = resolve(configPath); await stat(input);
  const config = await loadConfiguration(configuration), compiledConfig = compileConfig(config);
  const configHash = createHash("sha256").update(compiledConfig).digest("hex"), inputHash = createHash("sha256");
  status(`starting ${config.dataset} with ${workers} workers`);
  const pool = await WorkerPool.create(workers, assets.wasmPath, compiledConfig), store = await DiskPartitions.create(config.parameters.spoolPartitions);
  const log = [`${now()} webPORPID ${VERSION} started`, `${now()} execution: ${workers} WASM workers; disk-backed partition spool`,
    `${now()} parameters: error_rate=${config.parameters.errorRate}, lengths=(${config.parameters.minLength},${config.parameters.maxLength}), lda=${config.parameters.ldaThreshold}`];
  try {
    const pending = new Set(); let batches = 0, streamed = 0;
    for await (const batch of fastqBatches(input, inputHash)) {
      const task = pool.any({ type: "preprocess", text: batch.text, firstOrdinal: batch.firstOrdinal })
        .then((buffer) => store.appendFrames(new Uint8Array(buffer)));
      pending.add(task); void task.then(() => pending.delete(task), () => pending.delete(task)); batches++; streamed = batch.records;
      if (batches % 20 === 0) status(`preprocessing: ${streamed.toLocaleString()} reads streamed`);
      if (pending.size >= workers * 2) await Promise.race(pending);
    }
    await Promise.all(pending);
    const quality = mergeStats(await Promise.all(pool.clients.map((_, index) => pool.at(index, { type: "stats" }))), config.samples.length);
    log.push(`${now()} preprocessing: ${quality.totalReads} raw; ${quality.qualityReads} quality; ${quality.demultiplexedReads} demultiplexed; ${batches} bounded batches`);
    status(`preprocessing complete: ${quality.demultiplexedReads.toLocaleString()} demultiplexed reads`);

    const sampleCounts = quality.perSample.map(BigInt), cutoffValues = makeCutoffValues(sampleCounts, config.parameters.maxReadsPerSample);
    const cutoffs = makeCutoffs(sampleCounts, config.parameters.maxReadsPerSample), countParts = Array(config.parameters.spoolPartitions);
    await Promise.all(pool.clients.map(async (_, worker) => {
      for (let partition = worker; partition < countParts.length; partition += workers) {
        const bytes = await store.readSelected(partition, cutoffValues), data = bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength), cutoffCopy = cutoffs.slice().buffer;
        countParts[partition] = new Uint8Array(await pool.at(worker, { type: "countFamilies", bytes: data, cutoffs: cutoffCopy }, [data, cutoffCopy]));
      }
    }));
    const mergedCounts = mergeFamilyCounts(countParts), selectedReads = decodeFamilyCounts(mergedCounts).reduce((sum, entry) => sum + entry.count, 0);
    quality.downsampledReads = Math.max(0, quality.demultiplexedReads - selectedReads);
    const modelData = mergedCounts.buffer.slice(mergedCounts.byteOffset, mergedCounts.byteOffset + mergedCounts.byteLength);
    const familyModel = new Uint8Array(await pool.at(0, { type: "buildModel", bytes: modelData }, [modelData]));
    const umiFamilies = decodeFamilyModel(familyModel, config);
    quality.bpbRejects = umiFamilies.filter((row) => row.disposition === "BPB-rejects").reduce((sum, row) => sum + row.familySize, 0);
    await Promise.all(pool.clients.map((_, index) => { const copy = familyModel.slice().buffer; return pool.at(index, { type: "initModel", bytes: copy }, [copy]); }));
    log.push(`${now()} UMI model: ${umiFamilies.filter((row) => row.disposition !== "BPB-rejects").length} observed families; ${quality.bpbRejects} BPB rejects; ${umiFamilies.filter((row) => row.disposition === "likely_real").length} initially likely real`);
    status(`UMI model complete: ${umiFamilies.filter((row) => row.disposition !== "BPB-rejects").length.toLocaleString()} observed families`);

    const consensusParts = Array(countParts.length);
    await Promise.all(pool.clients.map(async (_, worker) => {
      for (let partition = worker; partition < consensusParts.length; partition += workers) {
        const bytes = await store.readSelected(partition, cutoffValues), data = bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength), cutoffCopy = cutoffs.slice().buffer;
        if (process.env.WEBPORPID_DEBUG) status(`debug: worker ${worker} consensus partition ${partition} bytes=${bytes.byteLength}`);
        const response = await pool.at(worker, { type: "consensus", bytes: data, cutoffs: cutoffCopy }, [data, cutoffCopy]);
        consensusParts[partition] = decodeConsensusOutput(new Uint8Array(response), config);
        if (process.env.WEBPORPID_DEBUG) status(`debug: worker ${worker} finished partition ${partition}`);
      }
    }));
    const consensuses = consensusParts.flatMap((part) => part.consensuses).sort((a, b) => a.sampleIndex - b.sampleIndex || a.umi.localeCompare(b.umi));
    const heteroduplexes = new Set(consensusParts.flatMap((part) => part.heteroduplexes));
    const consensusByFamily = new Map(consensuses.map((record) => [`${record.sampleIndex}\0${record.umi}`, record]));
    for (const family of umiFamilies) {
      const key = `${family.sampleIndex}\0${family.umi}`; if (heteroduplexes.has(key)) family.disposition = "heteroduplex";
      const consensus = consensusByFamily.get(key); if (consensus) family.minimumAgreement = consensus.minimumAgreement;
    }
    log.push(`${now()} consensus: ${consensuses.length} sequences; ${heteroduplexes.size} heteroduplex families`);
    status(`consensus complete: ${consensuses.length.toLocaleString()} sequences`); await store.close();

    const contamination = classifyContamination(consensuses, config);
    log.push(`${now()} contamination: ${contamination.filter((call) => call.discarded).length} discarded; ${contamination.filter((call) => call.suspectOnly).length} suspect calls`);
    const downstream = await postprocess(consensuses, contamination, config, undefined, createMsaRunner(assets.msaPath));
    downstream.summaries.forEach((summary, index) => {
      summary.demultiplexedReads = quality.perSample[index] ?? 0;
      summary.observedUmis = umiFamilies.filter((family) => family.sampleIndex === index && family.disposition !== "BPB-rejects").length;
      summary.likelyRealUmis = umiFamilies.filter((family) => family.sampleIndex === index && family.disposition === "likely_real").length;
    });
    const trees = {}, fastTree = createFastTreeRunner(assets.fastTreeJavascriptPath, assets.fastTreeWasmPath);
    for (const [name, alignment] of Object.entries(downstream.alignments).filter(([name]) => name.endsWith("/nucleotide"))) {
      status(`FastTree: ${name.split("/")[0]}`); trees[name] = await fastTree(alignment);
    }
    log.push(`${now()} postprocessing: ${downstream.records.filter((record) => record.artefactPass && record.agreementPass && record.contaminationPass && record.panelPass).length} sequences passed all non-functional filters`);
    const result = { schema: "webporpid-results/1", provenance: { webporpidVersion: VERSION, createdUtc: now(), engine: "C++20 WASM/WASI SIMD",
      workers, inputName: basename(input), inputSha256: inputHash.digest("hex"), configSha256: configHash,
      deterministicSeed: config.parameters.deterministicSeed.toString(), upstreamBranch: "nanopore", upstreamCommit: UPSTREAM_COMMIT },
      config: resultConfig(config), quality, summaries: downstream.summaries, umiFamilies, consensuses, contamination,
      records: downstream.records, alignments: downstream.alignments, trees, log };
    await mkdir(dirname(resolve(outputPath)), { recursive: true }); await writeFile(outputPath, encodeResult(result));
    status(`wrote ${outputPath}`); return result;
  } finally { await pool.close(); await store.close(); }
}

async function inspect(path) {
  const result = decodeResult(new Uint8Array(await readFile(path)));
  process.stdout.write(JSON.stringify({ schema: result.schema, provenance: result.provenance, quality: result.quality, summaries: result.summaries,
    components: { consensuses: result.consensuses.length, families: result.umiFamilies.length, contaminationCalls: result.contamination.length,
      records: result.records.length, alignments: Object.keys(result.alignments), trees: Object.keys(result.trees) } }, null, 2) + "\n");
}

async function exportResult(args) {
  const path = args[0], component = option(args, "--component"), sample = option(args, "--sample");
  if (!path || !component) throw new Error("export requires a results file and --component.");
  const result = decodeResult(new Uint8Array(await readFile(path))), exported = exportComponent(result, component, sample);
  const output = option(args, "--output") ?? `${safeDatasetName(result.config.dataset)}${sample ? `-${safeDatasetName(sample)}` : ""}.${exported.extension}`;
  await writeFile(output, exported.text); status(`wrote ${output}`);
}

export async function runCli(overrideAssets) {
  const args = process.argv.slice(2), command = args.shift();
  if (!command || command === "--help" || command === "-h" || command === "help") { process.stdout.write(`${usage()}\n`); return; }
  if (command === "--version" || command === "version") { process.stdout.write(`${VERSION}\n`); return; }
  if (command === "inspect") { if (!args[0]) throw new Error("inspect requires a .webporpid file."); await inspect(args[0]); return; }
  if (command === "export") { await exportResult(args); return; }
  if (command !== "run") throw new Error(`Unknown command ${command}.\n\n${usage()}`);
  const inputPath = args[0], configPath = option(args, "--config");
  if (!inputPath || inputPath.startsWith("--") || !configPath) throw new Error("run requires an input FASTQ and --config config.yaml.");
  const outputPath = option(args, "--output") ?? option(args, "--out") ?? `${basename(inputPath).replace(/\.(fastq|fq)(\.gz)?$/i, "")}.webporpid`;
  await runPipeline({ inputPath, configPath, outputPath, workers: integer(option(args, "--workers"), "--workers", availableParallelism()), assets: overrideAssets ?? defaultCliAssets() });
}
