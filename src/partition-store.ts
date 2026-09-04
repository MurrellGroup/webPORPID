import {
  concatenateSpoolRecords, parseSpoolRecordHeader, selectedSpoolRecord, selectSpoolChunks, SPOOL_HEADER_BYTES,
  type SpoolRecordHeader,
} from "./spool-record.ts";

const MAX_SAMPLING_HASH = (1n << 64n) - 1n;
const DEFAULT_COMPACTION_INTERVAL = 512 * 1024 * 1024;

function decodeBrowserFrames(bytes: Uint8Array): Array<{ partition: number; record: Uint8Array }> {
  const output: Array<{ partition: number; record: Uint8Array }> = [];
  const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength); let offset = 0;
  while (offset < bytes.byteLength) {
    if (offset + 5 > bytes.byteLength) throw new Error("Truncated spool routing frame.");
    const partition = bytes[offset++], length = view.getUint32(offset, true); offset += 4;
    if (offset + length > bytes.byteLength) throw new Error("Invalid spool routing frame length.");
    output.push({ partition, record: bytes.slice(offset, offset + length) }); offset += length;
  }
  return output;
}

interface PartitionBackend {
  readonly compactable: boolean;
  append(partition: number, chunks: Uint8Array[]): Promise<void>;
  compact(cutoffs: readonly bigint[]): Promise<number>;
  seal(): Promise<void>;
  readSelected(partition: number, cutoffs: readonly bigint[]): Promise<Uint8Array>;
  replaceWithResult(partition: number, bytes: Uint8Array): Promise<void>;
  readResult(partition: number): Promise<Uint8Array>;
  sizes(): number[];
  close(): Promise<void>;
}

export interface ExternalScratchWritable {
  write(data: Uint8Array): Promise<void>;
  close(): Promise<void>;
  abort(reason?: unknown): Promise<void>;
}

export interface ExternalScratchFileHandle {
  createWritable(options?: { keepExistingData?: boolean }): Promise<ExternalScratchWritable>;
  getFile(): Promise<Blob>;
}

export interface ExternalScratchDirectoryHandle {
  readonly name: string;
  getDirectoryHandle(name: string, options?: { create?: boolean }): Promise<ExternalScratchDirectoryHandle>;
  getFileHandle(name: string, options?: { create?: boolean }): Promise<ExternalScratchFileHandle>;
  removeEntry(name: string, options?: { recursive?: boolean }): Promise<void>;
  queryPermission?(descriptor: { mode: "readwrite" }): Promise<PermissionState>;
  requestPermission?(descriptor: { mode: "readwrite" }): Promise<PermissionState>;
}

export type PartitionStorageMode = "opfs" | "memory" | "external-directory";

export interface PartitionStoreOptions {
  sampleCount?: number;
  maximumReadsPerSample?: number;
  compactionIntervalBytes?: number;
  requestPersistence?: boolean;
  externalDirectory?: ExternalScratchDirectoryHandle;
}

export interface BrowserStorageInformation {
  quotaBytes?: number;
  usageBytes?: number;
  persisted?: boolean;
}

export interface PartitionStoreStatistics {
  observedRecords: number;
  observedPerSample: bigint[];
  admittedRecords: number;
  bypassedRecords: number;
  currentBytes: number;
  writtenBytes: number;
  reclaimedBytes: number;
  compactions: number;
}

export interface SyncByteWriter {
  write(buffer: Uint8Array, options: { at: number }): number;
}

/**
 * OPFS sync access handles may legally return a short write. Keep advancing until
 * every byte has landed; a zero-progress writer is treated as an I/O failure
 * instead of silently leaving a truncated record in the partition.
 */
export function writeAllSync(handle: SyncByteWriter, bytes: Uint8Array, at: number) {
  let offset = 0, noProgress = 0;
  while (offset < bytes.byteLength) {
    const written = handle.write(bytes.subarray(offset), { at: at + offset });
    if (!Number.isInteger(written) || written < 0 || written > bytes.byteLength - offset)
      throw new Error("OPFS returned an invalid spool write length.");
    if (written === 0) {
      if (++noProgress >= 3) throw new Error("OPFS could not make progress while writing the selected-read spool.");
      continue;
    }
    noProgress = 0; offset += written;
  }
}

export function adaptiveSpoolCutoff(observed: bigint, maximumReadsPerSample: number) {
  if (maximumReadsPerSample < 1 || observed <= BigInt(maximumReadsPerSample)) return MAX_SAMPLING_HASH;
  return (MAX_SAMPLING_HASH * BigInt(maximumReadsPerSample)) / observed;
}

export function admitSpoolRecord(header: SpoolRecordHeader, observed: bigint, maximumReadsPerSample: number) {
  return header.samplingHash <= adaptiveSpoolCutoff(observed, maximumReadsPerSample);
}

function formatBytes(bytes: number) {
  if (!Number.isFinite(bytes) || bytes < 0) return "unknown";
  const units = ["B", "KiB", "MiB", "GiB", "TiB"]; let value = bytes, unit = 0;
  while (value >= 1024 && unit < units.length - 1) { value /= 1024; unit++; }
  return `${value.toFixed(unit === 0 ? 0 : value >= 10 ? 1 : 2)} ${units[unit]}`;
}

async function browserStorageInformation(storage?: StorageManager): Promise<BrowserStorageInformation> {
  if (!storage) return {};
  const information: BrowserStorageInformation = {};
  try {
    const estimate = await storage.estimate(); information.quotaBytes = estimate.quota; information.usageBytes = estimate.usage;
  } catch { /* storage estimates are advisory */ }
  try { information.persisted = await storage.persisted(); } catch { /* not implemented by every browser */ }
  return information;
}

async function opfsWriteFailure(cause: unknown, currentBytes: number) {
  const storage = typeof navigator === "undefined" ? undefined : navigator.storage;
  const information = await browserStorageInformation(storage), available = information.quotaBytes == null || information.usageBytes == null
    ? undefined : Math.max(0, information.quotaBytes - information.usageBytes);
  const reason = cause instanceof Error ? cause.message : String(cause);
  return new Error(
    `Browser temporary storage could not continue writing the selected-read spool after ${formatBytes(currentBytes)} for this run`
    + `${available == null ? "" : ` (${formatBytes(available)} currently available to this site)`}. `
    + "Free storage, lower ‘Maximum reads / sample’, or run porpid-cli with WEBPORPID_TMPDIR on a disk with more space. "
    + `Storage error: ${reason}`,
    { cause },
  );
}

class MemoryBackend implements PartitionBackend {
  readonly compactable = true;
  private parts: Uint8Array[][];
  private lengths: number[];
  private maximumBytes: number;
  private resultReady: boolean[];
  constructor(count: number, maximumBytes = 512 * 1024 * 1024) {
    this.maximumBytes = maximumBytes; this.parts = Array.from({ length: count }, () => []); this.lengths = Array(count).fill(0);
    this.resultReady = Array(count).fill(false);
  }
  async append(partition: number, chunks: Uint8Array[]) {
    const added = chunks.reduce((sum, chunk) => sum + chunk.byteLength, 0);
    if (this.lengths.reduce((sum, value) => sum + value, 0) + added > this.maximumBytes)
      throw new Error("Browser-local temporary storage is unavailable and the selected-read in-memory spool exceeded 512 MiB. Use a browser with OPFS support, lower ‘Maximum reads / sample’, or use porpid-cli.");
    this.parts[partition].push(...chunks); this.lengths[partition] += added;
  }
  async compact(cutoffs: readonly bigint[]) {
    let reclaimed = 0;
    for (let partition = 0; partition < this.parts.length; partition++) {
      const retained = this.parts[partition].filter((record) => selectedSpoolRecord(parseSpoolRecordHeader(record), cutoffs));
      const length = retained.reduce((sum, record) => sum + record.byteLength, 0);
      reclaimed += this.lengths[partition] - length; this.parts[partition] = retained; this.lengths[partition] = length;
    }
    return reclaimed;
  }
  async seal() { /* memory is immediately readable */ }
  async readSelected(partition: number, cutoffs: readonly bigint[]) { return selectSpoolChunks(this.parts[partition], cutoffs); }
  async replaceWithResult(partition: number, bytes: Uint8Array) {
    const total = this.lengths.reduce((sum, value) => sum + value, 0) - this.lengths[partition] + bytes.byteLength;
    if (total > this.maximumBytes) throw new Error("Consensus results exceeded the bounded in-memory scratch limit.");
    this.parts[partition] = [bytes]; this.lengths[partition] = bytes.byteLength; this.resultReady[partition] = true;
  }
  async readResult(partition: number) {
    if (!this.resultReady[partition]) throw new Error("A consensus-result partition is not ready.");
    const bytes = concatenateChunks(this.parts[partition]);
    this.parts[partition] = []; this.lengths[partition] = 0; this.resultReady[partition] = false;
    return bytes;
  }
  sizes() { return [...this.lengths]; }
  async close() { this.parts = []; this.lengths = []; this.resultReady = []; }
}

class OpfsBackend implements PartitionBackend {
  readonly compactable = true;
  private closed = false;
  private root: FileSystemDirectoryHandle;
  private directoryName: string;
  private handles: FileSystemSyncAccessHandle[];
  private lengths: number[];
  private resultReady: boolean[];
  private constructor(
    root: FileSystemDirectoryHandle,
    directoryName: string,
    handles: FileSystemSyncAccessHandle[],
    lengths: number[],
  ) { this.root = root; this.directoryName = directoryName; this.handles = handles; this.lengths = lengths; this.resultReady = Array(lengths.length).fill(false); }

  static async create(count: number) {
    const root = await navigator.storage.getDirectory();
    const directoryName = `webporpid-${crypto.randomUUID()}`;
    const directory = await root.getDirectoryHandle(directoryName, { create: true });
    const handles: FileSystemSyncAccessHandle[] = [];
    try {
      for (let index = 0; index < count; index++) {
        const file = await directory.getFileHandle(`partition-${index}.bin`, { create: true });
        handles.push(await file.createSyncAccessHandle());
      }
      return new OpfsBackend(root, directoryName, handles, Array(count).fill(0));
    } catch (cause) {
      for (const handle of handles) try { handle.close(); } catch { /* best-effort rollback */ }
      try { await root.removeEntry(directoryName, { recursive: true }); } catch { /* best-effort rollback */ }
      throw cause;
    }
  }

  async append(partition: number, chunks: Uint8Array[]) {
    const handle = this.handles[partition]; let offset = this.lengths[partition];
    try {
      for (const chunk of chunks) { writeAllSync(handle, chunk, offset); offset += chunk.byteLength; }
      this.lengths[partition] = offset;
    } catch (cause) { throw await opfsWriteFailure(cause, this.lengths.reduce((sum, value) => sum + value, 0)); }
  }
  private readFully(handle: FileSystemSyncAccessHandle, output: Uint8Array, at: number) {
    let offset = 0;
    while (offset < output.byteLength) {
      const read = handle.read(output.subarray(offset), { at: at + offset });
      if (!read) throw new Error("OPFS returned an incomplete spool record."); offset += read;
    }
  }
  async compact(cutoffs: readonly bigint[]) {
    let reclaimed = 0;
    for (let partition = 0; partition < this.handles.length; partition++) {
      const handle = this.handles[partition], originalLength = this.lengths[partition];
      const headerBytes = new Uint8Array(SPOOL_HEADER_BYTES); let readOffset = 0, writeOffset = 0;
      try {
        while (readOffset < originalLength) {
          this.readFully(handle, headerBytes, readOffset); const header = parseSpoolRecordHeader(headerBytes);
          if (readOffset + header.recordLength > originalLength) throw new Error("OPFS contains a truncated spool record.");
          if (selectedSpoolRecord(header, cutoffs)) {
            if (writeOffset !== readOffset) {
              const record = new Uint8Array(header.recordLength); record.set(headerBytes);
              this.readFully(handle, record.subarray(SPOOL_HEADER_BYTES), readOffset + SPOOL_HEADER_BYTES);
              writeAllSync(handle, record, writeOffset);
            }
            writeOffset += header.recordLength;
          }
          readOffset += header.recordLength;
        }
        if (readOffset !== originalLength) throw new Error("OPFS contains trailing spool bytes.");
        if (writeOffset !== originalLength) handle.truncate(writeOffset);
        this.lengths[partition] = writeOffset; reclaimed += originalLength - writeOffset;
      } catch (cause) { throw await opfsWriteFailure(cause, this.lengths.reduce((sum, value) => sum + value, 0)); }
    }
    return reclaimed;
  }
  async seal() { for (const handle of this.handles) handle.flush(); }
  async readSelected(partition: number, cutoffs: readonly bigint[]) {
    const handle = this.handles[partition], records: Uint8Array[] = [], headerBytes = new Uint8Array(SPOOL_HEADER_BYTES);
    let offset = 0;
    while (offset < this.lengths[partition]) {
      this.readFully(handle, headerBytes, offset); const header = parseSpoolRecordHeader(headerBytes);
      if (offset + header.recordLength > this.lengths[partition]) throw new Error("OPFS contains a truncated spool record.");
      if (selectedSpoolRecord(header, cutoffs)) {
        const record = new Uint8Array(header.recordLength); record.set(headerBytes);
        this.readFully(handle, record.subarray(SPOOL_HEADER_BYTES), offset + SPOOL_HEADER_BYTES); records.push(record);
      }
      offset += header.recordLength;
    }
    if (offset !== this.lengths[partition]) throw new Error("OPFS contains trailing spool bytes.");
    return concatenateSpoolRecords(records);
  }
  async replaceWithResult(partition: number, bytes: Uint8Array) {
    const handle = this.handles[partition];
    try {
      handle.truncate(0); writeAllSync(handle, bytes, 0); handle.flush();
      this.lengths[partition] = bytes.byteLength; this.resultReady[partition] = true;
    } catch (cause) { throw await opfsWriteFailure(cause, this.lengths.reduce((sum, value) => sum + value, 0)); }
  }
  async readResult(partition: number) {
    if (!this.resultReady[partition]) throw new Error("A consensus-result partition is not ready.");
    const handle = this.handles[partition], bytes = new Uint8Array(this.lengths[partition]);
    this.readFully(handle, bytes, 0); handle.truncate(0);
    this.lengths[partition] = 0; this.resultReady[partition] = false; return bytes;
  }
  sizes() { return [...this.lengths]; }
  async close() {
    if (this.closed) return; this.closed = true;
    let firstFailure: unknown;
    for (const handle of this.handles) {
      try { handle.flush(); } catch (cause) { firstFailure ??= cause; }
      try { handle.close(); } catch (cause) { firstFailure ??= cause; }
    }
    this.handles = [];
    try { await this.root.removeEntry(this.directoryName, { recursive: true }); } catch (cause) { firstFailure ??= cause; }
    if (firstFailure) throw firstFailure;
  }
}

function concatenateChunks(chunks: readonly Uint8Array[]) {
  if (chunks.length === 1) return chunks[0];
  const length = chunks.reduce((sum, chunk) => sum + chunk.byteLength, 0), output = new Uint8Array(length);
  let offset = 0; for (const chunk of chunks) { output.set(chunk, offset); offset += chunk.byteLength; }
  return output;
}

function selectSpoolBuffer(bytes: Uint8Array, cutoffs: readonly bigint[]) {
  const selected: Uint8Array[] = []; let offset = 0;
  while (offset < bytes.byteLength) {
    const header = parseSpoolRecordHeader(bytes.subarray(offset));
    if (offset + header.recordLength > bytes.byteLength) throw new Error("An external scratch partition contains a truncated spool record.");
    if (selectedSpoolRecord(header, cutoffs)) selected.push(bytes.subarray(offset, offset + header.recordLength));
    offset += header.recordLength;
  }
  if (offset !== bytes.byteLength) throw new Error("An external scratch partition contains trailing spool bytes.");
  return selected.length === 1 && selected[0].byteLength === bytes.byteLength ? bytes : concatenateSpoolRecords(selected);
}

class ExternalDirectoryBackend implements PartitionBackend {
  readonly compactable = false;
  private readonly root: ExternalScratchDirectoryHandle;
  private readonly directoryName: string;
  private readonly directory: ExternalScratchDirectoryHandle;
  private readonly fileHandles: Array<ExternalScratchFileHandle | undefined>;
  private readonly writers: Array<ExternalScratchWritable | undefined>;
  private readonly lengths: number[];
  private readonly pending: Uint8Array[][];
  private readonly pendingBytes: number[];
  private readonly resultReady: boolean[];
  private readonly flushThreshold: number;
  private sealed = false;
  private closed = false;

  private constructor(
    root: ExternalScratchDirectoryHandle,
    directoryName: string,
    directory: ExternalScratchDirectoryHandle,
    count: number,
  ) {
    this.root = root; this.directoryName = directoryName; this.directory = directory;
    this.fileHandles = Array(count); this.writers = Array(count); this.lengths = Array(count).fill(0);
    this.pending = Array.from({ length: count }, () => []); this.pendingBytes = Array(count).fill(0);
    this.resultReady = Array(count).fill(false);
    this.flushThreshold = Math.max(256 * 1024, Math.min(4 * 1024 * 1024, Math.floor((64 * 1024 * 1024) / count)));
  }

  static async create(count: number, root: ExternalScratchDirectoryHandle) {
    const directoryName = `webporpid-scratch-${Date.now()}-${crypto.randomUUID()}`;
    try {
      const directory = await root.getDirectoryHandle(directoryName, { create: true });
      return new ExternalDirectoryBackend(root, directoryName, directory, count);
    } catch (cause) {
      throw new Error(`The selected scratch directory could not be prepared. ${cause instanceof Error ? cause.message : String(cause)}`, { cause });
    }
  }

  private async writer(partition: number) {
    if (this.sealed) throw new Error("The external scratch spool is already sealed for reading.");
    let writer = this.writers[partition];
    if (!writer) {
      const file = await this.directory.getFileHandle(`partition-${partition}.bin`, { create: true });
      writer = await file.createWritable({ keepExistingData: false });
      this.fileHandles[partition] = file; this.writers[partition] = writer;
    }
    return writer;
  }

  private async flush(partition: number) {
    const chunks = this.pending[partition]; if (!chunks.length) return;
    const bytes = concatenateChunks(chunks); this.pending[partition] = []; this.pendingBytes[partition] = 0;
    await (await this.writer(partition)).write(bytes);
  }

  async append(partition: number, chunks: Uint8Array[]) {
    if (!chunks.length) return;
    try {
      const added = chunks.reduce((sum, chunk) => sum + chunk.byteLength, 0);
      this.pending[partition].push(...chunks); this.pendingBytes[partition] += added;
      this.lengths[partition] += added;
      if (this.pendingBytes[partition] >= this.flushThreshold) await this.flush(partition);
    } catch (cause) {
      throw new Error(
        `The selected scratch disk could not continue writing after ${formatBytes(this.lengths.reduce((sum, value) => sum + value, 0))}. `
        + `Check its free space and write permission. ${cause instanceof Error ? cause.message : String(cause)}`,
        { cause },
      );
    }
  }

  async compact(_cutoffs: readonly bigint[]) { return 0; }

  async seal() {
    if (this.sealed) return;
    const flushes = await Promise.allSettled(this.pending.map((_, partition) => this.flush(partition)));
    const flushFailure = flushes.find((result): result is PromiseRejectedResult => result.status === "rejected");
    if (flushFailure) throw new Error(`The external scratch spool could not flush its final blocks. ${flushFailure.reason instanceof Error ? flushFailure.reason.message : String(flushFailure.reason)}`, { cause: flushFailure.reason });
    const active = this.writers.map((writer, index) => ({ writer, index })).filter((entry): entry is { writer: ExternalScratchWritable; index: number } => Boolean(entry.writer));
    const results = await Promise.allSettled(active.map((entry) => entry.writer.close()));
    results.forEach((result, index) => { if (result.status === "fulfilled") this.writers[active[index].index] = undefined; });
    const failed = results.find((result): result is PromiseRejectedResult => result.status === "rejected");
    if (failed) throw new Error(`The external scratch spool could not be finalized for reading. ${failed.reason instanceof Error ? failed.reason.message : String(failed.reason)}`, { cause: failed.reason });
    this.sealed = true;
  }

  async readSelected(partition: number, cutoffs: readonly bigint[]) {
    if (!this.sealed) throw new Error("The external scratch spool must be sealed before it can be read.");
    const handle = this.fileHandles[partition]; if (!handle) return new Uint8Array();
    const file = await handle.getFile();
    if (file.size !== this.lengths[partition]) throw new Error("An external scratch partition has an unexpected size.");
    const bytes = new Uint8Array(await file.arrayBuffer());
    return cutoffs.every((cutoff) => cutoff === MAX_SAMPLING_HASH) ? bytes : selectSpoolBuffer(bytes, cutoffs);
  }

  async replaceWithResult(partition: number, bytes: Uint8Array) {
    if (!this.sealed) throw new Error("The external scratch spool must be sealed before storing consensus results.");
    try {
      const handle = this.fileHandles[partition]
        ?? await this.directory.getFileHandle(`partition-${partition}.bin`, { create: true });
      const writer = await handle.createWritable({ keepExistingData: false });
      try { await writer.write(bytes); await writer.close(); }
      catch (cause) { await writer.abort(cause).catch(() => undefined); throw cause; }
      this.fileHandles[partition] = handle; this.lengths[partition] = bytes.byteLength; this.resultReady[partition] = true;
    } catch (cause) {
      throw new Error(`The selected scratch disk could not store a consensus-result block. ${cause instanceof Error ? cause.message : String(cause)}`, { cause });
    }
  }

  async readResult(partition: number) {
    if (!this.resultReady[partition]) throw new Error("A consensus-result partition is not ready.");
    const handle = this.fileHandles[partition]; if (!handle) throw new Error("A consensus-result scratch file is missing.");
    const file = await handle.getFile();
    if (file.size !== this.lengths[partition]) throw new Error("An external consensus-result partition has an unexpected size.");
    const bytes = new Uint8Array(await file.arrayBuffer());
    await this.directory.removeEntry(`partition-${partition}.bin`);
    this.fileHandles[partition] = undefined; this.lengths[partition] = 0; this.resultReady[partition] = false;
    return bytes;
  }

  sizes() { return [...this.lengths]; }

  async close() {
    if (this.closed) return; this.closed = true;
    const writers = this.writers.filter((writer): writer is ExternalScratchWritable => Boolean(writer));
    this.pending.forEach((chunks) => { chunks.length = 0; }); this.pendingBytes.fill(0);
    if (!this.sealed) await Promise.allSettled(writers.map((writer) => writer.abort("webPORPID analysis finished or stopped")));
    try { await this.root.removeEntry(this.directoryName, { recursive: true }); }
    catch (cause) { throw new Error(`The temporary webPORPID scratch directory could not be removed automatically. Remove “${this.directoryName}” manually. ${cause instanceof Error ? cause.message : String(cause)}`, { cause }); }
  }
}

export class PartitionStore {
  private backend: PartitionBackend;
  readonly persistent: boolean;
  readonly mode: PartitionStorageMode;
  readonly storage: BrowserStorageInformation;
  private tail: Promise<void> = Promise.resolve();
  private closed = false;
  private readonly seen: bigint[];
  private readonly partitionCount: number;
  private readonly maximumReadsPerSample: number;
  private readonly compactionIntervalBytes: number;
  private bytesSinceCompaction = 0;
  private observedRecords = 0;
  private admittedRecords = 0;
  private bypassedRecords = 0;
  private writtenBytes = 0;
  private reclaimedBytes = 0;
  private compactions = 0;

  private constructor(
    backend: PartitionBackend,
    mode: PartitionStorageMode,
    storage: BrowserStorageInformation,
    options: PartitionStoreOptions,
  ) {
    this.backend = backend; this.mode = mode; this.persistent = mode === "opfs"; this.storage = storage;
    this.seen = Array.from({ length: options.sampleCount ?? 0 }, () => 0n);
    this.partitionCount = backend.sizes().length;
    this.maximumReadsPerSample = Math.max(0, Math.floor(options.maximumReadsPerSample ?? 0));
    this.compactionIntervalBytes = Math.max(1, Math.floor(options.compactionIntervalBytes ?? DEFAULT_COMPACTION_INTERVAL));
  }

  static async create(count: number, options: PartitionStoreOptions = {}) {
    if (options.externalDirectory)
      return new PartitionStore(await ExternalDirectoryBackend.create(count, options.externalDirectory), "external-directory", {}, options);
    const storage = typeof navigator === "undefined" ? undefined : navigator.storage as StorageManager & { getDirectory?: () => Promise<FileSystemDirectoryHandle> };
    if (options.requestPersistence && storage) try { await storage.persist(); } catch { /* best effort */ }
    const information = await browserStorageInformation(storage);
    try {
      if (typeof storage?.getDirectory === "function" && typeof FileSystemFileHandle !== "undefined" && "createSyncAccessHandle" in FileSystemFileHandle.prototype)
        return new PartitionStore(await OpfsBackend.create(count), "opfs", information, options);
    } catch { /* quota/security failure falls back to a bounded in-memory spool */ }
    return new PartitionStore(new MemoryBackend(count), "memory", information, options);
  }

  private enqueue(operation: () => Promise<void>) {
    if (this.closed) return Promise.reject(new Error("The temporary partition store is closed."));
    const result = this.tail.then(operation); this.tail = result.catch(() => undefined); return result;
  }

  private provisionalCutoffs() { return this.seen.map((count) => adaptiveSpoolCutoff(count, this.maximumReadsPerSample)); }

  async appendFrames(bytes: Uint8Array) {
    return this.enqueue(async () => {
      const grouped = new Map<number, Uint8Array[]>(); let addedBytes = 0;
      for (const frame of decodeBrowserFrames(bytes)) {
        if (frame.partition < 0 || frame.partition >= this.partitionCount) throw new Error("A spool routing frame names an invalid partition.");
        const header = parseSpoolRecordHeader(frame.record);
        if (header.recordLength !== frame.record.byteLength) throw new Error("A spool routing frame contains more than one record.");
        if (header.sample >= this.seen.length) throw new Error("A spool record names an invalid sample.");
        const observed = ++this.seen[header.sample]; this.observedRecords++;
        if (!admitSpoolRecord(header, observed, this.maximumReadsPerSample)) { this.bypassedRecords++; continue; }
        const chunks = grouped.get(frame.partition) ?? []; chunks.push(frame.record); grouped.set(frame.partition, chunks);
        this.admittedRecords++; addedBytes += frame.record.byteLength;
      }
      await Promise.all([...grouped].map(([partition, chunks]) => this.backend.append(partition, chunks)));
      this.writtenBytes += addedBytes; this.bytesSinceCompaction += addedBytes;
      if (this.backend.compactable && this.maximumReadsPerSample > 0 && this.bytesSinceCompaction >= this.compactionIntervalBytes
          && this.seen.some((count) => count > BigInt(this.maximumReadsPerSample))) {
        const reclaimed = await this.backend.compact(this.provisionalCutoffs());
        this.reclaimedBytes += reclaimed; this.compactions++; this.bytesSinceCompaction = 0;
      }
    });
  }

  async compact(cutoffs: readonly bigint[]) {
    return this.enqueue(async () => {
      if (!this.backend.compactable) return;
      const reclaimed = await this.backend.compact(cutoffs);
      this.reclaimedBytes += reclaimed; this.compactions++; this.bytesSinceCompaction = 0;
    });
  }

  async seal() { return this.enqueue(() => this.backend.seal()); }

  async readSelected(partition: number, cutoffs: readonly bigint[]) { await this.tail; return this.backend.readSelected(partition, cutoffs); }
  async replaceWithResult(partition: number, bytes: Uint8Array) { await this.tail; return this.backend.replaceWithResult(partition, bytes); }
  async readResult(partition: number) { await this.tail; return this.backend.readResult(partition); }
  sizes() { return this.backend.sizes(); }
  statistics(): PartitionStoreStatistics {
    return {
      observedRecords: this.observedRecords, observedPerSample: [...this.seen],
      admittedRecords: this.admittedRecords, bypassedRecords: this.bypassedRecords,
      currentBytes: this.backend.sizes().reduce((sum, value) => sum + value, 0), writtenBytes: this.writtenBytes,
      reclaimedBytes: this.reclaimedBytes, compactions: this.compactions,
    };
  }
  async close() {
    if (this.closed) return; await this.tail; this.closed = true; await this.backend.close();
  }
}
