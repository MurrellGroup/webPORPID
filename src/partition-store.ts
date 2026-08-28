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
  append(partition: number, chunks: Uint8Array[]): Promise<void>;
  compact(cutoffs: readonly bigint[]): Promise<number>;
  readSelected(partition: number, cutoffs: readonly bigint[]): Promise<Uint8Array>;
  sizes(): number[];
  close(): Promise<void>;
}

export interface PartitionStoreOptions {
  sampleCount?: number;
  maximumReadsPerSample?: number;
  compactionIntervalBytes?: number;
  requestPersistence?: boolean;
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
  private parts: Uint8Array[][];
  private lengths: number[];
  private maximumBytes: number;
  constructor(count: number, maximumBytes = 512 * 1024 * 1024) {
    this.maximumBytes = maximumBytes; this.parts = Array.from({ length: count }, () => []); this.lengths = Array(count).fill(0);
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
  async readSelected(partition: number, cutoffs: readonly bigint[]) { return selectSpoolChunks(this.parts[partition], cutoffs); }
  sizes() { return [...this.lengths]; }
  async close() { this.parts = []; this.lengths = []; }
}

class OpfsBackend implements PartitionBackend {
  private closed = false;
  private root: FileSystemDirectoryHandle;
  private directoryName: string;
  private handles: FileSystemSyncAccessHandle[];
  private lengths: number[];
  private constructor(
    root: FileSystemDirectoryHandle,
    directoryName: string,
    handles: FileSystemSyncAccessHandle[],
    lengths: number[],
  ) { this.root = root; this.directoryName = directoryName; this.handles = handles; this.lengths = lengths; }

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

export class PartitionStore {
  private backend: PartitionBackend;
  readonly persistent: boolean;
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
    persistent: boolean,
    storage: BrowserStorageInformation,
    options: PartitionStoreOptions,
  ) {
    this.backend = backend; this.persistent = persistent; this.storage = storage;
    this.seen = Array.from({ length: options.sampleCount ?? 0 }, () => 0n);
    this.partitionCount = backend.sizes().length;
    this.maximumReadsPerSample = Math.max(0, Math.floor(options.maximumReadsPerSample ?? 0));
    this.compactionIntervalBytes = Math.max(1, Math.floor(options.compactionIntervalBytes ?? DEFAULT_COMPACTION_INTERVAL));
  }

  static async create(count: number, options: PartitionStoreOptions = {}) {
    const storage = typeof navigator === "undefined" ? undefined : navigator.storage as StorageManager & { getDirectory?: () => Promise<FileSystemDirectoryHandle> };
    if (options.requestPersistence && storage) try { await storage.persist(); } catch { /* best effort */ }
    const information = await browserStorageInformation(storage);
    try {
      if (typeof storage?.getDirectory === "function" && typeof FileSystemFileHandle !== "undefined" && "createSyncAccessHandle" in FileSystemFileHandle.prototype)
        return new PartitionStore(await OpfsBackend.create(count), true, information, options);
    } catch { /* quota/security failure falls back to a bounded in-memory spool */ }
    return new PartitionStore(new MemoryBackend(count), false, information, options);
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
      if (this.maximumReadsPerSample > 0 && this.bytesSinceCompaction >= this.compactionIntervalBytes
          && this.seen.some((count) => count > BigInt(this.maximumReadsPerSample))) {
        const reclaimed = await this.backend.compact(this.provisionalCutoffs());
        this.reclaimedBytes += reclaimed; this.compactions++; this.bytesSinceCompaction = 0;
      }
    });
  }

  async compact(cutoffs: readonly bigint[]) {
    return this.enqueue(async () => {
      const reclaimed = await this.backend.compact(cutoffs);
      this.reclaimedBytes += reclaimed; this.compactions++; this.bytesSinceCompaction = 0;
    });
  }

  async readSelected(partition: number, cutoffs: readonly bigint[]) { await this.tail; return this.backend.readSelected(partition, cutoffs); }
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
