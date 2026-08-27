import { decodeBrowserFrames } from "./wasm-runtime";
import { concatenateSpoolRecords, parseSpoolRecordHeader, selectedSpoolRecord, selectSpoolChunks, SPOOL_HEADER_BYTES } from "./spool-record";

interface PartitionBackend {
  append(partition: number, chunks: Uint8Array[]): Promise<void>;
  readSelected(partition: number, cutoffs: readonly bigint[]): Promise<Uint8Array>;
  sizes(): number[];
  close(): Promise<void>;
}

class MemoryBackend implements PartitionBackend {
  private parts: Uint8Array[][];
  private lengths: number[];
  constructor(count: number, private maximumBytes = 512 * 1024 * 1024) { this.parts = Array.from({ length: count }, () => []); this.lengths = Array(count).fill(0); }
  async append(partition: number, chunks: Uint8Array[]) {
    const added = chunks.reduce((sum, chunk) => sum + chunk.byteLength, 0);
    if (this.lengths.reduce((sum, value) => sum + value, 0) + added > this.maximumBytes)
      throw new Error("Browser-local temporary storage is unavailable and the in-memory spool exceeded 512 MiB. Use a browser with OPFS support or porpid-cli.");
    this.parts[partition].push(...chunks); this.lengths[partition] += added;
  }
  async readSelected(partition: number, cutoffs: readonly bigint[]) { return selectSpoolChunks(this.parts[partition], cutoffs); }
  sizes() { return [...this.lengths]; }
  async close() { this.parts = []; this.lengths = []; }
}

class OpfsBackend implements PartitionBackend {
  private constructor(
    private root: FileSystemDirectoryHandle,
    private directoryName: string,
    private handles: FileSystemSyncAccessHandle[],
    private lengths: number[],
  ) {}

  static async create(count: number) {
    const root = await navigator.storage.getDirectory();
    const directoryName = `webporpid-${crypto.randomUUID()}`;
    const directory = await root.getDirectoryHandle(directoryName, { create: true });
    const handles: FileSystemSyncAccessHandle[] = [];
    for (let index = 0; index < count; index++) {
      const file = await directory.getFileHandle(`partition-${index}.bin`, { create: true });
      handles.push(await file.createSyncAccessHandle());
    }
    return new OpfsBackend(root, directoryName, handles, Array(count).fill(0));
  }

  async append(partition: number, chunks: Uint8Array[]) {
    const handle = this.handles[partition]; let offset = this.lengths[partition];
    for (const chunk of chunks) { const written = handle.write(chunk, { at: offset }); if (written !== chunk.byteLength) throw new Error("OPFS wrote an incomplete spool frame."); offset += written; }
    this.lengths[partition] = offset;
  }
  private readFully(handle: FileSystemSyncAccessHandle, output: Uint8Array, at: number) {
    let offset = 0;
    while (offset < output.byteLength) {
      const read = handle.read(output.subarray(offset), { at: at + offset });
      if (!read) throw new Error("OPFS returned an incomplete spool record."); offset += read;
    }
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
    for (const handle of this.handles) { handle.flush(); handle.close(); }
    this.handles = []; await this.root.removeEntry(this.directoryName, { recursive: true });
  }
}

export class PartitionStore {
  private constructor(private backend: PartitionBackend, readonly persistent: boolean) {}
  static async create(count: number) {
    try {
      const storage = typeof navigator === "undefined" ? undefined : navigator.storage as StorageManager & { getDirectory?: () => Promise<FileSystemDirectoryHandle> };
      if (typeof storage?.getDirectory === "function" && "createSyncAccessHandle" in FileSystemFileHandle.prototype)
        return new PartitionStore(await OpfsBackend.create(count), true);
    } catch { /* quota/security failure falls back to a bounded in-memory spool */ }
    return new PartitionStore(new MemoryBackend(count), false);
  }
  async appendFrames(bytes: Uint8Array) {
    const grouped = new Map<number, Uint8Array[]>();
    for (const frame of decodeBrowserFrames(bytes)) {
      const chunks = grouped.get(frame.partition) ?? []; chunks.push(frame.record); grouped.set(frame.partition, chunks);
    }
    await Promise.all([...grouped].map(([partition, chunks]) => this.backend.append(partition, chunks)));
  }
  readSelected(partition: number, cutoffs: readonly bigint[]) { return this.backend.readSelected(partition, cutoffs); }
  sizes() { return this.backend.sizes(); }
  close() { return this.backend.close(); }
}
