export const SPOOL_HEADER_BYTES = 24;
const MAX_SPOOL_RECORD_BYTES = 256 * 1024 * 1024;

export interface SpoolRecordHeader {
  recordLength: number;
  sample: number;
  samplingHash: bigint;
}

export function parseSpoolRecordHeader(bytes: Uint8Array): SpoolRecordHeader {
  if (bytes.byteLength < SPOOL_HEADER_BYTES) throw new Error("A temporary spool record has a truncated header.");
  const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  const bodyLength = view.getUint32(0, true), sample = view.getUint16(4, true);
  const umiLength = view.getUint16(6, true), nameLength = view.getUint32(8, true), sequenceLength = view.getUint32(12, true);
  const recordLength = bodyLength + 4;
  const expectedBody = 20 + umiLength + nameLength + sequenceLength * 2;
  if (bodyLength < 20 || recordLength > MAX_SPOOL_RECORD_BYTES || bodyLength !== expectedBody || !Number.isSafeInteger(recordLength))
    throw new Error("A temporary spool record has inconsistent lengths.");
  return { recordLength, sample, samplingHash: view.getBigUint64(16, true) };
}

export function selectedSpoolRecord(header: SpoolRecordHeader, cutoffs: readonly bigint[]) {
  return header.sample >= cutoffs.length || header.samplingHash <= cutoffs[header.sample];
}

export function concatenateSpoolRecords(records: readonly Uint8Array[]) {
  const length = records.reduce((sum, record) => sum + record.byteLength, 0);
  const output = new Uint8Array(length); let offset = 0;
  for (const record of records) { output.set(record, offset); offset += record.byteLength; }
  return output;
}

export function selectSpoolChunks(records: readonly Uint8Array[], cutoffs: readonly bigint[]) {
  const selected: Uint8Array[] = [];
  for (const record of records) {
    const header = parseSpoolRecordHeader(record);
    if (header.recordLength !== record.byteLength) throw new Error("A temporary spool chunk contains more than one record.");
    if (selectedSpoolRecord(header, cutoffs)) selected.push(record);
  }
  return concatenateSpoolRecords(selected);
}
