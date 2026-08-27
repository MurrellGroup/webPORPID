const encoder = new TextEncoder();
const decoder = new TextDecoder();

export class BinaryWriter {
  private chunks: Uint8Array[] = [];
  private length = 0;

  private push(bytes: Uint8Array) { this.chunks.push(bytes); this.length += bytes.byteLength; }
  magic(value: string) { this.push(encoder.encode(value)); }
  u8(value: number) { this.push(Uint8Array.of(value & 255)); }
  u16(value: number) { const b = new Uint8Array(2); new DataView(b.buffer).setUint16(0, value, true); this.push(b); }
  i32(value: number) { const b = new Uint8Array(4); new DataView(b.buffer).setInt32(0, value, true); this.push(b); }
  u32(value: number) { const b = new Uint8Array(4); new DataView(b.buffer).setUint32(0, value, true); this.push(b); }
  u64(value: bigint) { const b = new Uint8Array(8); new DataView(b.buffer).setBigUint64(0, value, true); this.push(b); }
  f32(value: number) { const b = new Uint8Array(4); new DataView(b.buffer).setFloat32(0, value, true); this.push(b); }
  f64(value: number) { const b = new Uint8Array(8); new DataView(b.buffer).setFloat64(0, value, true); this.push(b); }
  bytes(value: Uint8Array) { this.push(value); }
  string(value: string) { const bytes = encoder.encode(value); this.u32(bytes.byteLength); this.push(bytes); }
  finish(): Uint8Array {
    const output = new Uint8Array(this.length); let offset = 0;
    for (const chunk of this.chunks) { output.set(chunk, offset); offset += chunk.byteLength; }
    return output;
  }
}

export class BinaryReader {
  private view: DataView;
  private offset = 0;
  constructor(private bytesValue: Uint8Array) { this.view = new DataView(bytesValue.buffer, bytesValue.byteOffset, bytesValue.byteLength); }
  private require(length: number) { if (this.offset + length > this.bytesValue.byteLength) throw new Error("Truncated webPORPID binary payload."); }
  magic(value: string) {
    const expected = encoder.encode(value); this.require(expected.length);
    for (let index = 0; index < expected.length; index++) if (this.bytesValue[this.offset + index] !== expected[index]) throw new Error(`Expected ${value} payload.`);
    this.offset += expected.length;
  }
  u8() { this.require(1); return this.view.getUint8(this.offset++); }
  u16() { this.require(2); const v = this.view.getUint16(this.offset, true); this.offset += 2; return v; }
  i32() { this.require(4); const v = this.view.getInt32(this.offset, true); this.offset += 4; return v; }
  u32() { this.require(4); const v = this.view.getUint32(this.offset, true); this.offset += 4; return v; }
  u64() { this.require(8); const v = this.view.getBigUint64(this.offset, true); this.offset += 8; return v; }
  f32() { this.require(4); const v = this.view.getFloat32(this.offset, true); this.offset += 4; return v; }
  f64() { this.require(8); const v = this.view.getFloat64(this.offset, true); this.offset += 8; return v; }
  bytes(length: number) { this.require(length); const v = this.bytesValue.subarray(this.offset, this.offset + length); this.offset += length; return v; }
  string() { const length = this.u32(); return decoder.decode(this.bytes(length)); }
  get done() { return this.offset === this.bytesValue.byteLength; }
  get remaining() { return this.bytesValue.byteLength - this.offset; }
}
