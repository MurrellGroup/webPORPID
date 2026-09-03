import { gzipSync } from "fflate";
import { encodeResult, exportComponent, safeDatasetName, type ExportKind } from "./result-file.ts";
import { buildOverviewExports } from "./overview-export.ts";
import { reportPlotSvgs } from "./report-svg.ts";
import type { ResultBundle } from "./types";

export const SAMPLE_EXPORT_KINDS: readonly ExportKind[] = [
  "consensus-fasta", "passed-consensus-fasta", "rejected-consensus-fasta", "trimmed-nt-fasta", "trimmed-aa-fasta",
  "family-csv", "low-agreement-csv", "contamination-csv", "postproc-csv", "apobec-csv", "collapse-csv",
  "nucleotide-alignment", "protein-alignment", "newick", "uncollapsed-nucleotide-alignment",
  "uncollapsed-protein-alignment", "uncollapsed-newick", "functional-nucleotide-alignment",
  "functional-protein-alignment", "functional-newick",
];

interface TarEntry { path: string; bytes: Uint8Array }
export interface ExportArchiveOptions { includeStaticTreeHighlighters?: boolean }
const encoder = new TextEncoder();

function writeAscii(target: Uint8Array, offset: number, length: number, value: string) {
  const bytes = encoder.encode(value);
  if (bytes.byteLength > length) throw new Error(`Archive path or header value is too long: ${value}`);
  target.set(bytes, offset);
}

function writeOctal(target: Uint8Array, offset: number, length: number, value: number) {
  const text = Math.max(0, Math.floor(value)).toString(8).padStart(length - 1, "0");
  writeAscii(target, offset, length - 1, text.slice(-(length - 1)));
  target[offset + length - 1] = 0;
}

function tarHeader(path: string, size: number, modifiedSeconds: number) {
  const header = new Uint8Array(512);
  if (encoder.encode(path).byteLength <= 100) writeAscii(header, 0, 100, path);
  else {
    const split = path.lastIndexOf("/");
    if (split < 1) throw new Error(`Archive path is too long: ${path}`);
    const prefix = path.slice(0, split), name = path.slice(split + 1);
    writeAscii(header, 0, 100, name); writeAscii(header, 345, 155, prefix);
  }
  writeOctal(header, 100, 8, 0o644);
  writeOctal(header, 108, 8, 0);
  writeOctal(header, 116, 8, 0);
  writeOctal(header, 124, 12, size);
  writeOctal(header, 136, 12, modifiedSeconds);
  header.fill(32, 148, 156);
  header[156] = 48;
  writeAscii(header, 257, 6, "ustar\0");
  writeAscii(header, 263, 2, "00");
  writeAscii(header, 265, 32, "webporpid");
  writeAscii(header, 297, 32, "webporpid");
  const checksum = header.reduce((sum, byte) => sum + byte, 0).toString(8).padStart(6, "0").slice(-6);
  writeAscii(header, 148, 6, checksum); header[154] = 0; header[155] = 32;
  return header;
}

function makeTar(entries: readonly TarEntry[], modifiedSeconds: number) {
  const total = entries.reduce((sum, entry) => sum + 512 + Math.ceil(entry.bytes.byteLength / 512) * 512, 1024);
  const output = new Uint8Array(total); let offset = 0;
  for (const entry of entries) {
    output.set(tarHeader(entry.path, entry.bytes.byteLength, modifiedSeconds), offset); offset += 512;
    output.set(entry.bytes, offset); offset += Math.ceil(entry.bytes.byteLength / 512) * 512;
  }
  return output;
}

function uniqueSampleDirectories(samples: readonly string[]) {
  const used = new Set<string>();
  return samples.map((sample, index) => {
    const base = (safeDatasetName(sample) || `sample-${index + 1}`).slice(0, 56); let candidate = base, suffix = 2;
    while (used.has(candidate.toLowerCase())) candidate = `${base.slice(0, 50)}-${suffix++}`;
    used.add(candidate.toLowerCase()); return candidate;
  });
}

/**
 * Build a standard ustar archive and gzip the complete container. Individual
 * donor outputs are generated through exportComponent, so Export all cannot
 * silently diverge from the corresponding one-at-a-time downloads.
 */
export function buildExportArchive(bundle: ResultBundle, options: ExportArchiveOptions = {}): Uint8Array {
  const includeStaticTreeHighlighters = options.includeStaticTreeHighlighters ?? true;
  const entries: TarEntry[] = [], dataset = safeDatasetName(bundle.config.dataset).slice(0, 72);
  const samples = bundle.config.samples.map((sample) => sample.name), directories = uniqueSampleDirectories(samples);
  const manifest = [
    `webPORPID export bundle for ${bundle.config.dataset}`,
    "",
    `Each sample directory contains the same sample-filtered components offered by the individual Export control, plus static SVG jitter and scatter plots${includeStaticTreeHighlighters ? ", and balanced-width phylogram + modal-highlighter figures" : ""}.`,
    "Every filename inside a sample directory begins with that sample's filesystem-safe name and an underscore.",
    "The complete editable .webporpid project and run log are stored at the archive root.",
    "cross-sample-overview/ contains the overview, parameter, provenance, stage-status, timing, and input-mapping tables.",
    "An empty alignment or Newick file means that component was unavailable or phylogeny inference was deferred.",
    "",
    ...samples.map((sample, index) => `${directories[index]}/\t${sample}`),
    "",
  ].join("\n");
  entries.push({ path: "README.txt", bytes: encoder.encode(manifest) });
  entries.push({ path: `${dataset}.webporpid`, bytes: encodeResult(bundle) });
  entries.push({ path: "run.log.txt", bytes: encoder.encode(exportComponent(bundle, "log").text) });
  for (const overview of buildOverviewExports(bundle))
    entries.push({ path: `cross-sample-overview/${overview.name}`, bytes: encoder.encode(overview.text) });
  samples.forEach((sample, index) => {
    const prefix = safeDatasetName(sample).slice(0, 72);
    const sampleFile = (extension: string) => `${prefix.slice(0, Math.max(1, 99 - extension.length))}_${extension}`;
    for (const kind of SAMPLE_EXPORT_KINDS) {
      const component = exportComponent(bundle, kind, sample);
      entries.push({ path: `${directories[index]}/${sampleFile(component.extension)}`, bytes: encoder.encode(component.text) });
    }
    for (const plot of reportPlotSvgs(bundle, sample, includeStaticTreeHighlighters))
      entries.push({ path: `${directories[index]}/${sampleFile(plot.extension)}`, bytes: encoder.encode(plot.text) });
  });
  const timestamp = Math.max(0, Math.floor((Date.parse(bundle.provenance.createdUtc) || 0) / 1000));
  return gzipSync(makeTar(entries, timestamp), { level: 6 });
}
