import type { ConsensusRecord, ContaminationCall, PipelineConfig } from "./types";

interface SparseKmers { kind: "sparse"; codes: Uint16Array; counts: Uint32Array; squaredNorm: number; total: number }
interface DenseKmers { kind: "dense"; values: Float64Array; squaredNorm: number; total: number }
type KmerVector = SparseKmers | DenseKmers;
interface DatabaseEntry { label: string; sample?: string; vector: KmerVector }
interface Cluster { center: KmerVector; members: number[] }

export interface ContaminationProgress {
  fraction: number;
  detail: string;
  completed: number;
  total: number;
  phase: "reference-vectors" | "sample-vectors" | "clustering" | "classification";
}

function hash(text: string, ordinal = 0) {
  let value = (0x811c9dc5 ^ ordinal) >>> 0;
  for (let index = 0; index < text.length; index++) value = Math.imul(value ^ text.charCodeAt(index), 0x01000193) >>> 0;
  return value || 1;
}

const IUPAC: Record<string, string> = {
  A: "A", C: "C", G: "G", T: "T", U: "T", R: "AG", Y: "CT", S: "CG", W: "AT",
  K: "GT", M: "AC", B: "CGT", D: "AGT", H: "ACT", V: "ACG", N: "ACGT",
};

function resolve(sequence: string, seed: number) {
  let state = seed >>> 0;
  return [...sequence.replaceAll("-", "").toUpperCase()].map((base) => {
    const choices = IUPAC[base] ?? "N"; state ^= state << 13; state ^= state >>> 17; state ^= state << 5;
    return choices[(state >>> 0) % choices.length];
  }).join("");
}

function kmers(sequence: string) {
  const dense = new Uint32Array(4096), index: Record<string, number> = { A: 0, C: 1, G: 2, T: 3 };
  let code = 0, valid = 0;
  for (const base of sequence) {
    const value = index[base]; if (value === undefined) { code = 0; valid = 0; continue; }
    code = ((code << 2) | value) & 4095; valid++; if (valid >= 6) dense[code]++;
  }
  let nonzero = 0, squaredNorm = 0, total = 0;
  for (const count of dense) if (count) { nonzero++; squaredNorm += count * count; total += count; }
  const codes = new Uint16Array(nonzero), counts = new Uint32Array(nonzero); let offset = 0;
  dense.forEach((count, kmer) => { if (count) { codes[offset] = kmer; counts[offset] = count; offset++; } });
  return { kind: "sparse", codes, counts, squaredNorm, total } satisfies SparseKmers;
}

function dot(left: KmerVector, right: KmerVector) {
  if (left.kind === "dense" && right.kind === "dense") {
    let output = 0; for (let index = 0; index < 4096; index++) output += left.values[index] * right.values[index]; return output;
  }
  if (left.kind === "dense") return dot(right, left);
  if (right.kind === "dense") {
    let output = 0; for (let index = 0; index < left.codes.length; index++) output += left.counts[index] * right.values[left.codes[index]]; return output;
  }
  let output = 0, a = 0, b = 0;
  while (a < left.codes.length && b < right.codes.length) {
    if (left.codes[a] < right.codes[b]) a++;
    else if (right.codes[b] < left.codes[a]) b++;
    else { output += left.counts[a] * right.counts[b]; a++; b++; }
  }
  return output;
}

function distance(left: KmerVector, right: KmerVector) {
  const squared = Math.max(0, left.squaredNorm + right.squaredNorm - 2 * dot(left, right));
  const total = left.total + right.total;
  return total === 0 ? 0 : squared / (6 * total);
}

function mean(vectors: KmerVector[], members: number[]): KmerVector {
  if (members.length === 1) return vectors[members[0]];
  const values = new Float64Array(4096);
  for (const member of members) {
    const vector = vectors[member];
    if (vector.kind === "dense") for (let index = 0; index < 4096; index++) values[index] += vector.values[index];
    else for (let index = 0; index < vector.codes.length; index++) values[vector.codes[index]] += vector.counts[index];
  }
  let squaredNorm = 0, total = 0;
  for (let index = 0; index < 4096; index++) { values[index] /= members.length; squaredNorm += values[index] ** 2; total += values[index]; }
  return { kind: "dense", values, squaredNorm, total };
}

function dpMeans(vectors: KmerVector[], radius: number): Cluster[] {
  if (!vectors.length) return [];
  let centers = [vectors[0]], assignments = Array(vectors.length).fill(-1);
  for (let iteration = 0; iteration < 30; iteration++) {
    const previous = assignments.join(",");
    for (let point = 0; point < vectors.length; point++) {
      let best = 0, bestDistance = Number.POSITIVE_INFINITY;
      for (let cluster = 0; cluster < centers.length; cluster++) {
        const candidate = distance(centers[cluster], vectors[point]); if (candidate < bestDistance) { best = cluster; bestDistance = candidate; }
      }
      if (bestDistance > radius) { centers.push(vectors[point]); assignments[point] = centers.length - 1; } else assignments[point] = best;
    }
    const members = Array.from({ length: centers.length }, () => [] as number[]);
    assignments.forEach((cluster, point) => members[cluster].push(point));
    centers = members.map((indices) => indices.length ? mean(vectors, indices)
      : { kind: "dense", values: new Float64Array(4096), squaredNorm: 0, total: 0 });
    if (assignments.join(",") === previous) break;
  }
  const members = Array.from({ length: centers.length }, () => [] as number[]);
  assignments.forEach((cluster, point) => members[cluster].push(point));
  return centers.map((center, index) => ({ center, members: members[index] })).filter((cluster) => cluster.members.length);
}

const abortError = () => new DOMException("Contamination checks skipped.", "AbortError");
const yieldToMessages = () => new Promise<void>((resolvePromise) => setTimeout(resolvePromise, 0));

async function dpMeansAsync(
  vectors: KmerVector[], radius: number, signal: AbortSignal | undefined,
  onIteration?: (iteration: number, assigned: number, total: number) => void,
): Promise<Cluster[]> {
  if (!vectors.length) return [];
  let centers = [vectors[0]], assignments = Array(vectors.length).fill(-1);
  for (let iteration = 0; iteration < 30; iteration++) {
    const previous = assignments.join(",");
    for (let point = 0; point < vectors.length; point++) {
      if (signal?.aborted) throw abortError();
      let best = 0, bestDistance = Number.POSITIVE_INFINITY;
      for (let cluster = 0; cluster < centers.length; cluster++) {
        const candidate = distance(centers[cluster], vectors[point]);
        if (candidate < bestDistance) { best = cluster; bestDistance = candidate; }
      }
      if (bestDistance > radius) { centers.push(vectors[point]); assignments[point] = centers.length - 1; }
      else assignments[point] = best;
      if ((point & 63) === 63) { onIteration?.(iteration, point + 1, vectors.length); await yieldToMessages(); }
    }
    const members = Array.from({ length: centers.length }, () => [] as number[]);
    assignments.forEach((cluster, point) => members[cluster].push(point));
    centers = members.map((indices) => indices.length ? mean(vectors, indices)
      : { kind: "dense", values: new Float64Array(4096), squaredNorm: 0, total: 0 });
    onIteration?.(iteration, vectors.length, vectors.length);
    await yieldToMessages();
    if (assignments.join(",") === previous) break;
  }
  const members = Array.from({ length: centers.length }, () => [] as number[]);
  assignments.forEach((cluster, point) => members[cluster].push(point));
  return centers.map((center, index) => ({ center, members: members[index] })).filter((cluster) => cluster.members.length);
}

function nearest(sample: string, vector: KmerVector, database: DatabaseEntry[], threshold: number) {
  let closestSelf = Number.POSITIVE_INFINITY, nonself: DatabaseEntry | undefined, nonselfDistance = Number.POSITIVE_INFINITY;
  for (const entry of database) {
    const value = distance(vector, entry.vector);
    if (entry.sample === sample) closestSelf = Math.min(closestSelf, value);
    else if (value < nonselfDistance) { nonself = entry; nonselfDistance = value; }
  }
  return nonself && nonselfDistance < threshold ? { label: nonself.label, distance: nonselfDistance, discard: closestSelf > threshold } : undefined;
}

/**
 * One consensus must have one displayed/stored contamination decision. Older
 * result files can contain both the primary call and the wider suspect-pass
 * call; prefer the primary decision and use the suspect call only when there
 * was no primary match.
 */
export function deduplicateContaminationCalls(calls: readonly ContaminationCall[]): ContaminationCall[] {
  const bySequence = new Map<string, ContaminationCall>();
  for (const call of calls) {
    const key = `${call.sample}\0${call.sequenceId}`, previous = bySequence.get(key);
    if (!previous || (previous.suspectOnly && !call.suspectOnly)
      || (previous.suspectOnly === call.suspectOnly && call.nearestNonselfDistance < previous.nearestNonselfDistance))
      bySequence.set(key, { ...call });
  }
  return [...bySequence.values()].sort((a, b) => a.sample.localeCompare(b.sample)
    || a.nearestNonselfDistance - b.nearestNonselfDistance || a.sequenceId.localeCompare(b.sequenceId));
}

export function classifyContamination(consensuses: ConsensusRecord[], config: PipelineConfig): ContaminationCall[] {
  if (!config.parameters.contaminationFilter) return [];
  const panel: DatabaseEntry[] = config.contaminationPanelSequences.map((record, index) => ({
    label: record.name, vector: kmers(resolve(record.sequence, hash(record.name, index))),
  }));
  const primary: DatabaseEntry[] = [], suspect: DatabaseEntry[] = [];
  for (const sample of config.samples) {
    const records = consensuses.filter((record) => record.sample === sample.name);
    if (!records.length) continue;
    const vectors = records.map((record, index) => kmers(resolve(record.sequence, hash(record.id, index))));
    const clusters = dpMeans(vectors, config.parameters.contaminationClusterThreshold);
    clusters.forEach((cluster, index) => {
      const percent = Math.round(1000 * cluster.members.length / vectors.length) / 10;
      const entry = { label: `${sample.name}_cluster${index + 1}_${percent}%`, sample: sample.name, vector: cluster.center };
      suspect.push(entry); if (cluster.members.length / vectors.length > config.parameters.contaminationProportionThreshold) primary.push(entry);
    });
    const all = { label: `${sample.name}_All`, sample: sample.name, vector: mean(vectors, vectors.map((_, index) => index)) };
    primary.push(all); suspect.push(all);
  }
  // The reference workflow concatenates the run-derived cluster database before the external
  // panel. Keeping this order preserves argmin tie selection and report labels.
  primary.push(...panel); suspect.push(...panel);
  const output: ContaminationCall[] = [], threshold = config.parameters.contaminationDistanceThreshold;
  for (const [index, record] of consensuses.entries()) {
    const vector = kmers(resolve(record.sequence, hash(record.id, index))), call = nearest(record.sample, vector, primary, threshold);
    if (call) output.push({ sample: record.sample, sequenceId: record.id, nearestNonselfVariant: call.label,
      nearestNonselfDistance: call.distance, flagged: true, discarded: call.discard, suspectOnly: false });
    const possible = nearest(record.sample, vector, suspect, threshold);
    if (!call && possible) output.push({ sample: record.sample, sequenceId: record.id, nearestNonselfVariant: possible.label,
      nearestNonselfDistance: possible.distance, flagged: true, discarded: false, suspectOnly: true });
  }
  return deduplicateContaminationCalls(output);
}

/**
 * Cancellable browser implementation of the same classifier. Work is yielded
 * in bounded chunks so progress and Skip messages remain responsive even for
 * large run-wide databases; all tie-breaking and thresholds match the
 * synchronous implementation above.
 */
export async function classifyContaminationAsync(
  consensuses: ConsensusRecord[], config: PipelineConfig, signal?: AbortSignal,
  onProgress?: (progress: ContaminationProgress) => void,
): Promise<ContaminationCall[]> {
  const report = (fraction: number, detail: string, completed: number, total: number, phase: ContaminationProgress["phase"]) =>
    onProgress?.({ fraction: Math.max(0, Math.min(1, fraction)), detail, completed, total, phase });
  if (!config.parameters.contaminationFilter) {
    report(1, "Contamination filtering is disabled in the configuration; no sequences were excluded.", consensuses.length, consensuses.length, "classification");
    return [];
  }
  if (signal?.aborted) throw abortError();
  const panel: DatabaseEntry[] = [];
  for (const [index, record] of config.contaminationPanelSequences.entries()) {
    panel.push({ label: record.name, vector: kmers(resolve(record.sequence, hash(record.name, index))) });
    if ((index & 31) === 31) {
      report(.02 * (index + 1) / Math.max(1, config.contaminationPanelSequences.length),
        `Prepared ${index + 1} of ${config.contaminationPanelSequences.length} contamination-reference sequences`, index + 1,
        config.contaminationPanelSequences.length, "reference-vectors");
      await yieldToMessages(); if (signal?.aborted) throw abortError();
    }
  }
  report(.02, `Prepared ${panel.length.toLocaleString()} contamination-reference sequences`, panel.length, panel.length, "reference-vectors");

  const primary: DatabaseEntry[] = [], suspect: DatabaseEntry[] = [];
  const sampleVectors = new Map<string, KmerVector[]>();
  let vectorsPrepared = 0;
  for (const sample of config.samples) {
    const records = consensuses.filter((record) => record.sample === sample.name);
    const vectors: KmerVector[] = [];
    for (const [index, record] of records.entries()) {
      vectors.push(kmers(resolve(record.sequence, hash(record.id, index)))); vectorsPrepared++;
      if ((vectorsPrepared & 63) === 0) {
        report(.02 + .18 * vectorsPrepared / Math.max(1, consensuses.length),
          `Prepared sequence signatures for ${vectorsPrepared.toLocaleString()} of ${consensuses.length.toLocaleString()} consensus sequences`,
          vectorsPrepared, consensuses.length, "sample-vectors");
        await yieldToMessages(); if (signal?.aborted) throw abortError();
      }
    }
    sampleVectors.set(sample.name, vectors);
  }
  report(.2, `Prepared all ${vectorsPrepared.toLocaleString()} consensus sequence signatures`, vectorsPrepared, consensuses.length, "sample-vectors");

  let clusteredSamples = 0;
  for (const sample of config.samples) {
    const vectors = sampleVectors.get(sample.name) ?? [];
    if (vectors.length) {
      const sampleIndex = clusteredSamples;
      const clusters = await dpMeansAsync(vectors, config.parameters.contaminationClusterThreshold, signal,
        (iteration, assigned, total) => report(.2 + .35 * (sampleIndex + Math.min(.95, (iteration + assigned / Math.max(1, total)) / 30)) / Math.max(1, config.samples.length),
          `Clustering ${sample.name}: iteration ${iteration + 1}, ${assigned.toLocaleString()} of ${total.toLocaleString()} sequence signatures assigned`,
          assigned, total, "clustering"));
      clusters.forEach((cluster, index) => {
        const percent = Math.round(1000 * cluster.members.length / vectors.length) / 10;
        const entry = { label: `${sample.name}_cluster${index + 1}_${percent}%`, sample: sample.name, vector: cluster.center };
        suspect.push(entry);
        if (cluster.members.length / vectors.length > config.parameters.contaminationProportionThreshold) primary.push(entry);
      });
      const all = { label: `${sample.name}_All`, sample: sample.name, vector: mean(vectors, vectors.map((_, index) => index)) };
      primary.push(all); suspect.push(all);
    }
    clusteredSamples++;
    report(.2 + .35 * clusteredSamples / Math.max(1, config.samples.length),
      `Built run-wide contamination signatures for ${clusteredSamples} of ${config.samples.length} samples`, clusteredSamples,
      config.samples.length, "clustering");
    await yieldToMessages(); if (signal?.aborted) throw abortError();
  }
  primary.push(...panel); suspect.push(...panel);

  const output: ContaminationCall[] = [], threshold = config.parameters.contaminationDistanceThreshold;
  for (const [index, record] of consensuses.entries()) {
    if (signal?.aborted) throw abortError();
    const vector = kmers(resolve(record.sequence, hash(record.id, index))), call = nearest(record.sample, vector, primary, threshold);
    if (call) output.push({ sample: record.sample, sequenceId: record.id, nearestNonselfVariant: call.label,
      nearestNonselfDistance: call.distance, flagged: true, discarded: call.discard, suspectOnly: false });
    const possible = nearest(record.sample, vector, suspect, threshold);
    if (!call && possible) output.push({ sample: record.sample, sequenceId: record.id, nearestNonselfVariant: possible.label,
      nearestNonselfDistance: possible.distance, flagged: true, discarded: false, suspectOnly: true });
    if ((index & 31) === 31 || index + 1 === consensuses.length) {
      report(.55 + .45 * (index + 1) / Math.max(1, consensuses.length),
        `Classified ${index + 1} of ${consensuses.length} consensus sequences; ${output.length.toLocaleString()} currently flagged`,
        index + 1, consensuses.length, "classification");
      await yieldToMessages();
    }
  }
  report(1, `Contamination checks complete for ${consensuses.length.toLocaleString()} consensus sequences`, consensuses.length,
    consensuses.length, "classification");
  return deduplicateContaminationCalls(output);
}
