import type { ConsensusRecord, ContaminationCall, PipelineConfig } from "./types";

interface VectorSummary { squaredNorm: number; rootNorm: number; total: number }
interface SparseKmers extends VectorSummary { kind: "sparse"; codes: Uint16Array; counts: Uint32Array }
interface DenseKmers extends VectorSummary { kind: "dense"; values: Float64Array }
type KmerVector = SparseKmers | DenseKmers;
interface DatabaseEntry { label: string; sample?: string; selfGroup?: string; vector: KmerVector }
interface DatabaseView { self: DatabaseEntry[]; nonself: DatabaseEntry[] }
interface DatabaseIndex {
  entries: DatabaseEntry[];
  offsets: Uint32Array;
  entryIds: Uint16Array | Uint32Array;
  values: Float64Array;
  dots: Float64Array;
}
interface Cluster { center: KmerVector; memberCount: number }
interface KmerScratch { counts: Uint32Array; touched: Uint16Array }

/** Samples without donor metadata remain separate biological self groups. */
export function selfGroupsBySample(config: PipelineConfig): Map<string, string> {
  return new Map(config.samples.map((sample) => [sample.name,
    sample.donorId?.trim() ? `donor:${sample.donorId.trim()}` : `sample:${sample.name}`]));
}

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

// Values are the same A/C/G/T choice ordering used by the original string-based
// resolver. Its literal fallback emits an unresolved `N` for unknown symbols,
// which resets rolling six-mer parsing rather than randomly resolving the base.
const IUPAC_CODES: Readonly<Record<number, readonly number[]>> = {
  65: [0], 67: [1], 71: [2], 84: [3], 85: [3], 82: [0, 2], 89: [1, 3],
  83: [1, 2], 87: [0, 3], 75: [2, 3], 77: [0, 1], 66: [1, 2, 3],
  68: [0, 2, 3], 72: [0, 1, 3], 86: [0, 1, 2], 78: [0, 1, 2, 3],
};
const INVALID_BASE = [-1] as const;

function createKmerScratch(): KmerScratch {
  return { counts: new Uint32Array(4096), touched: new Uint16Array(4096) };
}

/** Resolve IUPAC symbols and count six-mers without an intermediate string or a full-bin scan. */
function kmers(sequence: string, seed: number, scratch: KmerScratch): SparseKmers {
  const dense = scratch.counts, touched = scratch.touched;
  let state = seed >>> 0, code = 0, valid = 0, touchedCount = 0, squaredNorm = 0, total = 0;
  for (let offset = 0; offset < sequence.length; offset++) {
    let character = sequence.charCodeAt(offset);
    if (character === 45) continue;
    if (character >= 97 && character <= 122) character -= 32;
    state ^= state << 13; state ^= state >>> 17; state ^= state << 5;
    const choices = IUPAC_CODES[character] ?? INVALID_BASE;
    const base = choices[(state >>> 0) % choices.length];
    if (base < 0) { code = 0; valid = 0; continue; }
    code = ((code << 2) | base) & 4095; valid++;
    if (valid < 6) continue;
    const previous = dense[code];
    if (previous === 0) touched[touchedCount++] = code;
    dense[code] = previous + 1; squaredNorm += 2 * previous + 1; total++;
  }
  const codes = touched.slice(0, touchedCount); codes.sort();
  const counts = new Uint32Array(touchedCount);
  for (let index = 0; index < touchedCount; index++) {
    const observedCode = codes[index]; counts[index] = dense[observedCode]; dense[observedCode] = 0;
  }
  return { kind: "sparse", codes, counts, squaredNorm, rootNorm: Math.sqrt(squaredNorm), total };
}

function hasRandomAmbiguity(sequence: string) {
  for (let index = 0; index < sequence.length; index++) {
    let character = sequence.charCodeAt(index); if (character >= 97 && character <= 122) character -= 32;
    if ((IUPAC_CODES[character]?.length ?? 1) > 1) return true;
  }
  return false;
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

function distanceFromDot(left: KmerVector, right: KmerVector, dotProduct: number) {
  const squared = Math.max(0, left.squaredNorm + right.squaredNorm - 2 * dotProduct);
  const total = left.total + right.total;
  return total === 0 ? 0 : squared / (6 * total);
}

/** Reverse-triangle lower bound on the exact normalized squared distance. */
function distanceLowerBound(left: KmerVector, right: KmerVector) {
  const total = left.total + right.total; if (total === 0) return 0;
  const difference = left.rootNorm - right.rootNorm;
  return difference * difference / (6 * total);
}

const emptyCenter = (): DenseKmers => ({ kind: "dense", values: new Float64Array(4096), squaredNorm: 0, rootNorm: 0, total: 0 });

function meanAll(vectors: SparseKmers[]): KmerVector {
  if (vectors.length === 1) return vectors[0];
  if (!vectors.length) return emptyCenter();
  const values = new Float64Array(4096), touched: number[] = [];
  for (const vector of vectors) {
    for (let index = 0; index < vector.codes.length; index++) {
      const code = vector.codes[index]; if (values[code] === 0) touched.push(code);
      values[code] += vector.counts[index];
    }
  }
  let squaredNorm = 0, total = 0;
  touched.sort((left, right) => left - right);
  for (const code of touched) { values[code] /= vectors.length; squaredNorm += values[code] ** 2; total += values[code]; }
  return { kind: "dense", values, squaredNorm, rootNorm: Math.sqrt(squaredNorm), total };
}

function updateCenters(vectors: SparseKmers[], assignments: Int32Array, centerCount: number): KmerVector[] {
  const counts = new Int32Array(centerCount), first = new Int32Array(centerCount); first.fill(-1);
  for (let point = 0; point < assignments.length; point++) {
    const cluster = assignments[point]; counts[cluster]++; if (first[cluster] < 0) first[cluster] = point;
  }
  const values: Array<Float64Array | undefined> = Array(centerCount), touched: Array<number[] | undefined> = Array(centerCount);
  for (let cluster = 0; cluster < centerCount; cluster++) if (counts[cluster] > 1) {
    values[cluster] = new Float64Array(4096); touched[cluster] = [];
  }
  for (let point = 0; point < assignments.length; point++) {
    const cluster = assignments[point], target = values[cluster]; if (!target) continue;
    const observed = vectors[point], observedCodes = touched[cluster]!;
    for (let index = 0; index < observed.codes.length; index++) {
      const code = observed.codes[index]; if (target[code] === 0) observedCodes.push(code);
      target[code] += observed.counts[index];
    }
  }
  return Array.from({ length: centerCount }, (_, cluster): KmerVector => {
    if (counts[cluster] === 0) return emptyCenter();
    if (counts[cluster] === 1) return vectors[first[cluster]];
    const target = values[cluster]!, observedCodes = touched[cluster]!; observedCodes.sort((left, right) => left - right);
    let squaredNorm = 0, total = 0;
    for (const code of observedCodes) { target[code] /= counts[cluster]; squaredNorm += target[code] ** 2; total += target[code]; }
    return { kind: "dense", values: target, squaredNorm, rootNorm: Math.sqrt(squaredNorm), total };
  });
}

function assignPointIndexed(vector: SparseKmers, centers: KmerVector[], index: DatabaseIndex | undefined,
  indexedCount: number, radius: number) {
  let best = 0, bestDistance = Number.POSITIVE_INFINITY;
  if (index) {
    scoreDatabase(vector, index);
    for (let cluster = 0; cluster < indexedCount; cluster++) {
      const candidate = distanceFromDot(vector, centers[cluster], index.dots[cluster]);
      if (candidate < bestDistance) { best = cluster; bestDistance = candidate; }
    }
  }
  for (let cluster = indexedCount; cluster < centers.length; cluster++) {
    if (distanceLowerBound(centers[cluster], vector) >= bestDistance) continue;
    const candidate = distance(centers[cluster], vector);
    if (candidate < bestDistance) { best = cluster; bestDistance = candidate; }
  }
  return bestDistance > radius ? -1 : best;
}

const CENTER_INDEX_REBUILD_TAIL = 128;

function centerIndex(centers: KmerVector[]) {
  return buildDatabaseIndex(centers.map((vector, index) => ({ label: String(index), vector })));
}

function clusterMembers(assignments: Int32Array, centers: KmerVector[]): Cluster[] {
  const counts = new Int32Array(centers.length);
  for (let point = 0; point < assignments.length; point++) counts[assignments[point]]++;
  return centers.map((center, index) => ({ center, memberCount: counts[index] })).filter((cluster) => cluster.memberCount);
}

/** Fixed scientific/runtime trade-off requested for run-derived contamination signatures. */
export const CONTAMINATION_CLUSTER_PASSES = 3;

function dpMeans(vectors: SparseKmers[], radius: number): Cluster[] {
  if (!vectors.length) return [];
  let centers: KmerVector[] = [vectors[0]];
  const assignments = new Int32Array(vectors.length); assignments.fill(-1);
  for (let iteration = 0; iteration < CONTAMINATION_CLUSTER_PASSES; iteration++) {
    let changed = false, indexed = centerIndex(centers), indexedCount = indexed ? centers.length : 0, canRebuild = Boolean(indexed);
    for (let point = 0; point < vectors.length; point++) {
      if (canRebuild && centers.length - indexedCount >= CENTER_INDEX_REBUILD_TAIL) {
        const replacement = centerIndex(centers);
        if (replacement) { indexed = replacement; indexedCount = centers.length; } else canRebuild = false;
      }
      let assigned = assignPointIndexed(vectors[point], centers, indexed, indexedCount, radius);
      if (assigned < 0) { centers.push(vectors[point]); assigned = centers.length - 1; }
      if (assignments[point] !== assigned) { assignments[point] = assigned; changed = true; }
    }
    centers = updateCenters(vectors, assignments, centers.length);
    if (!changed) break;
  }
  return clusterMembers(assignments, centers);
}

const abortError = () => new DOMException("Contamination checks skipped.", "AbortError");
const yieldToMessages = () => new Promise<void>((resolvePromise) => setTimeout(resolvePromise, 0));

class CooperativeScheduler {
  private lastYield = performance.now();
  private readonly signal?: AbortSignal;
  private readonly budgetMs: number;
  constructor(signal?: AbortSignal, budgetMs = 32) { this.signal = signal; this.budgetMs = budgetMs; }
  due() { return performance.now() - this.lastYield >= this.budgetMs; }
  async yieldNow() {
    if (this.signal?.aborted) throw abortError();
    await yieldToMessages(); this.lastYield = performance.now();
    if (this.signal?.aborted) throw abortError();
  }
  check() { if (this.signal?.aborted) throw abortError(); }
}

async function dpMeansAsync(
  vectors: SparseKmers[], radius: number, scheduler: CooperativeScheduler,
  onIteration?: (iteration: number, assigned: number, total: number) => void,
): Promise<Cluster[]> {
  if (!vectors.length) return [];
  let centers: KmerVector[] = [vectors[0]];
  const assignments = new Int32Array(vectors.length); assignments.fill(-1);
  for (let iteration = 0; iteration < CONTAMINATION_CLUSTER_PASSES; iteration++) {
    let changed = false, indexed = centerIndex(centers), indexedCount = indexed ? centers.length : 0, canRebuild = Boolean(indexed);
    for (let point = 0; point < vectors.length; point++) {
      scheduler.check();
      if (canRebuild && centers.length - indexedCount >= CENTER_INDEX_REBUILD_TAIL) {
        const replacement = centerIndex(centers);
        if (replacement) { indexed = replacement; indexedCount = centers.length; } else canRebuild = false;
        if (scheduler.due()) { onIteration?.(iteration, point, vectors.length); await scheduler.yieldNow(); }
      }
      let assigned = assignPointIndexed(vectors[point], centers, indexed, indexedCount, radius);
      if (assigned < 0) { centers.push(vectors[point]); assigned = centers.length - 1; }
      if (assignments[point] !== assigned) { assignments[point] = assigned; changed = true; }
      if ((point & 31) === 31 && scheduler.due()) {
        onIteration?.(iteration, point + 1, vectors.length); await scheduler.yieldNow();
      }
    }
    centers = updateCenters(vectors, assignments, centers.length);
    onIteration?.(iteration, vectors.length, vectors.length);
    if (scheduler.due()) await scheduler.yieldNow();
    if (!changed) break;
  }
  return clusterMembers(assignments, centers);
}

function databaseView(database: DatabaseEntry[], selfGroup: string): DatabaseView {
  const self: DatabaseEntry[] = [], nonself: DatabaseEntry[] = [];
  for (const entry of database) (entry.selfGroup === selfGroup ? self : nonself).push(entry);
  return { self, nonself };
}

function nearestNonself(vector: KmerVector, entries: DatabaseEntry[], threshold: number) {
  let nearest: DatabaseEntry | undefined, nearestDistance = threshold;
  for (const entry of entries) {
    if (distanceLowerBound(vector, entry.vector) >= nearestDistance) continue;
    const candidate = distance(vector, entry.vector);
    if (candidate < nearestDistance) { nearest = entry; nearestDistance = candidate; }
  }
  return nearest ? { label: nearest.label, distance: nearestDistance } : undefined;
}

function nearestPrimary(vector: KmerVector, view: DatabaseView, threshold: number) {
  const nonself = nearestNonself(vector, view.nonself, threshold); if (!nonself) return undefined;
  let hasCloseSelf = false;
  for (const entry of view.self) {
    if (distanceLowerBound(vector, entry.vector) > threshold) continue;
    if (distance(vector, entry.vector) <= threshold) { hasCloseSelf = true; break; }
  }
  return { ...nonself, discard: !hasCloseSelf };
}

// A posting index turns query × database × observed-kmer scans into one pass
// over only shared six-mers.  The cap bounds the temporary index on unusually
// fragmented runs; those fall back to the exact sparse scan below.
const MAX_INDEX_POSTINGS = 4_000_000;

function buildDatabaseIndex(entries: DatabaseEntry[], postingLimit = MAX_INDEX_POSTINGS): DatabaseIndex | undefined {
  const binSizes = new Uint32Array(4096); let postingCount = 0;
  for (const entry of entries) {
    if (entry.vector.kind === "sparse") {
      postingCount += entry.vector.codes.length;
      if (postingCount > postingLimit) return undefined;
      for (const code of entry.vector.codes) binSizes[code]++;
    } else {
      for (let code = 0; code < 4096; code++) if (entry.vector.values[code] !== 0) {
        binSizes[code]++; if (++postingCount > postingLimit) return undefined;
      }
    }
  }
  const offsets = new Uint32Array(4097);
  for (let code = 0; code < 4096; code++) offsets[code + 1] = offsets[code] + binSizes[code];
  const cursor = offsets.slice(0, 4096);
  const entryIds = entries.length <= 0xffff ? new Uint16Array(postingCount) : new Uint32Array(postingCount);
  const values = new Float64Array(postingCount);
  const add = (code: number, entry: number, value: number) => {
    const position = cursor[code]++; entryIds[position] = entry; values[position] = value;
  };
  for (let entry = 0; entry < entries.length; entry++) {
    const vector = entries[entry].vector;
    if (vector.kind === "sparse") {
      for (let offset = 0; offset < vector.codes.length; offset++) add(vector.codes[offset], entry, vector.counts[offset]);
    } else {
      for (let code = 0; code < 4096; code++) if (vector.values[code] !== 0) add(code, entry, vector.values[code]);
    }
  }
  return { entries, offsets, entryIds, values, dots: new Float64Array(entries.length) };
}

function scoreDatabase(vector: SparseKmers, index: DatabaseIndex) {
  index.dots.fill(0);
  for (let queryOffset = 0; queryOffset < vector.codes.length; queryOffset++) {
    const code = vector.codes[queryOffset], queryCount = vector.counts[queryOffset];
    for (let position = index.offsets[code]; position < index.offsets[code + 1]; position++)
      index.dots[index.entryIds[position]] += queryCount * index.values[position];
  }
}

function nearestNonselfIndexed(vector: SparseKmers, selfGroup: string, index: DatabaseIndex, threshold: number) {
  scoreDatabase(vector, index);
  let nearest: DatabaseEntry | undefined, nearestDistance = threshold;
  for (let entry = 0; entry < index.entries.length; entry++) {
    const candidateEntry = index.entries[entry]; if (candidateEntry.selfGroup === selfGroup) continue;
    const candidate = distanceFromDot(vector, candidateEntry.vector, index.dots[entry]);
    if (candidate < nearestDistance) { nearest = candidateEntry; nearestDistance = candidate; }
  }
  return nearest ? { label: nearest.label, distance: nearestDistance } : undefined;
}

function nearestPrimaryIndexed(vector: SparseKmers, selfGroup: string, index: DatabaseIndex, threshold: number) {
  scoreDatabase(vector, index);
  let nearest: DatabaseEntry | undefined, nearestDistance = threshold, hasCloseSelf = false;
  for (let entry = 0; entry < index.entries.length; entry++) {
    const candidateEntry = index.entries[entry], candidate = distanceFromDot(vector, candidateEntry.vector, index.dots[entry]);
    if (candidateEntry.selfGroup === selfGroup) { if (candidate <= threshold) hasCloseSelf = true; }
    else if (candidate < nearestDistance) { nearest = candidateEntry; nearestDistance = candidate; }
  }
  return nearest ? { label: nearest.label, distance: nearestDistance, discard: !hasCloseSelf } : undefined;
}

/** One stored contamination decision per family, preferring a primary call over a suspect-only call. */
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

interface PreparedVectors {
  bySample: Map<string, SparseKmers[]>;
  byConsensus: Array<SparseKmers | undefined>;
}

function prepareConsensusVectors(consensuses: ConsensusRecord[], config: PipelineConfig, scratch: KmerScratch): PreparedVectors {
  const bySample = new Map(config.samples.map((sample) => [sample.name, [] as SparseKmers[]]));
  const byConsensus: Array<SparseKmers | undefined> = Array(consensuses.length);
  for (let globalIndex = 0; globalIndex < consensuses.length; globalIndex++) {
    const record = consensuses[globalIndex], sampleVectors = bySample.get(record.sample);
    if (!sampleVectors) continue;
    const vector = kmers(record.sequence, hash(record.id, sampleVectors.length), scratch);
    sampleVectors.push(vector); byConsensus[globalIndex] = vector;
  }
  return { bySample, byConsensus };
}

function addSampleDatabases(sampleName: string, selfGroup: string, vectors: SparseKmers[], config: PipelineConfig,
  primary: DatabaseEntry[], suspect: DatabaseEntry[], clusters: Cluster[]) {
  clusters.forEach((cluster, index) => {
    const proportion = cluster.memberCount / vectors.length;
    const percent = Math.round(1000 * proportion) / 10;
    const entry = { label: `${sampleName}_cluster${index + 1}_${percent}%`, sample: sampleName, selfGroup, vector: cluster.center };
    suspect.push(entry); if (proportion > config.parameters.contaminationProportionThreshold) primary.push(entry);
  });
  const all = { label: `${sampleName}_All`, sample: sampleName, selfGroup, vector: meanAll(vectors) };
  primary.push(all); suspect.push(all);
}

function classifyPrepared(consensuses: ConsensusRecord[], prepared: PreparedVectors, primary: DatabaseEntry[], suspect: DatabaseEntry[],
  config: PipelineConfig, scratch: KmerScratch): ContaminationCall[] {
  const threshold = config.parameters.contaminationDistanceThreshold, output: ContaminationCall[] = [];
  const primaryViews = new Map<string, DatabaseView>(), suspectViews = new Map<string, DatabaseView>();
  const primaryIndex = buildDatabaseIndex(primary), suspectIndex = buildDatabaseIndex(suspect);
  const selfGroups = selfGroupsBySample(config);
  for (const sample of new Set(consensuses.map((record) => record.sample))) {
    const selfGroup = selfGroups.get(sample) ?? `sample:${sample}`;
    primaryViews.set(sample, databaseView(primary, selfGroup)); suspectViews.set(sample, databaseView(suspect, selfGroup));
  }
  for (let index = 0; index < consensuses.length; index++) {
    const record = consensuses[index];
    const vector = !hasRandomAmbiguity(record.sequence) && prepared.byConsensus[index]
      ? prepared.byConsensus[index]! : kmers(record.sequence, hash(record.id, index), scratch);
    const selfGroup = selfGroups.get(record.sample) ?? `sample:${record.sample}`;
    const call = primaryIndex ? nearestPrimaryIndexed(vector, selfGroup, primaryIndex, threshold)
      : nearestPrimary(vector, primaryViews.get(record.sample)!, threshold);
    if (call) output.push({ sample: record.sample, sequenceId: record.id, nearestNonselfVariant: call.label,
      nearestNonselfDistance: call.distance, flagged: true, discarded: call.discard, suspectOnly: false });
    else {
      const possible = suspectIndex ? nearestNonselfIndexed(vector, selfGroup, suspectIndex, threshold)
        : nearestNonself(vector, suspectViews.get(record.sample)!.nonself, threshold);
      if (possible) output.push({ sample: record.sample, sequenceId: record.id, nearestNonselfVariant: possible.label,
        nearestNonselfDistance: possible.distance, flagged: true, discarded: false, suspectOnly: true });
    }
  }
  return deduplicateContaminationCalls(output);
}

export function classifyContamination(consensuses: ConsensusRecord[], config: PipelineConfig): ContaminationCall[] {
  if (!config.parameters.contaminationFilter) return [];
  const scratch = createKmerScratch();
  const panel: DatabaseEntry[] = config.contaminationPanelSequences.map((record, index) => ({
    label: record.name, vector: kmers(record.sequence, hash(record.name, index), scratch),
  }));
  const prepared = prepareConsensusVectors(consensuses, config, scratch), primary: DatabaseEntry[] = [], suspect: DatabaseEntry[] = [];
  const selfGroups = selfGroupsBySample(config);
  for (const sample of config.samples) {
    const vectors = prepared.bySample.get(sample.name) ?? []; if (!vectors.length) continue;
    addSampleDatabases(sample.name, selfGroups.get(sample.name)!, vectors, config, primary, suspect, dpMeans(vectors, config.parameters.contaminationClusterThreshold));
  }
  primary.push(...panel); suspect.push(...panel);
  return classifyPrepared(consensuses, prepared, primary, suspect, config, scratch);
}

/** Cancellable classifier with exact sparse kernels and time-budgeted cooperative yields. */
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
  const scheduler = new CooperativeScheduler(signal), scratch = createKmerScratch(); scheduler.check();
  const panel: DatabaseEntry[] = [];
  for (let index = 0; index < config.contaminationPanelSequences.length; index++) {
    const record = config.contaminationPanelSequences[index];
    panel.push({ label: record.name, vector: kmers(record.sequence, hash(record.name, index), scratch) });
    if (scheduler.due()) {
      report(.02 * (index + 1) / Math.max(1, config.contaminationPanelSequences.length),
        `Prepared ${index + 1} of ${config.contaminationPanelSequences.length} contamination-reference sequences`, index + 1,
        config.contaminationPanelSequences.length, "reference-vectors");
      await scheduler.yieldNow();
    }
  }
  report(.02, `Prepared ${panel.length.toLocaleString()} contamination-reference sequences`, panel.length, panel.length, "reference-vectors");

  const bySample = new Map(config.samples.map((sample) => [sample.name, [] as SparseKmers[]]));
  const byConsensus: Array<SparseKmers | undefined> = Array(consensuses.length);
  for (let index = 0; index < consensuses.length; index++) {
    const record = consensuses[index], sampleVectors = bySample.get(record.sample);
    if (sampleVectors) {
      const vector = kmers(record.sequence, hash(record.id, sampleVectors.length), scratch);
      sampleVectors.push(vector); byConsensus[index] = vector;
    }
    if (scheduler.due()) {
      report(.02 + .18 * (index + 1) / Math.max(1, consensuses.length),
        `Prepared sequence signatures for ${index + 1} of ${consensuses.length.toLocaleString()} consensus sequences`,
        index + 1, consensuses.length, "sample-vectors");
      await scheduler.yieldNow();
    }
  }
  const prepared = { bySample, byConsensus };
  report(.2, `Prepared all ${consensuses.length.toLocaleString()} consensus sequence signatures`, consensuses.length, consensuses.length, "sample-vectors");

  const primary: DatabaseEntry[] = [], suspect: DatabaseEntry[] = [];
  let clusteredVectors = 0;
  const selfGroups = selfGroupsBySample(config);
  for (const sample of config.samples) {
    const vectors = bySample.get(sample.name) ?? [], before = clusteredVectors; let clusterCount = 0;
    if (vectors.length) {
      const clusters = await dpMeansAsync(vectors, config.parameters.contaminationClusterThreshold, scheduler,
        (iteration, assigned, total) => report(.2 + .35 * (before + vectors.length * Math.min(.95, (iteration + assigned / Math.max(1, total)) / CONTAMINATION_CLUSTER_PASSES)) / Math.max(1, consensuses.length),
          `Clustering ${sample.name}: iteration ${iteration + 1}, ${assigned.toLocaleString()} of ${total.toLocaleString()} signatures assigned`,
          before + assigned, consensuses.length, "clustering"));
      clusterCount = clusters.length;
      addSampleDatabases(sample.name, selfGroups.get(sample.name)!, vectors, config, primary, suspect, clusters);
    }
    clusteredVectors += vectors.length;
    report(.2 + .35 * clusteredVectors / Math.max(1, consensuses.length),
      `Clustered ${sample.name} into ${clusterCount.toLocaleString()} run signatures; processed ${clusteredVectors.toLocaleString()} of ${consensuses.length.toLocaleString()} consensus sequences`,
      clusteredVectors, consensuses.length, "clustering");
    if (scheduler.due()) await scheduler.yieldNow();
  }
  primary.push(...panel); suspect.push(...panel);

  const output: ContaminationCall[] = [], threshold = config.parameters.contaminationDistanceThreshold;
  const primaryViews = new Map<string, DatabaseView>(), suspectViews = new Map<string, DatabaseView>();
  const primaryIndex = buildDatabaseIndex(primary), suspectIndex = buildDatabaseIndex(suspect);
  for (const sample of new Set(consensuses.map((record) => record.sample))) {
    const selfGroup = selfGroups.get(sample) ?? `sample:${sample}`;
    primaryViews.set(sample, databaseView(primary, selfGroup)); suspectViews.set(sample, databaseView(suspect, selfGroup));
  }
  for (let index = 0; index < consensuses.length; index++) {
    scheduler.check();
    const record = consensuses[index];
    const vector = !hasRandomAmbiguity(record.sequence) && prepared.byConsensus[index]
      ? prepared.byConsensus[index]! : kmers(record.sequence, hash(record.id, index), scratch);
    const selfGroup = selfGroups.get(record.sample) ?? `sample:${record.sample}`;
    const call = primaryIndex ? nearestPrimaryIndexed(vector, selfGroup, primaryIndex, threshold)
      : nearestPrimary(vector, primaryViews.get(record.sample)!, threshold);
    if (call) output.push({ sample: record.sample, sequenceId: record.id, nearestNonselfVariant: call.label,
      nearestNonselfDistance: call.distance, flagged: true, discarded: call.discard, suspectOnly: false });
    else {
      const possible = suspectIndex ? nearestNonselfIndexed(vector, selfGroup, suspectIndex, threshold)
        : nearestNonself(vector, suspectViews.get(record.sample)!.nonself, threshold);
      if (possible) output.push({ sample: record.sample, sequenceId: record.id, nearestNonselfVariant: possible.label,
        nearestNonselfDistance: possible.distance, flagged: true, discarded: false, suspectOnly: true });
    }
    if (scheduler.due() || index + 1 === consensuses.length) {
      report(.55 + .45 * (index + 1) / Math.max(1, consensuses.length),
        `Classified ${index + 1} of ${consensuses.length.toLocaleString()} consensus sequences; ${output.length.toLocaleString()} currently flagged`,
        index + 1, consensuses.length, "classification");
      if (scheduler.due()) await scheduler.yieldNow();
    }
  }
  report(1, `Contamination checks complete for ${consensuses.length.toLocaleString()} consensus sequences`, consensuses.length,
    consensuses.length, "classification");
  return deduplicateContaminationCalls(output);
}
