#!/usr/bin/env node
import { a as resolveReferenceFiles, i as parseConfigYaml, n as createFastTreeRunner, o as resultConfig, r as compileConfig } from "./chunks/direct-fasttree-BVyb7bUx.mjs";
import { a as encodeAlivibeMsaSequences, i as decodeAlivibeMsaSequences, n as createMsaRunner, r as assertAlivibeMsaResult } from "./chunks/direct-msa-DfCxIVAW.mjs";
import { a as makeCutoffValues, c as mergeStats, i as decodeFamilyModel, n as decodeConsensusOutput, o as makeCutoffs, r as decodeFamilyCounts, s as mergeFamilyCounts } from "./chunks/wasm-runtime-Dm9J7url.mjs";
import { createHash } from "node:crypto";
import { createReadStream } from "node:fs";
import { mkdir, mkdtemp, open, readFile, rm, stat, writeFile } from "node:fs/promises";
import { availableParallelism, tmpdir } from "node:os";
import { basename, dirname, isAbsolute, join, resolve } from "node:path";
import { createInterface } from "node:readline";
import { createGunzip } from "node:zlib";
import { Worker as Worker$1 } from "node:worker_threads";
//#region src/contamination.ts
function hash(text, ordinal = 0) {
	let value = (2166136261 ^ ordinal) >>> 0;
	for (let index = 0; index < text.length; index++) value = Math.imul(value ^ text.charCodeAt(index), 16777619) >>> 0;
	return value || 1;
}
const IUPAC = {
	A: "A",
	C: "C",
	G: "G",
	T: "T",
	U: "T",
	R: "AG",
	Y: "CT",
	S: "CG",
	W: "AT",
	K: "GT",
	M: "AC",
	B: "CGT",
	D: "AGT",
	H: "ACT",
	V: "ACG",
	N: "ACGT"
};
function resolve$1(sequence, seed) {
	let state = seed >>> 0;
	return [...sequence.replaceAll("-", "").toUpperCase()].map((base) => {
		const choices = IUPAC[base] ?? "N";
		state ^= state << 13;
		state ^= state >>> 17;
		state ^= state << 5;
		return choices[(state >>> 0) % choices.length];
	}).join("");
}
function kmers(sequence) {
	const dense = new Uint32Array(4096), index = {
		A: 0,
		C: 1,
		G: 2,
		T: 3
	};
	let code = 0, valid = 0;
	for (const base of sequence) {
		const value = index[base];
		if (value === void 0) {
			code = 0;
			valid = 0;
			continue;
		}
		code = (code << 2 | value) & 4095;
		valid++;
		if (valid >= 6) dense[code]++;
	}
	let nonzero = 0, squaredNorm = 0, total = 0;
	for (const count of dense) if (count) {
		nonzero++;
		squaredNorm += count * count;
		total += count;
	}
	const codes = new Uint16Array(nonzero), counts = new Uint32Array(nonzero);
	let offset = 0;
	dense.forEach((count, kmer) => {
		if (count) {
			codes[offset] = kmer;
			counts[offset] = count;
			offset++;
		}
	});
	return {
		kind: "sparse",
		codes,
		counts,
		squaredNorm,
		total
	};
}
function dot(left, right) {
	if (left.kind === "dense" && right.kind === "dense") {
		let output = 0;
		for (let index = 0; index < 4096; index++) output += left.values[index] * right.values[index];
		return output;
	}
	if (left.kind === "dense") return dot(right, left);
	if (right.kind === "dense") {
		let output = 0;
		for (let index = 0; index < left.codes.length; index++) output += left.counts[index] * right.values[left.codes[index]];
		return output;
	}
	let output = 0, a = 0, b = 0;
	while (a < left.codes.length && b < right.codes.length) if (left.codes[a] < right.codes[b]) a++;
	else if (right.codes[b] < left.codes[a]) b++;
	else {
		output += left.counts[a] * right.counts[b];
		a++;
		b++;
	}
	return output;
}
function distance(left, right) {
	const squared = Math.max(0, left.squaredNorm + right.squaredNorm - 2 * dot(left, right));
	const total = left.total + right.total;
	return total === 0 ? 0 : squared / (6 * total);
}
function mean(vectors, members) {
	if (members.length === 1) return vectors[members[0]];
	const values = new Float64Array(4096);
	for (const member of members) {
		const vector = vectors[member];
		if (vector.kind === "dense") for (let index = 0; index < 4096; index++) values[index] += vector.values[index];
		else for (let index = 0; index < vector.codes.length; index++) values[vector.codes[index]] += vector.counts[index];
	}
	let squaredNorm = 0, total = 0;
	for (let index = 0; index < 4096; index++) {
		values[index] /= members.length;
		squaredNorm += values[index] ** 2;
		total += values[index];
	}
	return {
		kind: "dense",
		values,
		squaredNorm,
		total
	};
}
function dpMeans(vectors, radius) {
	if (!vectors.length) return [];
	let centers = [vectors[0]], assignments = Array(vectors.length).fill(-1);
	for (let iteration = 0; iteration < 30; iteration++) {
		const previous = assignments.join(",");
		for (let point = 0; point < vectors.length; point++) {
			let best = 0, bestDistance = Number.POSITIVE_INFINITY;
			for (let cluster = 0; cluster < centers.length; cluster++) {
				const candidate = distance(centers[cluster], vectors[point]);
				if (candidate < bestDistance) {
					best = cluster;
					bestDistance = candidate;
				}
			}
			if (bestDistance > radius) {
				centers.push(vectors[point]);
				assignments[point] = centers.length - 1;
			} else assignments[point] = best;
		}
		const members = Array.from({ length: centers.length }, () => []);
		assignments.forEach((cluster, point) => members[cluster].push(point));
		centers = members.map((indices) => indices.length ? mean(vectors, indices) : {
			kind: "dense",
			values: new Float64Array(4096),
			squaredNorm: 0,
			total: 0
		});
		if (assignments.join(",") === previous) break;
	}
	const members = Array.from({ length: centers.length }, () => []);
	assignments.forEach((cluster, point) => members[cluster].push(point));
	return centers.map((center, index) => ({
		center,
		members: members[index]
	})).filter((cluster) => cluster.members.length);
}
function nearest(sample, vector, database, threshold) {
	let closestSelf = Number.POSITIVE_INFINITY, nonself, nonselfDistance = Number.POSITIVE_INFINITY;
	for (const entry of database) {
		const value = distance(vector, entry.vector);
		if (entry.sample === sample) closestSelf = Math.min(closestSelf, value);
		else if (value < nonselfDistance) {
			nonself = entry;
			nonselfDistance = value;
		}
	}
	return nonself && nonselfDistance < threshold ? {
		label: nonself.label,
		distance: nonselfDistance,
		discard: closestSelf > threshold
	} : void 0;
}
/**
* One consensus must have one displayed/stored contamination decision. Older
* result files can contain both the primary call and the wider suspect-pass
* call; prefer the primary decision and use the suspect call only when there
* was no primary match.
*/
function deduplicateContaminationCalls(calls) {
	const bySequence = /* @__PURE__ */ new Map();
	for (const call of calls) {
		const key = `${call.sample}\0${call.sequenceId}`, previous = bySequence.get(key);
		if (!previous || previous.suspectOnly && !call.suspectOnly || previous.suspectOnly === call.suspectOnly && call.nearestNonselfDistance < previous.nearestNonselfDistance) bySequence.set(key, { ...call });
	}
	return [...bySequence.values()].sort((a, b) => a.sample.localeCompare(b.sample) || a.nearestNonselfDistance - b.nearestNonselfDistance || a.sequenceId.localeCompare(b.sequenceId));
}
function classifyContamination(consensuses, config) {
	if (!config.parameters.contaminationFilter) return [];
	const panel = config.contaminationPanelSequences.map((record, index) => ({
		label: record.name,
		vector: kmers(resolve$1(record.sequence, hash(record.name, index)))
	}));
	const primary = [], suspect = [];
	for (const sample of config.samples) {
		const records = consensuses.filter((record) => record.sample === sample.name);
		if (!records.length) continue;
		const vectors = records.map((record, index) => kmers(resolve$1(record.sequence, hash(record.id, index))));
		dpMeans(vectors, config.parameters.contaminationClusterThreshold).forEach((cluster, index) => {
			const percent = Math.round(1e3 * cluster.members.length / vectors.length) / 10;
			const entry = {
				label: `${sample.name}_cluster${index + 1}_${percent}%`,
				sample: sample.name,
				vector: cluster.center
			};
			suspect.push(entry);
			if (cluster.members.length / vectors.length > config.parameters.contaminationProportionThreshold) primary.push(entry);
		});
		const all = {
			label: `${sample.name}_All`,
			sample: sample.name,
			vector: mean(vectors, vectors.map((_, index) => index))
		};
		primary.push(all);
		suspect.push(all);
	}
	primary.push(...panel);
	suspect.push(...panel);
	const output = [], threshold = config.parameters.contaminationDistanceThreshold;
	for (const [index, record] of consensuses.entries()) {
		const vector = kmers(resolve$1(record.sequence, hash(record.id, index))), call = nearest(record.sample, vector, primary, threshold);
		if (call) output.push({
			sample: record.sample,
			sequenceId: record.id,
			nearestNonselfVariant: call.label,
			nearestNonselfDistance: call.distance,
			flagged: true,
			discarded: call.discard,
			suspectOnly: false
		});
		const possible = nearest(record.sample, vector, suspect, threshold);
		if (!call && possible) output.push({
			sample: record.sample,
			sequenceId: record.id,
			nearestNonselfVariant: possible.label,
			nearestNonselfDistance: possible.distance,
			flagged: true,
			discarded: false,
			suspectOnly: true
		});
	}
	return deduplicateContaminationCalls(output);
}
//#endregion
//#region src/alivibe-msa-runtime.ts
async function runAlivibeMsa(sequences, signal, iterations = 3, scoringMode = "nucleotide") {
	if (sequences.length < 2) return [...sequences];
	const inputSequences = sequences.map(String), input = encodeAlivibeMsaSequences(inputSequences);
	const worker = new Worker(new URL("./alivibe-msa-worker.ts", import.meta.url), { type: "module" });
	return new Promise((resolve, reject) => {
		const abort = () => {
			worker.terminate();
			reject(new DOMException("MSA alignment cancelled.", "AbortError"));
		};
		if (signal?.aborted) return abort();
		signal?.addEventListener("abort", abort, { once: true });
		const finish = () => {
			signal?.removeEventListener("abort", abort);
			worker.terminate();
		};
		worker.onmessage = (event) => {
			if (event.data.type === "error") {
				finish();
				reject(new Error(event.data.message || "Alivibe MSA failed."));
				return;
			}
			try {
				const aligned = decodeAlivibeMsaSequences(event.data.result ?? /* @__PURE__ */ new ArrayBuffer(0));
				assertAlivibeMsaResult(inputSequences, aligned);
				finish();
				resolve(aligned);
			} catch (error) {
				finish();
				reject(error);
			}
		};
		worker.onerror = (event) => {
			finish();
			reject(new Error(event.message || "The Alivibe MSA worker stopped unexpectedly."));
		};
		worker.postMessage({
			input,
			iterations,
			scoringMode
		}, [input]);
	});
}
//#endregion
//#region src/alignment-utils.ts
const CODONS$1 = {
	TTT: "F",
	TTC: "F",
	TTA: "L",
	TTG: "L",
	TCT: "S",
	TCC: "S",
	TCA: "S",
	TCG: "S",
	TAT: "Y",
	TAC: "Y",
	TAA: "*",
	TAG: "*",
	TGT: "C",
	TGC: "C",
	TGA: "*",
	TGG: "W",
	CTT: "L",
	CTC: "L",
	CTA: "L",
	CTG: "L",
	CCT: "P",
	CCC: "P",
	CCA: "P",
	CCG: "P",
	CAT: "H",
	CAC: "H",
	CAA: "Q",
	CAG: "Q",
	CGT: "R",
	CGC: "R",
	CGA: "R",
	CGG: "R",
	ATT: "I",
	ATC: "I",
	ATA: "I",
	ATG: "M",
	ACT: "T",
	ACC: "T",
	ACA: "T",
	ACG: "T",
	AAT: "N",
	AAC: "N",
	AAA: "K",
	AAG: "K",
	AGT: "S",
	AGC: "S",
	AGA: "R",
	AGG: "R",
	GTT: "V",
	GTC: "V",
	GTA: "V",
	GTG: "V",
	GCT: "A",
	GCC: "A",
	GCA: "A",
	GCG: "A",
	GAT: "D",
	GAC: "D",
	GAA: "E",
	GAG: "E",
	GGT: "G",
	GGC: "G",
	GGA: "G",
	GGG: "G"
};
function exactFasta$1(records) {
	return records.map((record) => `>${record.name}\n${record.sequence}\n`).join("");
}
function parseAlignmentFasta(source) {
	const records = [];
	let name = "", sequence = "";
	for (const rawLine of source.split(/\r?\n/)) {
		const line = rawLine.trim();
		if (!line || line.startsWith(";")) continue;
		if (line.startsWith(">")) {
			if (name) records.push({
				name,
				sequence: sequence.toUpperCase()
			});
			name = line.slice(1).trim() || `sequence_${records.length + 1}`;
			sequence = "";
		} else {
			if (!name) throw new Error("FASTA sequence data appeared before its first header.");
			sequence += line.replace(/\s/g, "");
		}
	}
	if (name) records.push({
		name,
		sequence: sequence.toUpperCase()
	});
	if (!records.length) throw new Error("The FASTA file contains no records.");
	return records;
}
function assertFastaAlphabet(text) {
	for (const line of text.split(/\r?\n/)) {
		if (!line || line.startsWith(">")) continue;
		const unsupported = line.replace(/\s/g, "").match(/[^ACGTUNRYKMSWBDHV.\-]/i)?.[0];
		if (unsupported) throw new Error(`The alignment contains the unsupported nucleotide character ${JSON.stringify(unsupported)}.`);
	}
}
function fnv1a64(value) {
	let hash = 14695981039346656037n;
	const mask = 18446744073709551615n;
	for (const byte of new TextEncoder().encode(value)) {
		hash ^= BigInt(byte);
		hash = hash * 1099511628211n & mask;
	}
	return hash.toString(16).padStart(16, "0");
}
function inspectAlignment(text, minimumRows = 2) {
	assertFastaAlphabet(text);
	const records = parseAlignmentFasta(text.split(/\r?\n/).map((line) => line.startsWith(">") ? line : line.replace(/U/gi, "T")).join("\n"));
	if (records.length < minimumRows) throw new Error(`The alignment must contain at least ${minimumRows} sequences.`);
	const columns = records[0]?.sequence.length ?? 0;
	if (!columns || records.some((record) => record.sequence.length !== columns)) throw new Error("Every alignment record must have the same non-zero aligned length.");
	const names = /* @__PURE__ */ new Set();
	for (const record of records) {
		if (names.has(record.name)) throw new Error(`The alignment contains duplicate identifier ${record.name}.`);
		names.add(record.name);
	}
	const fasta = exactFasta$1(records);
	return {
		fasta,
		records,
		rows: records.length,
		columns,
		fingerprint: fnv1a64(fasta)
	};
}
function ungapped(sequence) {
	return sequence.replaceAll("-", "").replaceAll("U", "T");
}
/** Exact unit-cost Levenshtein decomposition with deterministic diagonal ties. */
function nucleotideEditCounts(original, replacement) {
	if (original === replacement) return {
		substitutedNucleotides: 0,
		insertedNucleotides: 0,
		removedNucleotides: 0
	};
	const width = replacement.length + 1;
	let previousCost = new Uint32Array(width), previousSubstitutions = new Uint32Array(width), previousInsertions = new Uint32Array(width), previousDeletions = new Uint32Array(width);
	let currentCost = new Uint32Array(width), currentSubstitutions = new Uint32Array(width), currentInsertions = new Uint32Array(width), currentDeletions = new Uint32Array(width);
	for (let column = 0; column < width; column += 1) {
		previousCost[column] = column;
		previousInsertions[column] = column;
	}
	for (let row = 1; row <= original.length; row += 1) {
		currentCost[0] = row;
		currentSubstitutions[0] = 0;
		currentInsertions[0] = 0;
		currentDeletions[0] = row;
		for (let column = 1; column < width; column += 1) {
			const mismatch = original[row - 1] === replacement[column - 1] ? 0 : 1;
			let cost = previousCost[column - 1] + mismatch;
			let substitutions = previousSubstitutions[column - 1] + mismatch;
			let insertions = previousInsertions[column - 1], deletions = previousDeletions[column - 1];
			const deletionCost = previousCost[column] + 1;
			if (deletionCost < cost) {
				cost = deletionCost;
				substitutions = previousSubstitutions[column];
				insertions = previousInsertions[column];
				deletions = previousDeletions[column] + 1;
			}
			const insertionCost = currentCost[column - 1] + 1;
			if (insertionCost < cost) {
				cost = insertionCost;
				substitutions = currentSubstitutions[column - 1];
				insertions = currentInsertions[column - 1] + 1;
				deletions = currentDeletions[column - 1];
			}
			currentCost[column] = cost;
			currentSubstitutions[column] = substitutions;
			currentInsertions[column] = insertions;
			currentDeletions[column] = deletions;
		}
		[previousCost, currentCost] = [currentCost, previousCost];
		[previousSubstitutions, currentSubstitutions] = [currentSubstitutions, previousSubstitutions];
		[previousInsertions, currentInsertions] = [currentInsertions, previousInsertions];
		[previousDeletions, currentDeletions] = [currentDeletions, previousDeletions];
	}
	const column = replacement.length;
	return {
		substitutedNucleotides: previousSubstitutions[column],
		insertedNucleotides: previousInsertions[column],
		removedNucleotides: previousDeletions[column]
	};
}
/**
* Validate the FASTA shape without restricting biological edits. Substitutions,
* insertions, deletions, renamed rows and removed rows are all retained in the
* edit audit so the UI can warn before saving them.
*/
function validateCorrectedAlignment(currentText, correctedText) {
	const current = inspectAlignment(currentText, 1);
	const corrected = inspectAlignment(correctedText, 1);
	const currentByName = new Map(current.records.map((record) => [record.name, record.sequence]));
	const correctedNames = new Set(corrected.records.map((record) => record.name));
	const removedRows = current.records.filter((record) => !correctedNames.has(record.name)).map((record) => record.name);
	const addedRows = corrected.records.filter((record) => !currentByName.has(record.name)).map((record) => record.name);
	const changedRows = [], rowChanges = [];
	let removedNucleotides = 0, insertedNucleotides = 0, substitutedNucleotides = 0;
	for (const record of corrected.records) {
		const originalAligned = currentByName.get(record.name);
		if (originalAligned == null) continue;
		const original = ungapped(originalAligned);
		const replacement = ungapped(record.sequence);
		if (original !== replacement || originalAligned !== record.sequence) changedRows.push(record.name);
		const edits = nucleotideEditCounts(original, replacement);
		const rowSubstitutions = edits.substitutedNucleotides, rowRemoved = edits.removedNucleotides, rowInserted = edits.insertedNucleotides;
		substitutedNucleotides += rowSubstitutions;
		removedNucleotides += rowRemoved;
		insertedNucleotides += rowInserted;
		if (original !== replacement || originalAligned !== record.sequence) rowChanges.push({
			name: record.name,
			substitutedNucleotides: rowSubstitutions,
			insertedNucleotides: rowInserted,
			removedNucleotides: rowRemoved,
			gapPlacementChanged: originalAligned !== record.sequence && original === replacement
		});
	}
	return {
		...corrected,
		removedRows,
		addedRows,
		changedRows,
		rowChanges,
		removedNucleotides,
		insertedNucleotides,
		substitutedNucleotides
	};
}
function summarizeAlignmentChanges(currentText, correctedText) {
	const current = inspectAlignment(currentText, 1), corrected = validateCorrectedAlignment(currentText, correctedText);
	const currentNames = new Set(current.records.map((record) => record.name)), correctedNames = new Set(corrected.records.map((record) => record.name));
	const sharedBefore = current.records.map((record) => record.name).filter((name) => correctedNames.has(name));
	const sharedAfter = corrected.records.map((record) => record.name).filter((name) => currentNames.has(name));
	return {
		rowsBefore: current.rows,
		rowsAfter: corrected.rows,
		columnsBefore: current.columns,
		columnsAfter: corrected.columns,
		rowOrderChanged: sharedBefore.join("\0") !== sharedAfter.join("\0"),
		rowOrderBefore: current.records.map((record) => record.name),
		rowOrderAfter: corrected.records.map((record) => record.name),
		removedRows: corrected.removedRows,
		addedRows: corrected.addedRows,
		changedRows: corrected.changedRows,
		rowChanges: corrected.rowChanges,
		removedNucleotides: corrected.removedNucleotides,
		insertedNucleotides: corrected.insertedNucleotides,
		substitutedNucleotides: corrected.substitutedNucleotides
	};
}
function translateAlignedNucleotides(sequence, frameOffset = 0) {
	const offset = frameOffset === 1 || frameOffset === 2 ? frameOffset : 0;
	let result = "";
	for (let index = offset; index < sequence.length; index += 3) {
		const codon = sequence.slice(index, index + 3).toUpperCase();
		result += codon === "---" ? "-" : codon.includes("-") || /[^ACGT]/.test(codon) || codon.length < 3 ? "X" : CODONS$1[codon] ?? "X";
	}
	return result;
}
function translateAlignmentFasta(fasta, frameOffset = 0) {
	return exactFasta$1(inspectAlignment(fasta, 1).records.map((record) => ({
		...record,
		sequence: translateAlignedNucleotides(record.sequence, frameOffset)
	})));
}
//#endregion
//#region src/collapse.ts
function exactFasta(rows) {
	return rows.map((row) => `>${row.name}\n${row.sequence}\n`).join("");
}
/**
* Collapse identical aligned haplotypes while retaining one count per UMI
* family. Read-family sizes are metadata only and never contribute to the
* collapse multiplicity. Per-family agreement is deliberately not projected
* onto the resulting haplotype.
*/
function collapseAlignment(fasta, sample) {
	const alignment = inspectAlignment(fasta, 1);
	const bySequence = /* @__PURE__ */ new Map();
	for (const row of alignment.records) {
		const haplotype = row.sequence.replaceAll("-", "");
		let group = bySequence.get(haplotype);
		if (!group) {
			group = {
				name: row.name,
				sequence: row.sequence,
				members: []
			};
			bySequence.set(haplotype, group);
		}
		group.members.push(row.name);
	}
	const groups = [...bySequence.values()].map((group) => ({
		sample,
		representativeId: group.name,
		memberIds: group.members,
		familyCount: group.members.length
	}));
	return {
		fasta: exactFasta([...bySequence.values()].map((group) => ({
			name: group.name,
			sequence: group.sequence
		}))),
		groups
	};
}
//#endregion
//#region src/panel-profile.ts
function profile(rows) {
	if (!rows.length) return [];
	const width = rows[0].length;
	if (!width || rows.some((row) => row.length !== width)) throw new Error("Panel and sample profiles require rectangular alignments.");
	return Array.from({ length: width }, (_, column) => {
		const counts = /* @__PURE__ */ new Map();
		for (const row of rows) {
			const value = row[column].toUpperCase();
			counts.set(value, (counts.get(value) ?? 0) + 1);
		}
		return new Map([...counts].map(([value, count]) => [value, count / rows.length]));
	});
}
function cost(left, right) {
	let equal = 0;
	for (const [base, probability] of left) equal += probability * (right.get(base) ?? 0);
	return 2 * equal - 1;
}
function alignProfiles(panel, sample) {
	const rows = panel.length + 1, columns = sample.length + 1, cells = rows * columns;
	if (!Number.isSafeInteger(cells) || cells > 512 * 1024 * 1024) throw new Error("The panel/profile alignment is too large for browser memory.");
	const trace = new Uint8Array(cells), gapOpen = -2, gapExtend = -.2, negative = Number.NEGATIVE_INFINITY;
	let previousM = new Float64Array(columns), previousX = new Float64Array(columns), previousY = new Float64Array(columns);
	let currentM = new Float64Array(columns), currentX = new Float64Array(columns), currentY = new Float64Array(columns);
	for (let column = 0; column < columns; column++) {
		previousM[column] = gapOpen + gapExtend * column;
		previousX[column] = negative;
		previousY[column] = gapOpen + gapExtend * column;
	}
	for (let row = 1; row < rows; row++) {
		currentM[0] = gapOpen + gapExtend * row;
		currentX[0] = gapOpen + gapExtend * row;
		currentY[0] = negative;
		for (let column = 1; column < columns; column++) {
			const diagonalCost = cost(panel[row - 1], sample[column - 1]);
			let bestM = previousM[column - 1] + diagonalCost, fromM = 0;
			const fromX = previousX[column - 1] + diagonalCost, fromY = previousY[column - 1] + diagonalCost;
			if (fromX > bestM) {
				bestM = fromX;
				fromM = 1;
			}
			if (fromY > bestM) {
				bestM = fromY;
				fromM = 2;
			}
			currentM[column] = bestM;
			const openX = previousM[column] + gapOpen, extendX = previousX[column] + gapExtend;
			const xContinues = extendX > openX;
			currentX[column] = xContinues ? extendX : openX;
			const openY = currentM[column - 1] + gapOpen, extendY = currentY[column - 1] + gapExtend;
			const yContinues = extendY > openY;
			currentY[column] = yContinues ? extendY : openY;
			trace[row * columns + column] = fromM | (xContinues ? 4 : 0) | (yContinues ? 8 : 0);
		}
		[previousM, currentM] = [currentM, previousM];
		[previousX, currentX] = [currentX, previousX];
		[previousY, currentY] = [currentY, previousY];
	}
	let state = 0, maximum = previousM[columns - 1];
	if (previousX[columns - 1] > maximum) {
		maximum = previousX[columns - 1];
		state = 1;
	}
	if (previousY[columns - 1] > maximum) state = 2;
	let row = panel.length, column = sample.length;
	const reverse = [];
	while (row > 0 && column > 0) {
		const packed = trace[row * columns + column];
		if (state === 0) {
			reverse.push(row - 1);
			row--;
			column--;
			state = packed & 3;
		} else if (state === 1) {
			row--;
			state = packed & 4 ? 1 : 0;
		} else {
			reverse.push(-1);
			column--;
			state = packed & 8 ? 2 : 0;
		}
	}
	while (column-- > 0) reverse.push(-1);
	reverse.reverse();
	if (reverse.length !== sample.length) throw new Error("The panel/profile alignment produced an invalid sample map.");
	return reverse;
}
function maximumSubarray(values) {
	let current = 0, best = 1;
	for (const value of values) {
		current += value;
		best = Math.max(best, current);
		if (current <= 0) current = 0;
	}
	return best;
}
function extractAndScorePanel(sampleRows, panelRows) {
	if (!sampleRows.length) return {
		sequences: [],
		scores: []
	};
	const sampleProfile = profile(sampleRows), panelProfile = profile(panelRows), sampleToPanel = alignProfiles(panelProfile, sampleProfile);
	const start = sampleToPanel.findIndex((index) => index >= 0);
	let end = sampleToPanel.length - 1;
	while (end >= 0 && sampleToPanel[end] < 0) end--;
	if (start < 0 || end < 0) throw new Error("The sample alignment does not overlap the reference panel.");
	const sequences = sampleRows.map((row) => row.slice(start, end + 1));
	return {
		sequences,
		scores: sequences.map((sequence) => {
			const transformed = [];
			for (let offset = 0; offset < sequence.length; offset++) {
				const base = sequence[offset].toUpperCase(), panelIndex = sampleToPanel[start + offset];
				if (base === "-" || panelIndex < 0) continue;
				const probability = panelProfile[panelIndex].get(base) ?? 0;
				transformed.push(-Math.log(probability + .01) + Math.log(.25));
			}
			return maximumSubarray(transformed);
		})
	};
}
//#endregion
//#region src/scalable-msa.ts
const MONOLITHIC_ROWS = 8e3;
const MONOLITHIC_BASES = 128 * 1024 * 1024;
const PROFILE_BATCH_ROWS = 2e3;
function decompose(anchor, alignedAnchor, alignedRow) {
	if (alignedAnchor.length !== alignedRow.length || alignedAnchor.replaceAll("-", "") !== anchor) throw new Error("A scalable MSA batch returned an invalid anchor alignment.");
	const insertions = Array.from({ length: anchor.length + 1 }, () => ""), bases = Array(anchor.length);
	let position = 0;
	for (let column = 0; column < alignedAnchor.length; column++) if (alignedAnchor[column] === "-") insertions[position] += alignedRow[column];
	else {
		if (alignedAnchor[column] !== anchor[position]) throw new Error("A scalable MSA batch changed the anchor sequence.");
		bases[position++] = alignedRow[column];
	}
	if (position !== anchor.length) throw new Error("A scalable MSA batch truncated the anchor sequence.");
	return {
		insertions,
		bases
	};
}
function padInsertion(value, width, slot) {
	return slot === 0 ? value.padStart(width, "-") : value.padEnd(width, "-");
}
async function runScalableMsa(sequences, runMsa, signal, iterations = 3, scoringMode = "literal") {
	if (sequences.length < 2) return [...sequences];
	const totalBases = sequences.reduce((sum, sequence) => sum + sequence.length, 0);
	if (sequences.length <= MONOLITHIC_ROWS && totalBases <= MONOLITHIC_BASES) return runMsa(sequences, signal, iterations, scoringMode);
	const anchor = sequences[0], widths = new Uint32Array(anchor.length + 1);
	const decomposed = Array(sequences.length);
	for (let start = 0; start < sequences.length; start += PROFILE_BATCH_ROWS) {
		if (signal?.aborted) throw new DOMException("Analysis cancelled.", "AbortError");
		const end = Math.min(sequences.length, start + PROFILE_BATCH_ROWS);
		const input = start === 0 ? sequences.slice(start, end) : [anchor, ...sequences.slice(start, end)];
		const aligned = await runMsa(input, signal, iterations, scoringMode);
		if (aligned.length !== input.length || aligned.some((row) => row.length !== aligned[0].length)) throw new Error("A scalable MSA batch returned a non-rectangular alignment.");
		const alignedAnchor = aligned[0], first = start === 0 ? 0 : 1;
		for (let row = first; row < aligned.length; row++) {
			const index = start === 0 ? row : start + row - 1;
			if (aligned[row].replaceAll("-", "") !== sequences[index]) throw new Error("A scalable MSA batch changed an input sequence.");
			const parts = decompose(anchor, alignedAnchor, aligned[row]);
			decomposed[index] = parts;
			parts.insertions.forEach((value, slot) => {
				widths[slot] = Math.max(widths[slot], value.length);
			});
		}
	}
	return decomposed.map((row) => {
		if (!row) throw new Error("A scalable MSA batch omitted a sequence.");
		let output = padInsertion(row.insertions[0], widths[0], 0);
		for (let position = 0; position < anchor.length; position++) output += row.bases[position] + padInsertion(row.insertions[position + 1], widths[position + 1], position + 1);
		return output;
	});
}
//#endregion
//#region src/postprocess.ts
const degap = (sequence) => sequence.replaceAll("-", "").toUpperCase();
const fasta$1 = (rows) => rows.map((row) => `>${row.name}\n${row.sequence.match(/.{1,80}/g)?.join("\n") ?? ""}`).join("\n") + (rows.length ? "\n" : "");
function quantile(values, probability) {
	if (!values.length) return 0;
	const sorted = [...values].sort((a, b) => a - b);
	if (sorted.length === 1) return sorted[0];
	const position = Math.max(0, Math.min(1, probability)) * (sorted.length - 1), lower = Math.floor(position), fraction = position - lower;
	return sorted[lower] * (1 - fraction) + sorted[Math.min(lower + 1, sorted.length - 1)] * fraction;
}
function alignmentConsensus(rows) {
	if (!rows.length) return "";
	let output = "";
	for (let position = 0; position < rows[0].length; position++) {
		const counts = /* @__PURE__ */ new Map();
		rows.forEach((row, index) => {
			const base = row[position].toUpperCase(), value = counts.get(base) ?? {
				count: 0,
				first: index
			};
			value.count++;
			counts.set(base, value);
		});
		output += [...counts].sort((a, b) => b[1].count - a[1].count || a[1].first - b[1].first)[0][0];
	}
	return output;
}
function matrixMultiply(left, right) {
	const output = Array(16).fill(0);
	for (let i = 0; i < 4; i++) for (let k = 0; k < 4; k++) for (let j = 0; j < 4; j++) output[i * 4 + j] += left[i * 4 + k] * right[k * 4 + j];
	return output;
}
const identity = () => [
	1,
	0,
	0,
	0,
	0,
	1,
	0,
	0,
	0,
	0,
	1,
	0,
	0,
	0,
	0,
	1
];
function matrixExp(matrix) {
	const norm = Math.max(...[
		0,
		1,
		2,
		3
	].map((row) => matrix.slice(row * 4, row * 4 + 4).reduce((sum, value) => sum + Math.abs(value), 0)));
	const scale = Math.max(0, Math.ceil(Math.log2(Math.max(norm, 1)))), divisor = 2 ** scale;
	matrix = matrix.map((value) => value / divisor);
	let result = identity(), term = identity();
	for (let order = 1; order <= 24; order++) {
		term = matrixMultiply(term, matrix).map((value) => value / order);
		result = result.map((value, index) => value + term[index]);
	}
	for (let index = 0; index < scale; index++) result = matrixMultiply(result, result);
	return result;
}
function normalPdf(value, mean, sd) {
	return Math.exp(-((value - mean) ** 2) / (2 * sd * sd)) / (sd * Math.sqrt(2 * Math.PI));
}
let gridCache;
function apobecGrid() {
	if (gridCache) return gridCache;
	const output = [];
	for (let ti = 0; ti < 23; ti++) for (let gi = 0; gi < 121; gi++) {
		const t = -12 + ti * .5, ga = -1 + gi * .05, mu = Math.exp(t), multiplier = Math.exp(ga), tv = 4.5;
		const q = [
			-(1 + tv + 1),
			1,
			tv,
			1,
			1,
			-(2 + tv),
			1,
			tv,
			tv * multiplier,
			1,
			-(tv * multiplier + 1 + 1),
			1,
			1,
			tv,
			1,
			-(1 + tv + 1)
		].map((value) => value * mu);
		const prior = Math.log(normalPdf(t, -5, 1)) + Math.log(.99 * normalPdf(ga, 0, .1) + .01 * normalPdf(ga, 0, 1));
		output.push({
			t,
			ga,
			transitions: matrixExp(q),
			prior
		});
	}
	return gridCache = output;
}
const baseIndex = {
	A: 0,
	C: 1,
	G: 2,
	T: 3
};
function apobec(consensus, sequence) {
	const counts = Array(16).fill(0);
	for (let position = 0; position < Math.min(consensus.length, sequence.length); position++) {
		const left = baseIndex[consensus[position].toUpperCase()], right = baseIndex[sequence[position].toUpperCase()];
		if (left !== void 0 && right !== void 0) counts[left * 4 + right]++;
	}
	const weights = apobecGrid().map((point) => {
		let value = point.prior;
		for (let index = 0; index < 16; index++) if (counts[index]) value += counts[index] * Math.log(Math.max(point.transitions[index], 1e-300));
		return value;
	});
	const maximum = Math.max(...weights), scaled = weights.map((value) => Math.exp(value - maximum)), total = scaled.reduce((a, b) => a + b, 0);
	let t = 0, ga = 0, inflated = 0;
	apobecGrid().forEach((point, index) => {
		const probability = scaled[index] / total;
		t += point.t * probability;
		ga += point.ga * probability;
		if (point.ga > 0) inflated += probability;
	});
	let totalMutations = 0;
	for (let a = 0; a < 4; a++) for (let b = 0; b < 4; b++) if (a !== b) totalMutations += counts[a * 4 + b];
	return {
		posteriorMeanGaMultiplier: Math.exp(ga),
		posteriorGaInflated: inflated,
		posteriorMeanMutationRate: Math.exp(t),
		gaMutations: counts[8],
		totalMutations
	};
}
const CODONS = {
	TTT: "F",
	TTC: "F",
	TTA: "L",
	TTG: "L",
	CTT: "L",
	CTC: "L",
	CTA: "L",
	CTG: "L",
	ATT: "I",
	ATC: "I",
	ATA: "I",
	ATG: "M",
	GTT: "V",
	GTC: "V",
	GTA: "V",
	GTG: "V",
	TCT: "S",
	TCC: "S",
	TCA: "S",
	TCG: "S",
	AGT: "S",
	AGC: "S",
	CCT: "P",
	CCC: "P",
	CCA: "P",
	CCG: "P",
	ACT: "T",
	ACC: "T",
	ACA: "T",
	ACG: "T",
	GCT: "A",
	GCC: "A",
	GCA: "A",
	GCG: "A",
	TAT: "Y",
	TAC: "Y",
	TAA: "*",
	TAG: "*",
	TGA: "*",
	CAT: "H",
	CAC: "H",
	CAA: "Q",
	CAG: "Q",
	AAT: "N",
	AAC: "N",
	AAA: "K",
	AAG: "K",
	GAT: "D",
	GAC: "D",
	GAA: "E",
	GAG: "E",
	TGT: "C",
	TGC: "C",
	TGG: "W",
	CGT: "R",
	CGC: "R",
	CGA: "R",
	CGG: "R",
	AGA: "R",
	AGG: "R",
	GGT: "G",
	GGC: "G",
	GGA: "G",
	GGG: "G"
};
function translate(sequence) {
	let output = "";
	for (let index = 0; index + 2 < sequence.length; index += 3) output += CODONS[sequence.slice(index, index + 3)] ?? "X";
	return output;
}
function longestOrf(sequence) {
	let best;
	for (let frame = 0; frame < 3; frame++) {
		const aa = translate(sequence.slice(frame, frame + Math.floor((sequence.length - frame) / 3) * 3));
		for (let start = 0; start < aa.length; start++) if (aa[start] === "M") {
			const offset = aa.indexOf("*", start);
			if (offset >= 0 && (!best || offset + 1 - start > best.length)) best = {
				length: offset + 1 - start,
				start: frame + start * 3,
				end: frame + (offset + 1) * 3
			};
		}
	}
	return best ? sequence.slice(best.start, best.end) : sequence.slice(0, Math.floor(sequence.length / 3) * 3);
}
function backtranslate(alignedAminoAcids, codingNucleotides) {
	let output = "", offset = 0;
	for (const residue of alignedAminoAcids) if (residue === "-") output += "---";
	else {
		output += codingNucleotides.slice(offset, offset + 3);
		offset += 3;
	}
	if (offset !== codingNucleotides.length) throw new Error("The amino-acid alignment did not preserve its nucleotide coding sequence.");
	return output;
}
/**
* SeededAlignment's functional filter is codon-aware and evaluates every
* sequence against one joint reference-anchored alignment.  Doing independent
* nucleotide NW alignments introduced arbitrary one-base gaps and incorrectly
* labeled almost every long-read sequence as a frameshift.  Translate first,
* align the complete batch, and project gaps back as codon triplets instead.
*/
async function functionalFilterBatch(reference, sequences, threshold, runMsa, signal) {
	const outcomes = Array(sequences.length), coding = [];
	sequences.forEach((raw, index) => {
		const sequence = degap(raw);
		if (/[^ACGT]/.test(sequence)) outcomes[index] = {
			passed: false,
			reasons: ["ambiguousSymbols-reject"]
		};
		else coding.push({
			index,
			sequence: longestOrf(sequence)
		});
	});
	const referenceCoding = degap(reference).slice(0, Math.floor(degap(reference).length / 3) * 3);
	if (!coding.length) return {
		outcomes,
		referenceNt: referenceCoding,
		referenceAa: translate(referenceCoding)
	};
	const nucleotideAlignment = (await runScalableMsa([translate(referenceCoding), ...coding.map((row) => translate(row.sequence))], runMsa, signal, 3, "amino-acid")).map((row, index) => backtranslate(row, index ? coding[index - 1].sequence : referenceCoding));
	const alignedReference = nucleotideAlignment[0], first = alignedReference.search(/[^-]/);
	let last = alignedReference.length - 1;
	while (last >= first && alignedReference[last] === "-") last--;
	const referenceRegion = alignedReference.slice(first, last + 1);
	coding.forEach((entry, position) => {
		const queryRegion = nucleotideAlignment[position + 1].slice(first, last + 1), reasons = [];
		const trimmed = degap(queryRegion), aa = translate(trimmed);
		if (queryRegion.slice(0, 3) !== "ATG") reasons.push("lateStart-reject");
		if (queryRegion.slice(-3) === "---") reasons.push("earlyStop-reject");
		let matches = 0;
		for (let column = 0; column < referenceRegion.length; column += 1) if (referenceRegion[column] === queryRegion[column]) matches++;
		const rawRatio = matches / Math.max(1, referenceRegion.length), digits = rawRatio === 0 ? 3 : 3 - Math.floor(Math.log10(Math.abs(rawRatio))) - 1;
		const ratio = Number(rawRatio.toFixed(Math.max(0, digits)));
		if (ratio < threshold) reasons.push(`badMatch-reject (match=${ratio})`);
		outcomes[entry.index] = {
			passed: !reasons.length,
			reasons,
			nt: trimmed,
			aa,
			alignedNt: queryRegion,
			alignedAa: translateAlignedNucleotides(queryRegion, 0)
		};
	});
	return {
		outcomes,
		referenceNt: referenceRegion,
		referenceAa: translateAlignedNucleotides(referenceRegion, 0)
	};
}
async function postprocess(consensuses, contamination, config, signal, runMsa = runAlivibeMsa, sampleConcurrency = 1, onProgress, options = {}) {
	const discarded = new Set(contamination.filter((call) => call.discarded).map((call) => call.sequenceId));
	const outputs = Array(config.samples.length);
	let cursor = 0;
	const sampleProgress = Array(config.samples.length).fill(0);
	const report = (sampleIndex, fraction, detail) => {
		sampleProgress[sampleIndex] = Math.max(sampleProgress[sampleIndex], Math.max(0, Math.min(1, fraction)));
		onProgress?.({
			fraction: sampleProgress.reduce((sum, value) => sum + value, 0) / Math.max(1, sampleProgress.length),
			detail
		});
	};
	await Promise.all(Array.from({ length: Math.min(config.samples.length, Math.max(1, Math.floor(sampleConcurrency))) }, async () => {
		while (true) {
			const sampleIndex = cursor++;
			if (sampleIndex >= config.samples.length) return;
			const sample = config.samples[sampleIndex], records = [], summaries = [], alignments = {};
			const referenceAlignments = {}, collapseGroups = {};
			if (signal?.aborted) throw new DOMException("Analysis cancelled.", "AbortError");
			report(sampleIndex, .03, `Preparing filters for sample ${sample.name}`);
			const source = consensuses.filter((record) => record.sample === sample.name), sizes = source.filter((record) => !discarded.has(record.id)).map((record) => record.familySize);
			const artefactCutoff = Math.ceil(quantile(sizes, sample.outlierQuantileOverride ?? config.parameters.outlierQuantile) * (sample.artefactFractionOverride ?? config.parameters.artefactFraction));
			const agreementThreshold = sample.agreementOverride ?? config.parameters.agreementThreshold;
			const preliminary = source.map((record, index) => ({
				record,
				index
			})).filter(({ record }) => record.familySize >= artefactCutoff && record.minimumAgreement >= agreementThreshold && !discarded.has(record.id));
			const scores = Array(source.length).fill(0), panelPass = Array(source.length).fill(true), extracted = /* @__PURE__ */ new Map();
			if (preliminary.length && sample.panelSequences.length) {
				report(sampleIndex, .18, `Screening candidate sequences against the reference panel for sample ${sample.name}`);
				const panelResult = extractAndScorePanel(preliminary.length > 1 ? await runScalableMsa(preliminary.map(({ record }) => degap(record.sequence)), runMsa, signal, 3, "nucleotide") : preliminary.map(({ record }) => degap(record.sequence)), sample.panelSequences.map((record) => record.sequence));
				preliminary.forEach(({ index }, candidate) => {
					extracted.set(index, degap(panelResult.sequences[candidate]));
					scores[index] = panelResult.scores[candidate];
					panelPass[index] = scores[index] < config.parameters.panelThreshold;
				});
			} else preliminary.forEach(({ record, index }) => extracted.set(index, degap(record.sequence)));
			report(sampleIndex, .4, `Reference-panel screening complete for sample ${sample.name}`);
			const accepted = preliminary.filter(({ index }) => panelPass[index]);
			const displayReference = sample.functionalReferenceSequence?.sequence ?? sample.panelSequences[0]?.sequence;
			let acceptedAlignment = [], alignedReference = "";
			if (accepted.length) {
				report(sampleIndex, .48, `Building the retained-sequence alignment for sample ${sample.name}`);
				const inputs = [...displayReference ? [degap(displayReference)] : [], ...accepted.map(({ index }) => extracted.get(index))];
				const aligned = inputs.length > 1 ? await runScalableMsa(inputs, runMsa, signal, 3, "nucleotide") : inputs;
				if (displayReference) {
					alignedReference = aligned[0];
					acceptedAlignment = aligned.slice(1);
				} else {
					acceptedAlignment = aligned;
					alignedReference = alignmentConsensus(aligned);
				}
			}
			report(sampleIndex, .66, `Retained-sequence alignment complete for sample ${sample.name}`);
			const alignmentByIndex = new Map(accepted.map(({ index }, position) => [index, acceptedAlignment[position]]));
			const consensus = alignmentConsensus(acceptedAlignment), nucleotideRows = [];
			const functionalNucleotideRows = [];
			const functionalProteinRows = [];
			const functionalByIndex = /* @__PURE__ */ new Map();
			let functionalReferenceNt = "";
			if (sample.functionalReferenceSequence && accepted.length) {
				report(sampleIndex, .72, `Checking coding-frame and functional-reference requirements for sample ${sample.name}`);
				const batch = await functionalFilterBatch(sample.functionalReferenceSequence.sequence, accepted.map(({ index }) => extracted.get(index)), sample.functionalMatchOverride ?? config.parameters.functionalMatchThreshold, runMsa, signal);
				functionalReferenceNt = batch.referenceNt;
				batch.referenceAa;
				accepted.forEach(({ index }, position) => functionalByIndex.set(index, batch.outcomes[position]));
			}
			report(sampleIndex, .84, `Calculating sequence annotations and filter decisions for sample ${sample.name}`);
			let functionalPassed = 0;
			for (const [index, record] of source.entries()) {
				if (signal?.aborted) throw new DOMException("Downstream filtering skipped.", "AbortError");
				const artefactPass = record.familySize >= artefactCutoff, agreementPass = record.minimumAgreement >= agreementThreshold;
				const contaminationPass = !discarded.has(record.id), acceptedRow = alignmentByIndex.get(index), rejectionReasons = [];
				if (!artefactPass) rejectionReasons.push(`ccs_count < artefact cutoff (${artefactCutoff})`);
				if (!agreementPass) rejectionReasons.push(`minimum_agreement < ${agreementThreshold}`);
				if (!contaminationPass) rejectionReasons.push("contamination filter");
				if (!panelPass[index]) rejectionReasons.push(`distance_from_panel >= ${config.parameters.panelThreshold}`);
				let trimmedNt, trimmedAa, functionalPass;
				if (sample.functionalReferenceSequence) {
					if (acceptedRow) {
						const outcome = functionalByIndex.get(index);
						functionalPass = outcome.passed;
						trimmedNt = outcome.nt;
						trimmedAa = outcome.aa;
						rejectionReasons.push(...outcome.reasons);
						if (outcome.passed) {
							functionalPassed++;
							if (outcome.alignedNt && outcome.alignedAa) {
								functionalNucleotideRows.push({
									name: record.id,
									sequence: outcome.alignedNt
								});
								functionalProteinRows.push({
									name: record.id,
									sequence: outcome.alignedAa
								});
							}
						}
					}
				}
				if (acceptedRow) nucleotideRows.push({
					name: record.id,
					sequence: acceptedRow
				});
				records.push({
					id: record.id,
					sample: sample.name,
					umi: record.umi,
					familySize: record.familySize,
					minimumAgreement: record.minimumAgreement,
					consensusNt: record.sequence,
					alignedNt: acceptedRow,
					trimmedNt,
					trimmedAa,
					panelScore: scores[index],
					artefactPass,
					agreementPass,
					contaminationPass,
					panelPass: panelPass[index],
					functionalPass,
					rejectionReasons,
					apobec: acceptedRow ? apobec(consensus, acceptedRow) : void 0
				});
				if ((index & 7) === 7 || index + 1 === source.length) {
					report(sampleIndex, .84 + .12 * (index + 1) / Math.max(1, source.length), `Annotated ${index + 1} of ${source.length} consensus-family records for ${sample.name}`);
					await new Promise((resolvePromise) => setTimeout(resolvePromise, 0));
				}
			}
			if (nucleotideRows.length) {
				const uncollapsed = fasta$1(nucleotideRows);
				alignments[`${sample.name}/uncollapsed-nucleotide`] = uncollapsed;
				alignments[`${sample.name}/uncollapsed-protein`] = fasta$1(nucleotideRows.map((row) => ({
					...row,
					sequence: translateAlignedNucleotides(row.sequence, 0)
				})));
				referenceAlignments[`${sample.name}/uncollapsed-nucleotide`] = fasta$1([{
					name: "reference",
					sequence: alignedReference
				}]);
			}
			if (functionalNucleotideRows.length) {
				alignments[`${sample.name}/functional-nucleotide`] = fasta$1(functionalNucleotideRows);
				alignments[`${sample.name}/functional-protein`] = fasta$1(functionalProteinRows);
				referenceAlignments[`${sample.name}/functional-nucleotide`] = fasta$1([{
					name: "functional_reference",
					sequence: functionalReferenceNt
				}]);
			}
			summaries.push({
				sample: sample.name,
				demultiplexedReads: 0,
				observedUmis: 0,
				likelyRealUmis: 0,
				consensusSequences: source.length,
				contaminationPassed: source.filter((record) => !discarded.has(record.id)).length,
				postprocPassed: accepted.length,
				functionalPassed: sample.functionalReferenceSequence ? functionalPassed : void 0,
				artefactCutoff
			});
			outputs[sampleIndex] = {
				records,
				summaries,
				alignments,
				referenceAlignments,
				collapseGroups,
				collapseSeconds: 0
			};
			report(sampleIndex, 1, `Downstream processing complete for sample ${sample.name}`);
		}
	}));
	const combined = {
		records: outputs.flatMap((output) => output.records),
		summaries: outputs.flatMap((output) => output.summaries),
		alignments: Object.assign({}, ...outputs.map((output) => output.alignments)),
		referenceAlignments: Object.assign({}, ...outputs.map((output) => output.referenceAlignments)),
		collapseGroups: Object.assign({}, ...outputs.map((output) => output.collapseGroups)),
		collapseSeconds: 0
	};
	return options.collapse === false ? combined : collapsePostprocess(combined, config, signal, options.onCollapseProgress);
}
/** Run family-count-preserving haplotype collapse as its own resumable stage. */
async function collapsePostprocess(output, config, signal, onProgress) {
	const started = performance.now(), alignments = { ...output.alignments }, referenceAlignments = { ...output.referenceAlignments };
	const collapseGroups = { ...output.collapseGroups }, summaries = output.summaries.map((summary) => ({ ...summary }));
	for (const [index, sample] of config.samples.entries()) {
		if (signal?.aborted) throw new DOMException("Haplotype collapse skipped.", "AbortError");
		onProgress?.({
			fraction: index / Math.max(1, config.samples.length),
			detail: `Collapsing identical retained UMI-family sequences for ${sample.name}`
		});
		await new Promise((resolvePromise) => setTimeout(resolvePromise, 0));
		const uncollapsed = alignments[`${sample.name}/uncollapsed-nucleotide`];
		let collapsedCount = 0;
		if (uncollapsed) {
			const collapsed = collapseAlignment(uncollapsed, sample.name);
			collapsedCount = collapsed.groups.length;
			alignments[`${sample.name}/nucleotide`] = collapsed.fasta;
			collapseGroups[sample.name] = collapsed.groups;
			const rows = new Map(uncollapsed.split(/^>/m).filter(Boolean).map((block) => {
				const [name, ...sequence] = block.trimEnd().split(/\r?\n/);
				return [name, sequence.join("")];
			}));
			alignments[`${sample.name}/protein`] = fasta$1(collapsed.groups.map((group) => ({
				name: group.representativeId,
				sequence: translateAlignedNucleotides(rows.get(group.representativeId) ?? "", 0)
			})));
			const reference = referenceAlignments[`${sample.name}/uncollapsed-nucleotide`];
			if (reference) referenceAlignments[`${sample.name}/nucleotide`] = reference;
		}
		const summary = summaries.find((row) => row.sample === sample.name);
		if (summary) summary.collapsedSequences = collapsedCount;
		onProgress?.({
			fraction: (index + 1) / Math.max(1, config.samples.length),
			detail: `Collapsed ${sample.name} into ${collapsedCount.toLocaleString()} haplotypes; counts represent UMI families`
		});
	}
	return {
		...output,
		summaries,
		alignments,
		referenceAlignments,
		collapseGroups,
		collapseSeconds: output.collapseSeconds + (performance.now() - started) / 1e3
	};
}
//#endregion
//#region src/optional-stages.ts
function downstreamResources(config) {
	return { samples: config.samples.map((sample) => ({
		name: sample.name,
		panelSequences: sample.panelSequences.map((record) => ({ ...record })),
		functionalReferenceSequence: sample.functionalReferenceSequence ? { ...sample.functionalReferenceSequence } : void 0
	})) };
}
function statusRecord(state, detail, updatedUtc = (/* @__PURE__ */ new Date()).toISOString()) {
	return {
		state,
		detail,
		updatedUtc
	};
}
//#endregion
//#region ../../work/webPORPID/node_modules/@msgpack/msgpack/dist.esm/utils/utf8.mjs
function utf8Count(str) {
	const strLength = str.length;
	let byteLength = 0;
	let pos = 0;
	while (pos < strLength) {
		let value = str.charCodeAt(pos++);
		if ((value & 4294967168) === 0) {
			byteLength++;
			continue;
		} else if ((value & 4294965248) === 0) byteLength += 2;
		else {
			if (value >= 55296 && value <= 56319) {
				if (pos < strLength) {
					const extra = str.charCodeAt(pos);
					if ((extra & 64512) === 56320) {
						++pos;
						value = ((value & 1023) << 10) + (extra & 1023) + 65536;
					}
				}
			}
			if ((value & 4294901760) === 0) byteLength += 3;
			else byteLength += 4;
		}
	}
	return byteLength;
}
function utf8EncodeJs(str, output, outputOffset) {
	const strLength = str.length;
	let offset = outputOffset;
	let pos = 0;
	while (pos < strLength) {
		let value = str.charCodeAt(pos++);
		if ((value & 4294967168) === 0) {
			output[offset++] = value;
			continue;
		} else if ((value & 4294965248) === 0) output[offset++] = value >> 6 & 31 | 192;
		else {
			if (value >= 55296 && value <= 56319) {
				if (pos < strLength) {
					const extra = str.charCodeAt(pos);
					if ((extra & 64512) === 56320) {
						++pos;
						value = ((value & 1023) << 10) + (extra & 1023) + 65536;
					}
				}
			}
			if ((value & 4294901760) === 0) {
				output[offset++] = value >> 12 & 15 | 224;
				output[offset++] = value >> 6 & 63 | 128;
			} else {
				output[offset++] = value >> 18 & 7 | 240;
				output[offset++] = value >> 12 & 63 | 128;
				output[offset++] = value >> 6 & 63 | 128;
			}
		}
		output[offset++] = value & 63 | 128;
	}
}
const sharedTextEncoder = new TextEncoder();
const TEXT_ENCODER_THRESHOLD = 50;
function utf8EncodeTE(str, output, outputOffset) {
	sharedTextEncoder.encodeInto(str, output.subarray(outputOffset));
}
function utf8Encode(str, output, outputOffset) {
	if (str.length > TEXT_ENCODER_THRESHOLD) utf8EncodeTE(str, output, outputOffset);
	else utf8EncodeJs(str, output, outputOffset);
}
const CHUNK_SIZE = 4096;
function utf8DecodeJs(bytes, inputOffset, byteLength) {
	let offset = inputOffset;
	const end = offset + byteLength;
	const units = [];
	let result = "";
	while (offset < end) {
		const byte1 = bytes[offset++];
		if ((byte1 & 128) === 0) units.push(byte1);
		else if ((byte1 & 224) === 192) {
			const byte2 = bytes[offset++] & 63;
			units.push((byte1 & 31) << 6 | byte2);
		} else if ((byte1 & 240) === 224) {
			const byte2 = bytes[offset++] & 63;
			const byte3 = bytes[offset++] & 63;
			units.push((byte1 & 31) << 12 | byte2 << 6 | byte3);
		} else if ((byte1 & 248) === 240) {
			const byte2 = bytes[offset++] & 63;
			const byte3 = bytes[offset++] & 63;
			const byte4 = bytes[offset++] & 63;
			let unit = (byte1 & 7) << 18 | byte2 << 12 | byte3 << 6 | byte4;
			if (unit > 65535) {
				unit -= 65536;
				units.push(unit >>> 10 & 1023 | 55296);
				unit = 56320 | unit & 1023;
			}
			units.push(unit);
		} else units.push(byte1);
		if (units.length >= CHUNK_SIZE) {
			result += String.fromCharCode(...units);
			units.length = 0;
		}
	}
	if (units.length > 0) result += String.fromCharCode(...units);
	return result;
}
const sharedTextDecoder = new TextDecoder();
const TEXT_DECODER_THRESHOLD = 200;
function utf8DecodeTD(bytes, inputOffset, byteLength) {
	const stringBytes = bytes.subarray(inputOffset, inputOffset + byteLength);
	return sharedTextDecoder.decode(stringBytes);
}
function utf8Decode(bytes, inputOffset, byteLength) {
	if (byteLength > TEXT_DECODER_THRESHOLD) return utf8DecodeTD(bytes, inputOffset, byteLength);
	else return utf8DecodeJs(bytes, inputOffset, byteLength);
}
//#endregion
//#region ../../work/webPORPID/node_modules/@msgpack/msgpack/dist.esm/ExtData.mjs
/**
* ExtData is used to handle Extension Types that are not registered to ExtensionCodec.
*/
var ExtData = class {
	constructor(type, data) {
		this.type = type;
		this.data = data;
	}
};
//#endregion
//#region ../../work/webPORPID/node_modules/@msgpack/msgpack/dist.esm/DecodeError.mjs
var DecodeError = class DecodeError extends Error {
	constructor(message) {
		super(message);
		const proto = Object.create(DecodeError.prototype);
		Object.setPrototypeOf(this, proto);
		Object.defineProperty(this, "name", {
			configurable: true,
			enumerable: false,
			value: DecodeError.name
		});
	}
};
function setUint64(view, offset, value) {
	const high = value / 4294967296;
	const low = value;
	view.setUint32(offset, high);
	view.setUint32(offset + 4, low);
}
function setInt64(view, offset, value) {
	const high = Math.floor(value / 4294967296);
	const low = value;
	view.setUint32(offset, high);
	view.setUint32(offset + 4, low);
}
function getInt64(view, offset) {
	const high = view.getInt32(offset);
	const low = view.getUint32(offset + 4);
	return high * 4294967296 + low;
}
function getUint64(view, offset) {
	const high = view.getUint32(offset);
	const low = view.getUint32(offset + 4);
	return high * 4294967296 + low;
}
const TIMESTAMP32_MAX_SEC = 4294967295;
const TIMESTAMP64_MAX_SEC = 17179869183;
function encodeTimeSpecToTimestamp({ sec, nsec }) {
	if (sec >= 0 && nsec >= 0 && sec <= TIMESTAMP64_MAX_SEC) if (nsec === 0 && sec <= TIMESTAMP32_MAX_SEC) {
		const rv = new Uint8Array(4);
		new DataView(rv.buffer).setUint32(0, sec);
		return rv;
	} else {
		const secHigh = sec / 4294967296;
		const secLow = sec & 4294967295;
		const rv = new Uint8Array(8);
		const view = new DataView(rv.buffer);
		view.setUint32(0, nsec << 2 | secHigh & 3);
		view.setUint32(4, secLow);
		return rv;
	}
	else {
		const rv = new Uint8Array(12);
		const view = new DataView(rv.buffer);
		view.setUint32(0, nsec);
		setInt64(view, 4, sec);
		return rv;
	}
}
function encodeDateToTimeSpec(date) {
	const msec = date.getTime();
	const sec = Math.floor(msec / 1e3);
	const nsec = (msec - sec * 1e3) * 1e6;
	const nsecInSec = Math.floor(nsec / 1e9);
	return {
		sec: sec + nsecInSec,
		nsec: nsec - nsecInSec * 1e9
	};
}
function encodeTimestampExtension(object) {
	if (object instanceof Date) return encodeTimeSpecToTimestamp(encodeDateToTimeSpec(object));
	else return null;
}
function decodeTimestampToTimeSpec(data) {
	const view = new DataView(data.buffer, data.byteOffset, data.byteLength);
	switch (data.byteLength) {
		case 4: return {
			sec: view.getUint32(0),
			nsec: 0
		};
		case 8: {
			const nsec30AndSecHigh2 = view.getUint32(0);
			const secLow32 = view.getUint32(4);
			return {
				sec: (nsec30AndSecHigh2 & 3) * 4294967296 + secLow32,
				nsec: nsec30AndSecHigh2 >>> 2
			};
		}
		case 12: return {
			sec: getInt64(view, 4),
			nsec: view.getUint32(0)
		};
		default: throw new DecodeError(`Unrecognized data size for timestamp (expected 4, 8, or 12): ${data.length}`);
	}
}
function decodeTimestampExtension(data) {
	const timeSpec = decodeTimestampToTimeSpec(data);
	return /* @__PURE__ */ new Date(timeSpec.sec * 1e3 + timeSpec.nsec / 1e6);
}
const timestampExtension = {
	type: -1,
	encode: encodeTimestampExtension,
	decode: decodeTimestampExtension
};
//#endregion
//#region ../../work/webPORPID/node_modules/@msgpack/msgpack/dist.esm/ExtensionCodec.mjs
var ExtensionCodec = class {
	constructor() {
		this.builtInEncoders = [];
		this.builtInDecoders = [];
		this.encoders = [];
		this.decoders = [];
		this.register(timestampExtension);
	}
	register({ type, encode, decode }) {
		if (type >= 0) {
			this.encoders[type] = encode;
			this.decoders[type] = decode;
		} else {
			const index = -1 - type;
			this.builtInEncoders[index] = encode;
			this.builtInDecoders[index] = decode;
		}
	}
	tryToEncode(object, context) {
		for (let i = 0; i < this.builtInEncoders.length; i++) {
			const encodeExt = this.builtInEncoders[i];
			if (encodeExt != null) {
				const data = encodeExt(object, context);
				if (data != null) return new ExtData(-1 - i, data);
			}
		}
		for (let i = 0; i < this.encoders.length; i++) {
			const encodeExt = this.encoders[i];
			if (encodeExt != null) {
				const data = encodeExt(object, context);
				if (data != null) return new ExtData(i, data);
			}
		}
		if (object instanceof ExtData) return object;
		return null;
	}
	decode(data, type, context) {
		const decodeExt = type < 0 ? this.builtInDecoders[-1 - type] : this.decoders[type];
		if (decodeExt) return decodeExt(data, type, context);
		else return new ExtData(type, data);
	}
};
ExtensionCodec.defaultCodec = new ExtensionCodec();
//#endregion
//#region ../../work/webPORPID/node_modules/@msgpack/msgpack/dist.esm/utils/typedArrays.mjs
function isArrayBufferLike(buffer) {
	return buffer instanceof ArrayBuffer || typeof SharedArrayBuffer !== "undefined" && buffer instanceof SharedArrayBuffer;
}
function ensureUint8Array(buffer) {
	if (buffer instanceof Uint8Array) return buffer;
	else if (ArrayBuffer.isView(buffer)) return new Uint8Array(buffer.buffer, buffer.byteOffset, buffer.byteLength);
	else if (isArrayBufferLike(buffer)) return new Uint8Array(buffer);
	else return Uint8Array.from(buffer);
}
var Encoder = class Encoder {
	constructor(options) {
		this.entered = false;
		this.extensionCodec = options?.extensionCodec ?? ExtensionCodec.defaultCodec;
		this.context = options?.context;
		this.useBigInt64 = options?.useBigInt64 ?? false;
		this.maxDepth = options?.maxDepth ?? 100;
		this.initialBufferSize = options?.initialBufferSize ?? 2048;
		this.sortKeys = options?.sortKeys ?? false;
		this.forceFloat32 = options?.forceFloat32 ?? false;
		this.ignoreUndefined = options?.ignoreUndefined ?? false;
		this.forceIntegerToFloat = options?.forceIntegerToFloat ?? false;
		this.pos = 0;
		this.view = new DataView(new ArrayBuffer(this.initialBufferSize));
		this.bytes = new Uint8Array(this.view.buffer);
	}
	clone() {
		return new Encoder({
			extensionCodec: this.extensionCodec,
			context: this.context,
			useBigInt64: this.useBigInt64,
			maxDepth: this.maxDepth,
			initialBufferSize: this.initialBufferSize,
			sortKeys: this.sortKeys,
			forceFloat32: this.forceFloat32,
			ignoreUndefined: this.ignoreUndefined,
			forceIntegerToFloat: this.forceIntegerToFloat
		});
	}
	reinitializeState() {
		this.pos = 0;
	}
	/**
	* This is almost equivalent to {@link Encoder#encode}, but it returns an reference of the encoder's internal buffer and thus much faster than {@link Encoder#encode}.
	*
	* @returns Encodes the object and returns a shared reference the encoder's internal buffer.
	*/
	encodeSharedRef(object) {
		if (this.entered) return this.clone().encodeSharedRef(object);
		try {
			this.entered = true;
			this.reinitializeState();
			this.doEncode(object, 1);
			return this.bytes.subarray(0, this.pos);
		} finally {
			this.entered = false;
		}
	}
	/**
	* @returns Encodes the object and returns a copy of the encoder's internal buffer.
	*/
	encode(object) {
		if (this.entered) return this.clone().encode(object);
		try {
			this.entered = true;
			this.reinitializeState();
			this.doEncode(object, 1);
			return this.bytes.slice(0, this.pos);
		} finally {
			this.entered = false;
		}
	}
	doEncode(object, depth) {
		if (depth > this.maxDepth) throw new Error(`Too deep objects in depth ${depth}`);
		if (object == null) this.encodeNil();
		else if (typeof object === "boolean") this.encodeBoolean(object);
		else if (typeof object === "number") if (!this.forceIntegerToFloat) this.encodeNumber(object);
		else this.encodeNumberAsFloat(object);
		else if (typeof object === "string") this.encodeString(object);
		else if (this.useBigInt64 && typeof object === "bigint") this.encodeBigInt64(object);
		else this.encodeObject(object, depth);
	}
	ensureBufferSizeToWrite(sizeToWrite) {
		const requiredSize = this.pos + sizeToWrite;
		if (this.view.byteLength < requiredSize) this.resizeBuffer(requiredSize * 2);
	}
	resizeBuffer(newSize) {
		const newBuffer = new ArrayBuffer(newSize);
		const newBytes = new Uint8Array(newBuffer);
		const newView = new DataView(newBuffer);
		newBytes.set(this.bytes);
		this.view = newView;
		this.bytes = newBytes;
	}
	encodeNil() {
		this.writeU8(192);
	}
	encodeBoolean(object) {
		if (object === false) this.writeU8(194);
		else this.writeU8(195);
	}
	encodeNumber(object) {
		if (!this.forceIntegerToFloat && Number.isSafeInteger(object)) if (object >= 0) if (object < 128) this.writeU8(object);
		else if (object < 256) {
			this.writeU8(204);
			this.writeU8(object);
		} else if (object < 65536) {
			this.writeU8(205);
			this.writeU16(object);
		} else if (object < 4294967296) {
			this.writeU8(206);
			this.writeU32(object);
		} else if (!this.useBigInt64) {
			this.writeU8(207);
			this.writeU64(object);
		} else this.encodeNumberAsFloat(object);
		else if (object >= -32) this.writeU8(224 | object + 32);
		else if (object >= -128) {
			this.writeU8(208);
			this.writeI8(object);
		} else if (object >= -32768) {
			this.writeU8(209);
			this.writeI16(object);
		} else if (object >= -2147483648) {
			this.writeU8(210);
			this.writeI32(object);
		} else if (!this.useBigInt64) {
			this.writeU8(211);
			this.writeI64(object);
		} else this.encodeNumberAsFloat(object);
		else this.encodeNumberAsFloat(object);
	}
	encodeNumberAsFloat(object) {
		if (this.forceFloat32) {
			this.writeU8(202);
			this.writeF32(object);
		} else {
			this.writeU8(203);
			this.writeF64(object);
		}
	}
	encodeBigInt64(object) {
		if (object >= BigInt(0)) {
			this.writeU8(207);
			this.writeBigUint64(object);
		} else {
			this.writeU8(211);
			this.writeBigInt64(object);
		}
	}
	writeStringHeader(byteLength) {
		if (byteLength < 32) this.writeU8(160 + byteLength);
		else if (byteLength < 256) {
			this.writeU8(217);
			this.writeU8(byteLength);
		} else if (byteLength < 65536) {
			this.writeU8(218);
			this.writeU16(byteLength);
		} else if (byteLength < 4294967296) {
			this.writeU8(219);
			this.writeU32(byteLength);
		} else throw new Error(`Too long string: ${byteLength} bytes in UTF-8`);
	}
	encodeString(object) {
		const maxHeaderSize = 5;
		const byteLength = utf8Count(object);
		this.ensureBufferSizeToWrite(maxHeaderSize + byteLength);
		this.writeStringHeader(byteLength);
		utf8Encode(object, this.bytes, this.pos);
		this.pos += byteLength;
	}
	encodeObject(object, depth) {
		const ext = this.extensionCodec.tryToEncode(object, this.context);
		if (ext != null) this.encodeExtension(ext);
		else if (Array.isArray(object)) this.encodeArray(object, depth);
		else if (ArrayBuffer.isView(object)) this.encodeBinary(object);
		else if (typeof object === "object") this.encodeMap(object, depth);
		else throw new Error(`Unrecognized object: ${Object.prototype.toString.apply(object)}`);
	}
	encodeBinary(object) {
		const size = object.byteLength;
		if (size < 256) {
			this.writeU8(196);
			this.writeU8(size);
		} else if (size < 65536) {
			this.writeU8(197);
			this.writeU16(size);
		} else if (size < 4294967296) {
			this.writeU8(198);
			this.writeU32(size);
		} else throw new Error(`Too large binary: ${size}`);
		const bytes = ensureUint8Array(object);
		this.writeU8a(bytes);
	}
	encodeArray(object, depth) {
		const size = object.length;
		if (size < 16) this.writeU8(144 + size);
		else if (size < 65536) {
			this.writeU8(220);
			this.writeU16(size);
		} else if (size < 4294967296) {
			this.writeU8(221);
			this.writeU32(size);
		} else throw new Error(`Too large array: ${size}`);
		for (const item of object) this.doEncode(item, depth + 1);
	}
	countWithoutUndefined(object, keys) {
		let count = 0;
		for (const key of keys) if (object[key] !== void 0) count++;
		return count;
	}
	encodeMap(object, depth) {
		const keys = Object.keys(object);
		if (this.sortKeys) keys.sort();
		const size = this.ignoreUndefined ? this.countWithoutUndefined(object, keys) : keys.length;
		if (size < 16) this.writeU8(128 + size);
		else if (size < 65536) {
			this.writeU8(222);
			this.writeU16(size);
		} else if (size < 4294967296) {
			this.writeU8(223);
			this.writeU32(size);
		} else throw new Error(`Too large map object: ${size}`);
		for (const key of keys) {
			const value = object[key];
			if (!(this.ignoreUndefined && value === void 0)) {
				this.encodeString(key);
				this.doEncode(value, depth + 1);
			}
		}
	}
	encodeExtension(ext) {
		if (typeof ext.data === "function") {
			const data = ext.data(this.pos + 6);
			const size = data.length;
			if (size >= 4294967296) throw new Error(`Too large extension object: ${size}`);
			this.writeU8(201);
			this.writeU32(size);
			this.writeI8(ext.type);
			this.writeU8a(data);
			return;
		}
		const size = ext.data.length;
		if (size === 1) this.writeU8(212);
		else if (size === 2) this.writeU8(213);
		else if (size === 4) this.writeU8(214);
		else if (size === 8) this.writeU8(215);
		else if (size === 16) this.writeU8(216);
		else if (size < 256) {
			this.writeU8(199);
			this.writeU8(size);
		} else if (size < 65536) {
			this.writeU8(200);
			this.writeU16(size);
		} else if (size < 4294967296) {
			this.writeU8(201);
			this.writeU32(size);
		} else throw new Error(`Too large extension object: ${size}`);
		this.writeI8(ext.type);
		this.writeU8a(ext.data);
	}
	writeU8(value) {
		this.ensureBufferSizeToWrite(1);
		this.view.setUint8(this.pos, value);
		this.pos++;
	}
	writeU8a(values) {
		const size = values.length;
		this.ensureBufferSizeToWrite(size);
		this.bytes.set(values, this.pos);
		this.pos += size;
	}
	writeI8(value) {
		this.ensureBufferSizeToWrite(1);
		this.view.setInt8(this.pos, value);
		this.pos++;
	}
	writeU16(value) {
		this.ensureBufferSizeToWrite(2);
		this.view.setUint16(this.pos, value);
		this.pos += 2;
	}
	writeI16(value) {
		this.ensureBufferSizeToWrite(2);
		this.view.setInt16(this.pos, value);
		this.pos += 2;
	}
	writeU32(value) {
		this.ensureBufferSizeToWrite(4);
		this.view.setUint32(this.pos, value);
		this.pos += 4;
	}
	writeI32(value) {
		this.ensureBufferSizeToWrite(4);
		this.view.setInt32(this.pos, value);
		this.pos += 4;
	}
	writeF32(value) {
		this.ensureBufferSizeToWrite(4);
		this.view.setFloat32(this.pos, value);
		this.pos += 4;
	}
	writeF64(value) {
		this.ensureBufferSizeToWrite(8);
		this.view.setFloat64(this.pos, value);
		this.pos += 8;
	}
	writeU64(value) {
		this.ensureBufferSizeToWrite(8);
		setUint64(this.view, this.pos, value);
		this.pos += 8;
	}
	writeI64(value) {
		this.ensureBufferSizeToWrite(8);
		setInt64(this.view, this.pos, value);
		this.pos += 8;
	}
	writeBigUint64(value) {
		this.ensureBufferSizeToWrite(8);
		this.view.setBigUint64(this.pos, value);
		this.pos += 8;
	}
	writeBigInt64(value) {
		this.ensureBufferSizeToWrite(8);
		this.view.setBigInt64(this.pos, value);
		this.pos += 8;
	}
};
//#endregion
//#region ../../work/webPORPID/node_modules/@msgpack/msgpack/dist.esm/encode.mjs
/**
* It encodes `value` in the MessagePack format and
* returns a byte buffer.
*
* The returned buffer is a slice of a larger `ArrayBuffer`, so you have to use its `#byteOffset` and `#byteLength` in order to convert it to another typed arrays including NodeJS `Buffer`.
*/
function encode(value, options) {
	return new Encoder(options).encodeSharedRef(value);
}
//#endregion
//#region ../../work/webPORPID/node_modules/@msgpack/msgpack/dist.esm/utils/prettyByte.mjs
function prettyByte(byte) {
	return `${byte < 0 ? "-" : ""}0x${Math.abs(byte).toString(16).padStart(2, "0")}`;
}
//#endregion
//#region ../../work/webPORPID/node_modules/@msgpack/msgpack/dist.esm/CachedKeyDecoder.mjs
const DEFAULT_MAX_KEY_LENGTH = 16;
const DEFAULT_MAX_LENGTH_PER_KEY = 16;
var CachedKeyDecoder = class {
	constructor(maxKeyLength = DEFAULT_MAX_KEY_LENGTH, maxLengthPerKey = DEFAULT_MAX_LENGTH_PER_KEY) {
		this.hit = 0;
		this.miss = 0;
		this.maxKeyLength = maxKeyLength;
		this.maxLengthPerKey = maxLengthPerKey;
		this.caches = [];
		for (let i = 0; i < this.maxKeyLength; i++) this.caches.push([]);
	}
	canBeCached(byteLength) {
		return byteLength > 0 && byteLength <= this.maxKeyLength;
	}
	find(bytes, inputOffset, byteLength) {
		const records = this.caches[byteLength - 1];
		FIND_CHUNK: for (const record of records) {
			const recordBytes = record.bytes;
			for (let j = 0; j < byteLength; j++) if (recordBytes[j] !== bytes[inputOffset + j]) continue FIND_CHUNK;
			return record.str;
		}
		return null;
	}
	store(bytes, value) {
		const records = this.caches[bytes.length - 1];
		const record = {
			bytes,
			str: value
		};
		if (records.length >= this.maxLengthPerKey) records[Math.random() * records.length | 0] = record;
		else records.push(record);
	}
	decode(bytes, inputOffset, byteLength) {
		const cachedValue = this.find(bytes, inputOffset, byteLength);
		if (cachedValue != null) {
			this.hit++;
			return cachedValue;
		}
		this.miss++;
		const str = utf8DecodeJs(bytes, inputOffset, byteLength);
		const slicedCopyOfBytes = Uint8Array.prototype.slice.call(bytes, inputOffset, inputOffset + byteLength);
		this.store(slicedCopyOfBytes, str);
		return str;
	}
};
//#endregion
//#region ../../work/webPORPID/node_modules/@msgpack/msgpack/dist.esm/Decoder.mjs
const STATE_ARRAY = "array";
const STATE_MAP_KEY = "map_key";
const STATE_MAP_VALUE = "map_value";
const mapKeyConverter = (key) => {
	if (typeof key === "string" || typeof key === "number") return key;
	throw new DecodeError("The type of key must be string or number but " + typeof key);
};
var StackPool = class {
	constructor() {
		this.stack = [];
		this.stackHeadPosition = -1;
	}
	get length() {
		return this.stackHeadPosition + 1;
	}
	top() {
		return this.stack[this.stackHeadPosition];
	}
	pushArrayState(size) {
		const state = this.getUninitializedStateFromPool();
		state.type = STATE_ARRAY;
		state.position = 0;
		state.size = size;
		state.array = new Array(size);
	}
	pushMapState(size) {
		const state = this.getUninitializedStateFromPool();
		state.type = STATE_MAP_KEY;
		state.readCount = 0;
		state.size = size;
		state.map = {};
	}
	getUninitializedStateFromPool() {
		this.stackHeadPosition++;
		if (this.stackHeadPosition === this.stack.length) this.stack.push({
			type: void 0,
			size: 0,
			array: void 0,
			position: 0,
			readCount: 0,
			map: void 0,
			key: null
		});
		return this.stack[this.stackHeadPosition];
	}
	release(state) {
		if (this.stack[this.stackHeadPosition] !== state) throw new Error("Invalid stack state. Released state is not on top of the stack.");
		if (state.type === STATE_ARRAY) {
			const partialState = state;
			partialState.size = 0;
			partialState.array = void 0;
			partialState.position = 0;
			partialState.type = void 0;
		}
		if (state.type === STATE_MAP_KEY || state.type === STATE_MAP_VALUE) {
			const partialState = state;
			partialState.size = 0;
			partialState.map = void 0;
			partialState.readCount = 0;
			partialState.type = void 0;
		}
		this.stackHeadPosition--;
	}
	reset() {
		this.stack.length = 0;
		this.stackHeadPosition = -1;
	}
};
const HEAD_BYTE_REQUIRED = -1;
const EMPTY_VIEW = /* @__PURE__ */ new DataView(/* @__PURE__ */ new ArrayBuffer(0));
const EMPTY_BYTES = new Uint8Array(EMPTY_VIEW.buffer);
try {
	EMPTY_VIEW.getInt8(0);
} catch (e) {
	if (!(e instanceof RangeError)) throw new Error("This module is not supported in the current JavaScript engine because DataView does not throw RangeError on out-of-bounds access");
}
const MORE_DATA = /* @__PURE__ */ new RangeError("Insufficient data");
const sharedCachedKeyDecoder = new CachedKeyDecoder();
var Decoder = class Decoder {
	constructor(options) {
		this.totalPos = 0;
		this.pos = 0;
		this.view = EMPTY_VIEW;
		this.bytes = EMPTY_BYTES;
		this.headByte = HEAD_BYTE_REQUIRED;
		this.stack = new StackPool();
		this.entered = false;
		this.extensionCodec = options?.extensionCodec ?? ExtensionCodec.defaultCodec;
		this.context = options?.context;
		this.useBigInt64 = options?.useBigInt64 ?? false;
		this.rawStrings = options?.rawStrings ?? false;
		this.maxStrLength = options?.maxStrLength ?? 4294967295;
		this.maxBinLength = options?.maxBinLength ?? 4294967295;
		this.maxArrayLength = options?.maxArrayLength ?? 4294967295;
		this.maxMapLength = options?.maxMapLength ?? 4294967295;
		this.maxExtLength = options?.maxExtLength ?? 4294967295;
		this.keyDecoder = options?.keyDecoder !== void 0 ? options.keyDecoder : sharedCachedKeyDecoder;
		this.mapKeyConverter = options?.mapKeyConverter ?? mapKeyConverter;
	}
	clone() {
		return new Decoder({
			extensionCodec: this.extensionCodec,
			context: this.context,
			useBigInt64: this.useBigInt64,
			rawStrings: this.rawStrings,
			maxStrLength: this.maxStrLength,
			maxBinLength: this.maxBinLength,
			maxArrayLength: this.maxArrayLength,
			maxMapLength: this.maxMapLength,
			maxExtLength: this.maxExtLength,
			keyDecoder: this.keyDecoder
		});
	}
	reinitializeState() {
		this.totalPos = 0;
		this.headByte = HEAD_BYTE_REQUIRED;
		this.stack.reset();
	}
	setBuffer(buffer) {
		const bytes = ensureUint8Array(buffer);
		this.bytes = bytes;
		this.view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
		this.pos = 0;
	}
	appendBuffer(buffer) {
		if (this.headByte === HEAD_BYTE_REQUIRED && !this.hasRemaining(1)) this.setBuffer(buffer);
		else {
			const remainingData = this.bytes.subarray(this.pos);
			const newData = ensureUint8Array(buffer);
			const newBuffer = new Uint8Array(remainingData.length + newData.length);
			newBuffer.set(remainingData);
			newBuffer.set(newData, remainingData.length);
			this.setBuffer(newBuffer);
		}
	}
	hasRemaining(size) {
		return this.view.byteLength - this.pos >= size;
	}
	createExtraByteError(posToShow) {
		const { view, pos } = this;
		return /* @__PURE__ */ new RangeError(`Extra ${view.byteLength - pos} of ${view.byteLength} byte(s) found at buffer[${posToShow}]`);
	}
	/**
	* @throws {@link DecodeError}
	* @throws {@link RangeError}
	*/
	decode(buffer) {
		if (this.entered) return this.clone().decode(buffer);
		try {
			this.entered = true;
			this.reinitializeState();
			this.setBuffer(buffer);
			const object = this.doDecodeSync();
			if (this.hasRemaining(1)) throw this.createExtraByteError(this.pos);
			return object;
		} finally {
			this.entered = false;
		}
	}
	*decodeMulti(buffer) {
		if (this.entered) {
			yield* this.clone().decodeMulti(buffer);
			return;
		}
		try {
			this.entered = true;
			this.reinitializeState();
			this.setBuffer(buffer);
			while (this.hasRemaining(1)) yield this.doDecodeSync();
		} finally {
			this.entered = false;
		}
	}
	async decodeAsync(stream) {
		if (this.entered) return this.clone().decodeAsync(stream);
		try {
			this.entered = true;
			let decoded = false;
			let object;
			for await (const buffer of stream) {
				if (decoded) {
					this.entered = false;
					throw this.createExtraByteError(this.totalPos);
				}
				this.appendBuffer(buffer);
				try {
					object = this.doDecodeSync();
					decoded = true;
				} catch (e) {
					if (!(e instanceof RangeError)) throw e;
				}
				this.totalPos += this.pos;
			}
			if (decoded) {
				if (this.hasRemaining(1)) throw this.createExtraByteError(this.totalPos);
				return object;
			}
			const { headByte, pos, totalPos } = this;
			throw new RangeError(`Insufficient data in parsing ${prettyByte(headByte)} at ${totalPos} (${pos} in the current buffer)`);
		} finally {
			this.entered = false;
		}
	}
	decodeArrayStream(stream) {
		return this.decodeMultiAsync(stream, true);
	}
	decodeStream(stream) {
		return this.decodeMultiAsync(stream, false);
	}
	async *decodeMultiAsync(stream, isArray) {
		if (this.entered) {
			yield* this.clone().decodeMultiAsync(stream, isArray);
			return;
		}
		try {
			this.entered = true;
			let isArrayHeaderRequired = isArray;
			let arrayItemsLeft = -1;
			for await (const buffer of stream) {
				if (isArray && arrayItemsLeft === 0) throw this.createExtraByteError(this.totalPos);
				this.appendBuffer(buffer);
				if (isArrayHeaderRequired) {
					arrayItemsLeft = this.readArraySize();
					isArrayHeaderRequired = false;
					this.complete();
				}
				try {
					while (true) {
						yield this.doDecodeSync();
						if (--arrayItemsLeft === 0) break;
					}
				} catch (e) {
					if (!(e instanceof RangeError)) throw e;
				}
				this.totalPos += this.pos;
			}
		} finally {
			this.entered = false;
		}
	}
	doDecodeSync() {
		DECODE: while (true) {
			const headByte = this.readHeadByte();
			let object;
			if (headByte >= 224) object = headByte - 256;
			else if (headByte < 192) if (headByte < 128) object = headByte;
			else if (headByte < 144) {
				const size = headByte - 128;
				if (size !== 0) {
					this.pushMapState(size);
					this.complete();
					continue DECODE;
				} else object = {};
			} else if (headByte < 160) {
				const size = headByte - 144;
				if (size !== 0) {
					this.pushArrayState(size);
					this.complete();
					continue DECODE;
				} else object = [];
			} else {
				const byteLength = headByte - 160;
				object = this.decodeString(byteLength, 0);
			}
			else if (headByte === 192) object = null;
			else if (headByte === 194) object = false;
			else if (headByte === 195) object = true;
			else if (headByte === 202) object = this.readF32();
			else if (headByte === 203) object = this.readF64();
			else if (headByte === 204) object = this.readU8();
			else if (headByte === 205) object = this.readU16();
			else if (headByte === 206) object = this.readU32();
			else if (headByte === 207) if (this.useBigInt64) object = this.readU64AsBigInt();
			else object = this.readU64();
			else if (headByte === 208) object = this.readI8();
			else if (headByte === 209) object = this.readI16();
			else if (headByte === 210) object = this.readI32();
			else if (headByte === 211) if (this.useBigInt64) object = this.readI64AsBigInt();
			else object = this.readI64();
			else if (headByte === 217) {
				const byteLength = this.lookU8();
				object = this.decodeString(byteLength, 1);
			} else if (headByte === 218) {
				const byteLength = this.lookU16();
				object = this.decodeString(byteLength, 2);
			} else if (headByte === 219) {
				const byteLength = this.lookU32();
				object = this.decodeString(byteLength, 4);
			} else if (headByte === 220) {
				const size = this.readU16();
				if (size !== 0) {
					this.pushArrayState(size);
					this.complete();
					continue DECODE;
				} else object = [];
			} else if (headByte === 221) {
				const size = this.readU32();
				if (size !== 0) {
					this.pushArrayState(size);
					this.complete();
					continue DECODE;
				} else object = [];
			} else if (headByte === 222) {
				const size = this.readU16();
				if (size !== 0) {
					this.pushMapState(size);
					this.complete();
					continue DECODE;
				} else object = {};
			} else if (headByte === 223) {
				const size = this.readU32();
				if (size !== 0) {
					this.pushMapState(size);
					this.complete();
					continue DECODE;
				} else object = {};
			} else if (headByte === 196) {
				const size = this.lookU8();
				object = this.decodeBinary(size, 1);
			} else if (headByte === 197) {
				const size = this.lookU16();
				object = this.decodeBinary(size, 2);
			} else if (headByte === 198) {
				const size = this.lookU32();
				object = this.decodeBinary(size, 4);
			} else if (headByte === 212) object = this.decodeExtension(1, 0);
			else if (headByte === 213) object = this.decodeExtension(2, 0);
			else if (headByte === 214) object = this.decodeExtension(4, 0);
			else if (headByte === 215) object = this.decodeExtension(8, 0);
			else if (headByte === 216) object = this.decodeExtension(16, 0);
			else if (headByte === 199) {
				const size = this.lookU8();
				object = this.decodeExtension(size, 1);
			} else if (headByte === 200) {
				const size = this.lookU16();
				object = this.decodeExtension(size, 2);
			} else if (headByte === 201) {
				const size = this.lookU32();
				object = this.decodeExtension(size, 4);
			} else throw new DecodeError(`Unrecognized type byte: ${prettyByte(headByte)}`);
			this.complete();
			const stack = this.stack;
			while (stack.length > 0) {
				const state = stack.top();
				if (state.type === STATE_ARRAY) {
					state.array[state.position] = object;
					state.position++;
					if (state.position === state.size) {
						object = state.array;
						stack.release(state);
					} else continue DECODE;
				} else if (state.type === STATE_MAP_KEY) {
					if (object === "__proto__") throw new DecodeError("The key __proto__ is not allowed");
					state.key = this.mapKeyConverter(object);
					state.type = STATE_MAP_VALUE;
					continue DECODE;
				} else {
					state.map[state.key] = object;
					state.readCount++;
					if (state.readCount === state.size) {
						object = state.map;
						stack.release(state);
					} else {
						state.key = null;
						state.type = STATE_MAP_KEY;
						continue DECODE;
					}
				}
			}
			return object;
		}
	}
	readHeadByte() {
		if (this.headByte === HEAD_BYTE_REQUIRED) this.headByte = this.readU8();
		return this.headByte;
	}
	complete() {
		this.headByte = HEAD_BYTE_REQUIRED;
	}
	readArraySize() {
		const headByte = this.readHeadByte();
		switch (headByte) {
			case 220: return this.readU16();
			case 221: return this.readU32();
			default: if (headByte < 160) return headByte - 144;
			else throw new DecodeError(`Unrecognized array type byte: ${prettyByte(headByte)}`);
		}
	}
	pushMapState(size) {
		if (size > this.maxMapLength) throw new DecodeError(`Max length exceeded: map length (${size}) > maxMapLengthLength (${this.maxMapLength})`);
		this.stack.pushMapState(size);
	}
	pushArrayState(size) {
		if (size > this.maxArrayLength) throw new DecodeError(`Max length exceeded: array length (${size}) > maxArrayLength (${this.maxArrayLength})`);
		this.stack.pushArrayState(size);
	}
	decodeString(byteLength, headerOffset) {
		if (!this.rawStrings || this.stateIsMapKey()) return this.decodeUtf8String(byteLength, headerOffset);
		return this.decodeBinary(byteLength, headerOffset);
	}
	/**
	* @throws {@link RangeError}
	*/
	decodeUtf8String(byteLength, headerOffset) {
		if (byteLength > this.maxStrLength) throw new DecodeError(`Max length exceeded: UTF-8 byte length (${byteLength}) > maxStrLength (${this.maxStrLength})`);
		if (this.bytes.byteLength < this.pos + headerOffset + byteLength) throw MORE_DATA;
		const offset = this.pos + headerOffset;
		let object;
		if (this.stateIsMapKey() && this.keyDecoder?.canBeCached(byteLength)) object = this.keyDecoder.decode(this.bytes, offset, byteLength);
		else object = utf8Decode(this.bytes, offset, byteLength);
		this.pos += headerOffset + byteLength;
		return object;
	}
	stateIsMapKey() {
		if (this.stack.length > 0) return this.stack.top().type === STATE_MAP_KEY;
		return false;
	}
	/**
	* @throws {@link RangeError}
	*/
	decodeBinary(byteLength, headOffset) {
		if (byteLength > this.maxBinLength) throw new DecodeError(`Max length exceeded: bin length (${byteLength}) > maxBinLength (${this.maxBinLength})`);
		if (!this.hasRemaining(byteLength + headOffset)) throw MORE_DATA;
		const offset = this.pos + headOffset;
		const object = this.bytes.subarray(offset, offset + byteLength);
		this.pos += headOffset + byteLength;
		return object;
	}
	decodeExtension(size, headOffset) {
		if (size > this.maxExtLength) throw new DecodeError(`Max length exceeded: ext length (${size}) > maxExtLength (${this.maxExtLength})`);
		const extType = this.view.getInt8(this.pos + headOffset);
		const data = this.decodeBinary(size, headOffset + 1);
		return this.extensionCodec.decode(data, extType, this.context);
	}
	lookU8() {
		return this.view.getUint8(this.pos);
	}
	lookU16() {
		return this.view.getUint16(this.pos);
	}
	lookU32() {
		return this.view.getUint32(this.pos);
	}
	readU8() {
		const value = this.view.getUint8(this.pos);
		this.pos++;
		return value;
	}
	readI8() {
		const value = this.view.getInt8(this.pos);
		this.pos++;
		return value;
	}
	readU16() {
		const value = this.view.getUint16(this.pos);
		this.pos += 2;
		return value;
	}
	readI16() {
		const value = this.view.getInt16(this.pos);
		this.pos += 2;
		return value;
	}
	readU32() {
		const value = this.view.getUint32(this.pos);
		this.pos += 4;
		return value;
	}
	readI32() {
		const value = this.view.getInt32(this.pos);
		this.pos += 4;
		return value;
	}
	readU64() {
		const value = getUint64(this.view, this.pos);
		this.pos += 8;
		return value;
	}
	readI64() {
		const value = getInt64(this.view, this.pos);
		this.pos += 8;
		return value;
	}
	readU64AsBigInt() {
		const value = this.view.getBigUint64(this.pos);
		this.pos += 8;
		return value;
	}
	readI64AsBigInt() {
		const value = this.view.getBigInt64(this.pos);
		this.pos += 8;
		return value;
	}
	readF32() {
		const value = this.view.getFloat32(this.pos);
		this.pos += 4;
		return value;
	}
	readF64() {
		const value = this.view.getFloat64(this.pos);
		this.pos += 8;
		return value;
	}
};
//#endregion
//#region ../../work/webPORPID/node_modules/@msgpack/msgpack/dist.esm/decode.mjs
/**
* It decodes a single MessagePack object in a buffer.
*
* This is a synchronous decoding function.
* See other variants for asynchronous decoding: {@link decodeAsync}, {@link decodeMultiStream}, or {@link decodeArrayStream}.
*
* @throws {@link RangeError} if the buffer is incomplete, including the case where the buffer is empty.
* @throws {@link DecodeError} if the buffer contains invalid data.
*/
function decode(buffer, options) {
	return new Decoder(options).decode(buffer);
}
//#endregion
//#region ../../work/webPORPID/node_modules/fflate/esm/browser.js
var u8 = Uint8Array, u16 = Uint16Array, i32 = Int32Array;
var fleb = new u8([
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	1,
	1,
	1,
	1,
	2,
	2,
	2,
	2,
	3,
	3,
	3,
	3,
	4,
	4,
	4,
	4,
	5,
	5,
	5,
	5,
	0,
	0,
	0,
	0
]);
var fdeb = new u8([
	0,
	0,
	0,
	0,
	1,
	1,
	2,
	2,
	3,
	3,
	4,
	4,
	5,
	5,
	6,
	6,
	7,
	7,
	8,
	8,
	9,
	9,
	10,
	10,
	11,
	11,
	12,
	12,
	13,
	13,
	0,
	0
]);
var clim = new u8([
	16,
	17,
	18,
	0,
	8,
	7,
	9,
	6,
	10,
	5,
	11,
	4,
	12,
	3,
	13,
	2,
	14,
	1,
	15
]);
var freb = function(eb, start) {
	var b = new u16(31);
	for (var i = 0; i < 31; ++i) b[i] = start += 1 << eb[i - 1];
	var r = new i32(b[30]);
	for (var i = 1; i < 30; ++i) for (var j = b[i]; j < b[i + 1]; ++j) r[j] = j - b[i] << 5 | i;
	return {
		b,
		r
	};
};
var _a = freb(fleb, 2), fl = _a.b, revfl = _a.r;
fl[28] = 258, revfl[258] = 28;
var _b = freb(fdeb, 0), fd = _b.b, revfd = _b.r;
var rev = new u16(32768);
for (var i = 0; i < 32768; ++i) {
	var x = (i & 43690) >> 1 | (i & 21845) << 1;
	x = (x & 52428) >> 2 | (x & 13107) << 2;
	x = (x & 61680) >> 4 | (x & 3855) << 4;
	rev[i] = ((x & 65280) >> 8 | (x & 255) << 8) >> 1;
}
var hMap = (function(cd, mb, r) {
	var s = cd.length;
	var i = 0;
	var l = new u16(mb);
	for (; i < s; ++i) if (cd[i]) ++l[cd[i] - 1];
	var le = new u16(mb);
	for (i = 1; i < mb; ++i) le[i] = le[i - 1] + l[i - 1] << 1;
	var co;
	if (r) {
		co = new u16(1 << mb);
		var rvb = 15 - mb;
		for (i = 0; i < s; ++i) if (cd[i]) {
			var sv = i << 4 | cd[i];
			var r_1 = mb - cd[i];
			var v = le[cd[i] - 1]++ << r_1;
			for (var m = v | (1 << r_1) - 1; v <= m; ++v) co[rev[v] >> rvb] = sv;
		}
	} else {
		co = new u16(s);
		for (i = 0; i < s; ++i) if (cd[i]) co[i] = rev[le[cd[i] - 1]++] >> 15 - cd[i];
	}
	return co;
});
var flt = new u8(288);
for (var i = 0; i < 144; ++i) flt[i] = 8;
for (var i = 144; i < 256; ++i) flt[i] = 9;
for (var i = 256; i < 280; ++i) flt[i] = 7;
for (var i = 280; i < 288; ++i) flt[i] = 8;
var fdt = new u8(32);
for (var i = 0; i < 32; ++i) fdt[i] = 5;
var flm = /* @__PURE__ */ hMap(flt, 9, 0), flrm = /* @__PURE__ */ hMap(flt, 9, 1);
var fdm = /* @__PURE__ */ hMap(fdt, 5, 0), fdrm = /* @__PURE__ */ hMap(fdt, 5, 1);
var max = function(a) {
	var m = a[0];
	for (var i = 1; i < a.length; ++i) if (a[i] > m) m = a[i];
	return m;
};
var bits = function(d, p, m) {
	var o = p / 8 | 0;
	return (d[o] | d[o + 1] << 8) >> (p & 7) & m;
};
var bits16 = function(d, p) {
	var o = p / 8 | 0;
	return (d[o] | d[o + 1] << 8 | d[o + 2] << 16) >> (p & 7);
};
var shft = function(p) {
	return (p + 7) / 8 | 0;
};
var slc = function(v, s, e) {
	if (s == null || s < 0) s = 0;
	if (e == null || e > v.length) e = v.length;
	return new u8(v.subarray(s, e));
};
var ec = [
	"unexpected EOF",
	"invalid block type",
	"invalid length/literal",
	"invalid distance",
	"stream finished",
	"no stream handler",
	,
	"no callback",
	"invalid UTF-8 data",
	"extra field too long",
	"date not in range 1980-2099",
	"filename too long",
	"stream finishing",
	"invalid zip data"
];
var err = function(ind, msg, nt) {
	var e = new Error(msg || ec[ind]);
	e.code = ind;
	if (Error.captureStackTrace) Error.captureStackTrace(e, err);
	if (!nt) throw e;
	return e;
};
var inflt = function(dat, st, buf, dict) {
	var sl = dat.length, dl = dict ? dict.length : 0;
	if (!sl || st.f && !st.l) return buf || new u8(0);
	var noBuf = !buf;
	var resize = noBuf || st.i != 2;
	var noSt = st.i;
	if (noBuf) buf = new u8(sl * 3);
	var cbuf = function(l) {
		var bl = buf.length;
		if (l > bl) {
			var nbuf = new u8(Math.max(bl * 2, l));
			nbuf.set(buf);
			buf = nbuf;
		}
	};
	var final = st.f || 0, pos = st.p || 0, bt = st.b || 0, lm = st.l, dm = st.d, lbt = st.m, dbt = st.n;
	var tbts = sl * 8;
	do {
		if (!lm) {
			final = bits(dat, pos, 1);
			var type = bits(dat, pos + 1, 3);
			pos += 3;
			if (!type) {
				var s = shft(pos) + 4, l = dat[s - 4] | dat[s - 3] << 8, t = s + l;
				if (t > sl) {
					if (noSt) err(0);
					break;
				}
				if (resize) cbuf(bt + l);
				buf.set(dat.subarray(s, t), bt);
				st.b = bt += l, st.p = pos = t * 8, st.f = final;
				continue;
			} else if (type == 1) lm = flrm, dm = fdrm, lbt = 9, dbt = 5;
			else if (type == 2) {
				var hLit = bits(dat, pos, 31) + 257, hcLen = bits(dat, pos + 10, 15) + 4;
				var tl = hLit + bits(dat, pos + 5, 31) + 1;
				pos += 14;
				var ldt = new u8(tl);
				var clt = new u8(19);
				for (var i = 0; i < hcLen; ++i) clt[clim[i]] = bits(dat, pos + i * 3, 7);
				pos += hcLen * 3;
				var clb = max(clt), clbmsk = (1 << clb) - 1;
				var clm = hMap(clt, clb, 1);
				for (var i = 0; i < tl;) {
					var r = clm[bits(dat, pos, clbmsk)];
					pos += r & 15;
					var s = r >> 4;
					if (s < 16) ldt[i++] = s;
					else {
						var c = 0, n = 0;
						if (s == 16) n = 3 + bits(dat, pos, 3), pos += 2, c = ldt[i - 1];
						else if (s == 17) n = 3 + bits(dat, pos, 7), pos += 3;
						else if (s == 18) n = 11 + bits(dat, pos, 127), pos += 7;
						while (n--) ldt[i++] = c;
					}
				}
				var lt = ldt.subarray(0, hLit), dt = ldt.subarray(hLit);
				lbt = max(lt);
				dbt = max(dt);
				lm = hMap(lt, lbt, 1);
				dm = hMap(dt, dbt, 1);
			} else err(1);
			if (pos > tbts) {
				if (noSt) err(0);
				break;
			}
		}
		if (resize) cbuf(bt + 131072);
		var lms = (1 << lbt) - 1, dms = (1 << dbt) - 1;
		var lpos = pos;
		for (;; lpos = pos) {
			var c = lm[bits16(dat, pos) & lms], sym = c >> 4;
			pos += c & 15;
			if (pos > tbts) {
				if (noSt) err(0);
				break;
			}
			if (!c) err(2);
			if (sym < 256) buf[bt++] = sym;
			else if (sym == 256) {
				lpos = pos, lm = null;
				break;
			} else {
				var add = sym - 254;
				if (sym > 264) {
					var i = sym - 257, b = fleb[i];
					add = bits(dat, pos, (1 << b) - 1) + fl[i];
					pos += b;
				}
				var d = dm[bits16(dat, pos) & dms], dsym = d >> 4;
				if (!d) err(3);
				pos += d & 15;
				var dt = fd[dsym];
				if (dsym > 3) {
					var b = fdeb[dsym];
					dt += bits16(dat, pos) & (1 << b) - 1, pos += b;
				}
				if (pos > tbts) {
					if (noSt) err(0);
					break;
				}
				if (resize) cbuf(bt + 131072);
				var end = bt + add;
				if (bt < dt) {
					var shift = dl - dt, dend = Math.min(dt, end);
					if (shift + bt < 0) err(3);
					for (; bt < dend; ++bt) buf[bt] = dict[shift + bt];
				}
				for (; bt < end; ++bt) buf[bt] = buf[bt - dt];
			}
		}
		st.l = lm, st.p = lpos, st.b = bt, st.f = final;
		if (lm) final = 1, st.m = lbt, st.d = dm, st.n = dbt;
	} while (!final);
	return bt != buf.length && noBuf ? slc(buf, 0, bt) : buf.subarray(0, bt);
};
var wbits = function(d, p, v) {
	v <<= p & 7;
	var o = p / 8 | 0;
	d[o] |= v;
	d[o + 1] |= v >> 8;
};
var wbits16 = function(d, p, v) {
	v <<= p & 7;
	var o = p / 8 | 0;
	d[o] |= v;
	d[o + 1] |= v >> 8;
	d[o + 2] |= v >> 16;
};
var hTree = function(d, mb) {
	var t = [];
	for (var i = 0; i < d.length; ++i) if (d[i]) t.push({
		s: i,
		f: d[i]
	});
	var s = t.length;
	var t2 = t.slice();
	if (!s) return {
		t: et,
		l: 0
	};
	if (s == 1) {
		var v = new u8(t[0].s + 1);
		v[t[0].s] = 1;
		return {
			t: v,
			l: 1
		};
	}
	t.sort(function(a, b) {
		return a.f - b.f;
	});
	t.push({
		s: -1,
		f: 25001
	});
	var l = t[0], r = t[1], i0 = 0, i1 = 1, i2 = 2;
	t[0] = {
		s: -1,
		f: l.f + r.f,
		l,
		r
	};
	while (i1 != s - 1) {
		l = t[t[i0].f < t[i2].f ? i0++ : i2++];
		r = t[i0 != i1 && t[i0].f < t[i2].f ? i0++ : i2++];
		t[i1++] = {
			s: -1,
			f: l.f + r.f,
			l,
			r
		};
	}
	var maxSym = t2[0].s;
	for (var i = 1; i < s; ++i) if (t2[i].s > maxSym) maxSym = t2[i].s;
	var tr = new u16(maxSym + 1);
	var mbt = ln(t[i1 - 1], tr, 0);
	if (mbt > mb) {
		var i = 0, dt = 0;
		var lft = mbt - mb, cst = 1 << lft;
		t2.sort(function(a, b) {
			return tr[b.s] - tr[a.s] || a.f - b.f;
		});
		for (; i < s; ++i) {
			var i2_1 = t2[i].s;
			if (tr[i2_1] > mb) {
				dt += cst - (1 << mbt - tr[i2_1]);
				tr[i2_1] = mb;
			} else break;
		}
		dt >>= lft;
		while (dt > 0) {
			var i2_2 = t2[i].s;
			if (tr[i2_2] < mb) dt -= 1 << mb - tr[i2_2]++ - 1;
			else ++i;
		}
		for (; i >= 0 && dt; --i) {
			var i2_3 = t2[i].s;
			if (tr[i2_3] == mb) {
				--tr[i2_3];
				++dt;
			}
		}
		mbt = mb;
	}
	return {
		t: new u8(tr),
		l: mbt
	};
};
var ln = function(n, l, d) {
	return n.s == -1 ? Math.max(ln(n.l, l, d + 1), ln(n.r, l, d + 1)) : l[n.s] = d;
};
var lc = function(c) {
	var s = c.length;
	while (s && !c[--s]);
	var cl = new u16(++s);
	var cli = 0, cln = c[0], cls = 1;
	var w = function(v) {
		cl[cli++] = v;
	};
	for (var i = 1; i <= s; ++i) if (c[i] == cln && i != s) ++cls;
	else {
		if (!cln && cls > 2) {
			for (; cls > 138; cls -= 138) w(32754);
			if (cls > 2) {
				w(cls > 10 ? cls - 11 << 5 | 28690 : cls - 3 << 5 | 12305);
				cls = 0;
			}
		} else if (cls > 3) {
			w(cln), --cls;
			for (; cls > 6; cls -= 6) w(8304);
			if (cls > 2) w(cls - 3 << 5 | 8208), cls = 0;
		}
		while (cls--) w(cln);
		cls = 1;
		cln = c[i];
	}
	return {
		c: cl.subarray(0, cli),
		n: s
	};
};
var clen = function(cf, cl) {
	var l = 0;
	for (var i = 0; i < cl.length; ++i) l += cf[i] * cl[i];
	return l;
};
var wfblk = function(out, pos, dat) {
	var s = dat.length;
	var o = shft(pos + 2);
	out[o] = s & 255;
	out[o + 1] = s >> 8;
	out[o + 2] = out[o] ^ 255;
	out[o + 3] = out[o + 1] ^ 255;
	for (var i = 0; i < s; ++i) out[o + i + 4] = dat[i];
	return (o + 4 + s) * 8;
};
var wblk = function(dat, out, final, syms, lf, df, eb, li, bs, bl, p) {
	wbits(out, p++, final);
	++lf[256];
	var _a = hTree(lf, 15), dlt = _a.t, mlb = _a.l;
	var _b = hTree(df, 15), ddt = _b.t, mdb = _b.l;
	var _c = lc(dlt), lclt = _c.c, nlc = _c.n;
	var _d = lc(ddt), lcdt = _d.c, ndc = _d.n;
	var lcfreq = new u16(19);
	for (var i = 0; i < lclt.length; ++i) ++lcfreq[lclt[i] & 31];
	for (var i = 0; i < lcdt.length; ++i) ++lcfreq[lcdt[i] & 31];
	var _e = hTree(lcfreq, 7), lct = _e.t, mlcb = _e.l;
	var nlcc = 19;
	for (; nlcc > 4 && !lct[clim[nlcc - 1]]; --nlcc);
	var flen = bl + 5 << 3;
	var ftlen = clen(lf, flt) + clen(df, fdt) + eb;
	var dtlen = clen(lf, dlt) + clen(df, ddt) + eb + 14 + 3 * nlcc + clen(lcfreq, lct) + 2 * lcfreq[16] + 3 * lcfreq[17] + 7 * lcfreq[18];
	if (bs >= 0 && flen <= ftlen && flen <= dtlen) return wfblk(out, p, dat.subarray(bs, bs + bl));
	var lm, ll, dm, dl;
	wbits(out, p, 1 + (dtlen < ftlen)), p += 2;
	if (dtlen < ftlen) {
		lm = hMap(dlt, mlb, 0), ll = dlt, dm = hMap(ddt, mdb, 0), dl = ddt;
		var llm = hMap(lct, mlcb, 0);
		wbits(out, p, nlc - 257);
		wbits(out, p + 5, ndc - 1);
		wbits(out, p + 10, nlcc - 4);
		p += 14;
		for (var i = 0; i < nlcc; ++i) wbits(out, p + 3 * i, lct[clim[i]]);
		p += 3 * nlcc;
		var lcts = [lclt, lcdt];
		for (var it = 0; it < 2; ++it) {
			var clct = lcts[it];
			for (var i = 0; i < clct.length; ++i) {
				var len = clct[i] & 31;
				wbits(out, p, llm[len]), p += lct[len];
				if (len > 15) wbits(out, p, clct[i] >> 5 & 127), p += clct[i] >> 12;
			}
		}
	} else lm = flm, ll = flt, dm = fdm, dl = fdt;
	for (var i = 0; i < li; ++i) {
		var sym = syms[i];
		if (sym > 255) {
			var len = sym >> 18 & 31;
			wbits16(out, p, lm[len + 257]), p += ll[len + 257];
			if (len > 7) wbits(out, p, sym >> 23 & 31), p += fleb[len];
			var dst = sym & 31;
			wbits16(out, p, dm[dst]), p += dl[dst];
			if (dst > 3) wbits16(out, p, sym >> 5 & 8191), p += fdeb[dst];
		} else wbits16(out, p, lm[sym]), p += ll[sym];
	}
	wbits16(out, p, lm[256]);
	return p + ll[256];
};
var deo = /* @__PURE__ */ new i32([
	65540,
	131080,
	131088,
	131104,
	262176,
	1048704,
	1048832,
	2114560,
	2117632
]);
var et = /* @__PURE__ */ new u8(0);
var dflt = function(dat, lvl, plvl, pre, post, st) {
	var s = st.z || dat.length;
	var o = new u8(pre + s + 5 * (1 + Math.ceil(s / 7e3)) + post);
	var w = o.subarray(pre, o.length - post);
	var lst = st.l;
	var pos = (st.r || 0) & 7;
	if (lvl) {
		if (pos) w[0] = st.r >> 3;
		var opt = deo[lvl - 1];
		var n = opt >> 13, c = opt & 8191;
		var msk_1 = (1 << plvl) - 1;
		var prev = st.p || new u16(32768), head = st.h || new u16(msk_1 + 1);
		var bs1_1 = Math.ceil(plvl / 3), bs2_1 = 2 * bs1_1;
		var hsh = function(i) {
			return (dat[i] ^ dat[i + 1] << bs1_1 ^ dat[i + 2] << bs2_1) & msk_1;
		};
		var syms = new i32(25e3);
		var lf = new u16(288), df = new u16(32);
		var lc_1 = 0, eb = 0, i = st.i || 0, li = 0, wi = st.w || 0, bs = 0;
		for (; i + 2 < s; ++i) {
			var hv = hsh(i);
			var imod = i & 32767, pimod = head[hv];
			prev[imod] = pimod;
			head[hv] = imod;
			if (wi <= i) {
				var rem = s - i;
				if ((lc_1 > 7e3 || li > 24576) && (rem > 423 || !lst)) {
					pos = wblk(dat, w, 0, syms, lf, df, eb, li, bs, i - bs, pos);
					li = lc_1 = eb = 0, bs = i;
					for (var j = 0; j < 286; ++j) lf[j] = 0;
					for (var j = 0; j < 30; ++j) df[j] = 0;
				}
				var l = 2, d = 0, ch_1 = c, dif = imod - pimod & 32767;
				if (rem > 2 && hv == hsh(i - dif)) {
					var maxn = Math.min(n, rem) - 1;
					var maxd = Math.min(32767, i);
					var ml = Math.min(258, rem);
					while (dif <= maxd && --ch_1 && imod != pimod) {
						if (dat[i + l] == dat[i + l - dif]) {
							var nl = 0;
							for (; nl < ml && dat[i + nl] == dat[i + nl - dif]; ++nl);
							if (nl > l) {
								l = nl, d = dif;
								if (nl > maxn) break;
								var mmd = Math.min(dif, nl - 2);
								var md = 0;
								for (var j = 0; j < mmd; ++j) {
									var ti = i - dif + j & 32767;
									var cd = ti - prev[ti] & 32767;
									if (cd > md) md = cd, pimod = ti;
								}
							}
						}
						imod = pimod, pimod = prev[imod];
						dif += imod - pimod & 32767;
					}
				}
				if (d) {
					syms[li++] = 268435456 | revfl[l] << 18 | revfd[d];
					var lin = revfl[l] & 31, din = revfd[d] & 31;
					eb += fleb[lin] + fdeb[din];
					++lf[257 + lin];
					++df[din];
					wi = i + l;
					++lc_1;
				} else {
					syms[li++] = dat[i];
					++lf[dat[i]];
				}
			}
		}
		for (i = Math.max(i, wi); i < s; ++i) {
			syms[li++] = dat[i];
			++lf[dat[i]];
		}
		pos = wblk(dat, w, lst, syms, lf, df, eb, li, bs, i - bs, pos);
		if (!lst) {
			st.r = pos & 7 | w[pos / 8 | 0] << 3;
			pos -= 7;
			st.h = head, st.p = prev, st.i = i, st.w = wi;
		}
	} else {
		for (var i = st.w || 0; i < s + lst; i += 65535) {
			var e = i + 65535;
			if (e >= s) {
				w[pos / 8 | 0] = lst;
				e = s;
			}
			pos = wfblk(w, pos + 1, dat.subarray(i, e));
		}
		st.i = s;
	}
	return slc(o, 0, pre + shft(pos) + post);
};
var crct = /* @__PURE__ */ (function() {
	var t = new Int32Array(256);
	for (var i = 0; i < 256; ++i) {
		var c = i, k = 9;
		while (--k) c = (c & 1 && -306674912) ^ c >>> 1;
		t[i] = c;
	}
	return t;
})();
var crc = function() {
	var c = -1;
	return {
		p: function(d) {
			var cr = c;
			for (var i = 0; i < d.length; ++i) cr = crct[cr & 255 ^ d[i]] ^ cr >>> 8;
			c = cr;
		},
		d: function() {
			return ~c;
		}
	};
};
var dopt = function(dat, opt, pre, post, st) {
	if (!st) {
		st = { l: 1 };
		if (opt.dictionary) {
			var dict = opt.dictionary.subarray(-32768);
			var newDat = new u8(dict.length + dat.length);
			newDat.set(dict);
			newDat.set(dat, dict.length);
			dat = newDat;
			st.w = dict.length;
		}
	}
	return dflt(dat, opt.level == null ? 6 : opt.level, opt.mem == null ? st.l ? Math.ceil(Math.max(8, Math.min(13, Math.log(dat.length))) * 1.5) : 20 : 12 + opt.mem, pre, post, st);
};
var wbytes = function(d, b, v) {
	for (; v; ++b) d[b] = v, v >>>= 8;
};
var gzh = function(c, o) {
	var fn = o.filename;
	c[0] = 31, c[1] = 139, c[2] = 8, c[8] = o.level < 2 ? 4 : o.level == 9 ? 2 : 0, c[9] = 3;
	if (o.mtime != 0) wbytes(c, 4, Math.floor(new Date(o.mtime || Date.now()) / 1e3));
	if (fn) {
		c[3] = 8;
		for (var i = 0; i <= fn.length; ++i) c[i + 10] = fn.charCodeAt(i);
	}
};
var gzs = function(d) {
	if (d[0] != 31 || d[1] != 139 || d[2] != 8) err(6, "invalid gzip data");
	var flg = d[3];
	var st = 10;
	if (flg & 4) st += (d[10] | d[11] << 8) + 2;
	for (var zs = (flg >> 3 & 1) + (flg >> 4 & 1); zs > 0; zs -= !d[st++]);
	return st + (flg & 2);
};
var gzl = function(d) {
	var l = d.length;
	return (d[l - 4] | d[l - 3] << 8 | d[l - 2] << 16 | d[l - 1] << 24) >>> 0;
};
var gzhl = function(o) {
	return 10 + (o.filename ? o.filename.length + 1 : 0);
};
/**
* Compresses data with GZIP
* @param data The data to compress
* @param opts The compression options
* @returns The gzipped version of the data
*/
function gzipSync(data, opts) {
	if (!opts) opts = {};
	var c = crc(), l = data.length;
	c.p(data);
	var d = dopt(data, opts, gzhl(opts), 8), s = d.length;
	return gzh(d, opts), wbytes(d, s - 8, c.d()), wbytes(d, s - 4, l), d;
}
/**
* Expands GZIP data
* @param data The data to decompress
* @param opts The decompression options
* @returns The decompressed version of the data
*/
function gunzipSync(data, opts) {
	var st = gzs(data);
	if (st + 8 > data.length) err(6, "invalid gzip data");
	return inflt(data.subarray(st, -8), { i: 2 }, opts && opts.out || new u8(gzl(data)), opts && opts.dictionary);
}
var td = typeof TextDecoder != "undefined" && /* @__PURE__ */ new TextDecoder();
try {
	td.decode(et, { stream: true });
} catch (e) {}
//#endregion
//#region src/result-file.ts
const MAGIC = Uint8Array.of(87, 80, 82, 0, 1, 13, 10, 26);
const MAX_COMPRESSED_BYTES = 1024 * 1024 * 1024;
const MAX_UNCOMPRESSED_BYTES = 2 * 1024 * 1024 * 1024;
const DISPOSITIONS = new Set([
	"likely_real",
	"BPB-rejects",
	"heteroduplex",
	"LDA-rejects",
	"UMI_len != 8",
	"family-size-reject"
]);
const object = (value, label) => {
	if (!value || typeof value !== "object" || Array.isArray(value)) throw new Error(`${label} must be an object.`);
	return value;
};
const array = (value, label) => {
	if (!Array.isArray(value)) throw new Error(`${label} must be an array.`);
	return value;
};
const text = (value, label) => {
	if (typeof value !== "string") throw new Error(`${label} must be text.`);
	return value;
};
const bool = (value, label) => {
	if (typeof value !== "boolean") throw new Error(`${label} must be boolean.`);
	return value;
};
const numeric = (value, label, finite = true) => {
	if (typeof value !== "number" || Number.isNaN(value) || finite && !Number.isFinite(value)) throw new Error(`${label} must be numeric.`);
	return value;
};
const count = (value, label) => {
	const result = numeric(value, label);
	if (!Number.isSafeInteger(result) || result < 0) throw new Error(`${label} must be a non-negative integer.`);
	return result;
};
const optionalText = (value, label) => {
	if (value != null) text(value, label);
};
const optionalNumber = (value, label) => {
	if (value != null) numeric(value, label);
};
const optionalBool = (value, label) => {
	if (value != null) bool(value, label);
};
function validateResult(value) {
	const bundle = object(value, "Results payload");
	if (bundle.schema !== "webporpid-results/1") throw new Error("Unsupported webPORPID result schema.");
	const provenance = object(bundle.provenance, "provenance");
	for (const key of [
		"webporpidVersion",
		"createdUtc",
		"engine",
		"inputName",
		"inputSha256",
		"configSha256",
		"deterministicSeed",
		"upstreamBranch",
		"upstreamCommit"
	]) text(provenance[key], `provenance.${key}`);
	if (count(provenance.workers, "provenance.workers") < 1) throw new Error("provenance.workers must be at least one.");
	const config = object(bundle.config, "config");
	text(config.dataset, "config.dataset");
	const parameters = object(config.parameters, "config.parameters");
	for (const key of [
		"errorRate",
		"minLength",
		"maxLength",
		"primerTolerance",
		"primerWindow",
		"primerChop",
		"maxReadsPerSample",
		"familySizeThreshold",
		"ldaThreshold",
		"contaminationClusterThreshold",
		"contaminationProportionThreshold",
		"contaminationDistanceThreshold",
		"agreementThreshold",
		"artefactFraction",
		"outlierQuantile",
		"panelThreshold",
		"functionalMatchThreshold",
		"spoolPartitions"
	]) numeric(parameters[key], `config.parameters.${key}`);
	bool(parameters.contaminationFilter, "config.parameters.contaminationFilter");
	text(parameters.deterministicSeed, "config.parameters.deterministicSeed");
	const samples = array(config.samples, "config.samples").map((entry, index) => {
		const sample = object(entry, `config.samples[${index}]`);
		const name = text(sample.name, `config.samples[${index}].name`);
		text(sample.cdnaPrimer, `config.samples[${index}].cdnaPrimer`);
		text(sample.secondStrandPrimer, `config.samples[${index}].secondStrandPrimer`);
		text(sample.panel, `config.samples[${index}].panel`);
		optionalText(sample.functionalReference, `config.samples[${index}].functionalReference`);
		return name;
	});
	if (new Set(samples).size !== samples.length) throw new Error("Result sample names must be unique.");
	const sampleSet = new Set(samples);
	const sampleIndices = new Map(samples.map((name, index) => [name, index]));
	const knownSample = (sample, label) => {
		if (!sampleSet.has(sample)) throw new Error(`${label} references an unknown sample.`);
	};
	const quality = object(bundle.quality, "quality");
	for (const key of [
		"totalReads",
		"qualityReads",
		"badReads",
		"shortReads",
		"longReads",
		"primerRejects",
		"idRejects",
		"demultiplexedReads",
		"bpbRejects",
		"malformedRecords",
		"downsampledReads"
	]) count(quality[key], `quality.${key}`);
	const perSample = array(quality.perSample, "quality.perSample");
	if (perSample.length !== samples.length) throw new Error("quality.perSample has the wrong sample count.");
	perSample.forEach((entry, index) => count(entry, `quality.perSample[${index}]`));
	const summarySamples = /* @__PURE__ */ new Set();
	array(bundle.summaries, "summaries").forEach((entry, index) => {
		const row = object(entry, `summaries[${index}]`), sample = text(row.sample, `summaries[${index}].sample`);
		if (!sampleSet.has(sample) || summarySamples.has(sample)) throw new Error("Result summaries contain an unknown or duplicate sample.");
		summarySamples.add(sample);
		for (const key of [
			"demultiplexedReads",
			"observedUmis",
			"likelyRealUmis",
			"consensusSequences"
		]) count(row[key], `summaries[${index}].${key}`);
		for (const key of [
			"contaminationPassed",
			"postprocPassed",
			"artefactCutoff"
		]) if (row[key] != null) count(row[key], `summaries[${index}].${key}`);
		if (row.selectedReads != null) count(row.selectedReads, `summaries[${index}].selectedReads`);
		if (row.downsampledReads != null) count(row.downsampledReads, `summaries[${index}].downsampledReads`);
		if (row.selectedReads != null && row.downsampledReads != null && count(row.selectedReads, `summaries[${index}].selectedReads`) + count(row.downsampledReads, `summaries[${index}].downsampledReads`) !== count(row.demultiplexedReads, `summaries[${index}].demultiplexedReads`)) throw new Error("A sample summary has inconsistent selected and subsampled read counts.");
		if (row.collapsedSequences != null) count(row.collapsedSequences, `summaries[${index}].collapsedSequences`);
		if (row.functionalPassed != null) count(row.functionalPassed, `summaries[${index}].functionalPassed`);
	});
	if (summarySamples.size !== samples.length) throw new Error("Result summaries are missing a configured sample.");
	const familyKeys = /* @__PURE__ */ new Set(), familyReadsBySample = /* @__PURE__ */ new Map();
	array(bundle.umiFamilies, "umiFamilies").forEach((entry, index) => {
		const row = object(entry, `umiFamilies[${index}]`), sample = text(row.sample, `umiFamilies[${index}].sample`), sampleIndex = count(row.sampleIndex, `umiFamilies[${index}].sampleIndex`);
		knownSample(sample, `umiFamilies[${index}]`);
		if (sampleIndices.get(sample) !== sampleIndex) throw new Error("A UMI family has an inconsistent sample index.");
		const familyKey = `${sampleIndex}\0${text(row.umi, `umiFamilies[${index}].umi`)}`;
		if (familyKeys.has(familyKey)) throw new Error("UMI family identifiers must be unique within a sample.");
		familyKeys.add(familyKey);
		const familySize = count(row.familySize, `umiFamilies[${index}].familySize`);
		familyReadsBySample.set(sample, (familyReadsBySample.get(sample) ?? 0) + familySize);
		text(row.mostLikelyParent, `umiFamilies[${index}].mostLikelyParent`);
		numeric(row.posteriorProbability, `umiFamilies[${index}].posteriorProbability`);
		numeric(row.logOffspringProbability, `umiFamilies[${index}].logOffspringProbability`, false);
		if (!DISPOSITIONS.has(text(row.disposition, `umiFamilies[${index}].disposition`))) throw new Error("A UMI family has an unknown disposition.");
		optionalNumber(row.minimumAgreement, `umiFamilies[${index}].minimumAgreement`);
	});
	array(bundle.summaries, "summaries").forEach((entry, index) => {
		const row = object(entry, `summaries[${index}]`), sample = text(row.sample, `summaries[${index}].sample`);
		if (row.selectedReads != null && count(row.selectedReads, `summaries[${index}].selectedReads`) !== (familyReadsBySample.get(sample) ?? 0)) throw new Error("A sample summary selected-read count does not match its stored family calls.");
	});
	const consensusIds = /* @__PURE__ */ new Set();
	const consensusById = /* @__PURE__ */ new Map();
	array(bundle.consensuses, "consensuses").forEach((entry, index) => {
		const row = object(entry, `consensuses[${index}]`), id = text(row.id, `consensuses[${index}].id`);
		if (consensusIds.has(id)) throw new Error("Consensus identifiers must be unique.");
		consensusIds.add(id);
		consensusById.set(id, row);
		const sample = text(row.sample, `consensuses[${index}].sample`), sampleIndex = count(row.sampleIndex, `consensuses[${index}].sampleIndex`);
		knownSample(sample, `consensuses[${index}]`);
		if (sampleIndices.get(sample) !== sampleIndex) throw new Error("A consensus has an inconsistent sample index.");
		text(row.umi, `consensuses[${index}].umi`);
		count(row.familySize, `consensuses[${index}].familySize`);
		numeric(row.minimumAgreement, `consensuses[${index}].minimumAgreement`);
		text(row.sequence, `consensuses[${index}].sequence`);
		array(row.lowAgreementSites, `consensuses[${index}].lowAgreementSites`).forEach((site, siteIndex) => {
			const low = object(site, `consensuses[${index}].lowAgreementSites[${siteIndex}]`);
			count(low.position, "low-agreement position");
			numeric(low.agreement, "low-agreement value");
			text(low.modalReadBase, "low-agreement modal base");
			count(low.modalRunLength, "low-agreement run length");
		});
	});
	array(bundle.contamination, "contamination").forEach((entry, index) => {
		const row = object(entry, `contamination[${index}]`), sample = text(row.sample, "contamination sample"), sequenceId = text(row.sequenceId, "contamination sequence ID");
		knownSample(sample, `contamination[${index}]`);
		if (consensusById.get(sequenceId)?.sample !== sample) throw new Error("A contamination call references an unknown consensus or sample.");
		text(row.nearestNonselfVariant, "nearest non-self variant");
		numeric(row.nearestNonselfDistance, "nearest non-self distance");
		bool(row.flagged, "contamination flagged");
		bool(row.discarded, "contamination discarded");
		bool(row.suspectOnly, "contamination suspectOnly");
	});
	if (bundle.contaminationReferences != null) array(bundle.contaminationReferences, "contaminationReferences").forEach((entry, index) => {
		const record = object(entry, `contaminationReferences[${index}]`);
		text(record.name, `contaminationReferences[${index}].name`);
		text(record.sequence, `contaminationReferences[${index}].sequence`);
	});
	if (bundle.downstreamResources != null) {
		const resources = object(bundle.downstreamResources, "downstreamResources"), resourceSamples = /* @__PURE__ */ new Set();
		array(resources.samples, "downstreamResources.samples").forEach((rawSample, index) => {
			const resource = object(rawSample, `downstreamResources.samples[${index}]`), name = text(resource.name, "downstream resource sample");
			knownSample(name, `downstreamResources.samples[${index}]`);
			if (resourceSamples.has(name)) throw new Error("Downstream resource sample names must be unique.");
			resourceSamples.add(name);
			array(resource.panelSequences, "downstream panel sequences").forEach((rawRecord) => {
				const record = object(rawRecord, "downstream panel record");
				text(record.name, "downstream panel name");
				text(record.sequence, "downstream panel sequence");
			});
			if (resource.functionalReferenceSequence != null) {
				const record = object(resource.functionalReferenceSequence, "downstream functional reference");
				text(record.name, "downstream functional reference name");
				text(record.sequence, "downstream functional reference sequence");
			}
		});
		if (resourceSamples.size !== samples.length) throw new Error("Downstream resources are missing a configured sample.");
	}
	const rawStageStatuses = bundle.optionalStages == null ? void 0 : object(bundle.optionalStages, "optionalStages");
	const postprocessingComplete = rawStageStatuses == null || object(rawStageStatuses.postprocessing, "optionalStages.postprocessing").state === "completed";
	const recordIds = /* @__PURE__ */ new Set(), recordMetadata = /* @__PURE__ */ new Map();
	array(bundle.records, "records").forEach((entry, index) => {
		const row = object(entry, `records[${index}]`), id = text(row.id, `records[${index}].id`);
		if (recordIds.has(id)) throw new Error("Post-processing identifiers must be unique.");
		recordIds.add(id);
		const sample = text(row.sample, "postproc sample"), source = consensusById.get(id);
		knownSample(sample, `records[${index}]`);
		if (!source || source.sample !== sample) throw new Error("A post-processing record references an unknown consensus or sample.");
		const umi = text(row.umi, "postproc UMI"), familySize = count(row.familySize, "postproc family size"), minimumAgreement = numeric(row.minimumAgreement, "postproc agreement");
		if (umi !== source.umi || familySize !== source.familySize || minimumAgreement !== source.minimumAgreement) throw new Error("A post-processing record has inconsistent consensus metadata.");
		if (text(row.consensusNt, "postproc consensus") !== source.sequence) throw new Error("A post-processing record has inconsistent consensus sequence data.");
		const alignedNt = row.alignedNt == null ? void 0 : text(row.alignedNt, "postproc aligned sequence");
		optionalText(row.trimmedNt, "postproc trimmed nucleotide");
		optionalText(row.trimmedAa, "postproc trimmed protein");
		numeric(row.panelScore, "postproc panel score");
		for (const key of [
			"artefactPass",
			"agreementPass",
			"contaminationPass",
			"panelPass"
		]) bool(row[key], `postproc ${key}`);
		optionalBool(row.functionalPass, "postproc functionalPass");
		array(row.rejectionReasons, "postproc rejectionReasons").forEach((reason) => text(reason, "postproc rejection reason"));
		if (row.apobec != null) {
			const model = object(row.apobec, "postproc APOBEC");
			for (const key of [
				"posteriorMeanGaMultiplier",
				"posteriorGaInflated",
				"posteriorMeanMutationRate",
				"gaMutations",
				"totalMutations"
			]) numeric(model[key], `APOBEC ${key}`);
		}
		recordMetadata.set(id, {
			sample,
			alignedNt,
			minimumAgreement
		});
	});
	if (postprocessingComplete && (recordIds.size !== consensusIds.size || [...consensusIds].some((id) => !recordIds.has(id)))) throw new Error("Consensus and post-processing records are inconsistent.");
	if (!postprocessingComplete && recordIds.size) throw new Error("An uncomputed post-processing stage cannot contain partial post-processing records.");
	for (const [label, entries] of [["alignments", object(bundle.alignments, "alignments")], ["trees", object(bundle.trees, "trees")]]) for (const [name, contents] of Object.entries(entries)) {
		text(name, `${label} name`);
		text(contents, `${label}.${name}`);
		const sample = name.split("/", 1)[0];
		knownSample(sample, `${label}.${name}`);
	}
	if (bundle.referenceAlignments != null) for (const [name, contents] of Object.entries(object(bundle.referenceAlignments, "referenceAlignments"))) {
		const sample = name.split("/", 1)[0];
		knownSample(sample, `referenceAlignments.${name}`);
		const reference = inspectAlignment(text(contents, `referenceAlignments.${name}`), 1);
		const nucleotide = object(bundle.alignments, "alignments")[name];
		if (nucleotide != null && inspectAlignment(text(nucleotide, `alignments.${name}`), 1).columns !== reference.columns) throw new Error("A stored reference row does not match its nucleotide alignment width.");
	}
	if (bundle.collapseGroups != null) for (const [sample, rawGroups] of Object.entries(object(bundle.collapseGroups, "collapseGroups"))) {
		knownSample(sample, `collapseGroups.${sample}`);
		const representatives = /* @__PURE__ */ new Set(), membersSeen = /* @__PURE__ */ new Set();
		const collapsed = inspectAlignment(text(object(bundle.alignments, "alignments")[`${sample}/nucleotide`], `alignments.${sample}/nucleotide`), 1);
		const uncollapsed = inspectAlignment(text(object(bundle.alignments, "alignments")[`${sample}/uncollapsed-nucleotide`], `alignments.${sample}/uncollapsed-nucleotide`), 1);
		const collapsedByName = new Map(collapsed.records.map((record) => [record.name, record.sequence.replaceAll("-", "")]));
		const uncollapsedByName = new Map(uncollapsed.records.map((record) => [record.name, record.sequence.replaceAll("-", "")]));
		array(rawGroups, `collapseGroups.${sample}`).forEach((rawGroup, index) => {
			const group = object(rawGroup, `collapseGroups.${sample}[${index}]`);
			if (text(group.sample, "collapse group sample") !== sample) throw new Error("A collapse group has an inconsistent sample.");
			const representative = text(group.representativeId, "collapse representative");
			if (representatives.has(representative)) throw new Error("Collapse representative identifiers must be unique.");
			representatives.add(representative);
			const members = array(group.memberIds, "collapse members").map((entry) => text(entry, "collapse member"));
			if (!members.includes(representative)) throw new Error("A collapse group must include its representative.");
			if (count(group.familyCount, "collapse family count") !== members.length) throw new Error("Collapse counts must count UMI families.");
			const representativeSequence = collapsedByName.get(representative);
			if (representativeSequence == null) throw new Error("A collapse representative is missing from the collapsed alignment.");
			const agreements = [];
			for (const member of members) {
				if (membersSeen.has(member)) throw new Error("A retained UMI family occurs in more than one collapse group.");
				membersSeen.add(member);
				if (uncollapsedByName.get(member) !== representativeSequence) throw new Error("A collapse group contains different nucleotide haplotypes.");
				const metadata = recordMetadata.get(member);
				if (!metadata || metadata.sample !== sample || metadata.alignedNt == null) throw new Error("A collapse member is not a retained UMI-family consensus.");
				agreements.push(metadata.minimumAgreement);
			}
			if (group.minimumAgreement != null && numeric(group.minimumAgreement, "legacy collapse minimum agreement") !== Math.min(...agreements)) throw new Error("A legacy collapse group has inconsistent family-agreement metadata.");
		});
		if (representatives.size !== collapsed.records.length || membersSeen.size !== uncollapsed.records.length) throw new Error("Collapse membership does not cover the stored nucleotide alignments.");
		const summary = array(bundle.summaries, "summaries").map((entry) => object(entry, "summary")).find((entry) => entry.sample === sample);
		if (summary?.collapsedSequences != null && count(summary.collapsedSequences, "summary collapsed count") !== representatives.size) throw new Error("A summary has an inconsistent collapsed haplotype count.");
	}
	if (bundle.inputMappings != null) array(bundle.inputMappings, "inputMappings").forEach((rawMapping, index) => {
		const mapping = object(rawMapping, `inputMappings[${index}]`);
		text(mapping.slot, "input slot");
		if (![
			"reads",
			"configuration",
			"panel",
			"functional-reference",
			"contamination-panel"
		].includes(text(mapping.role, "input role"))) throw new Error("An input mapping has an unknown role.");
		optionalText(mapping.expectedName, "expected filename");
		text(mapping.uploadedName, "uploaded filename");
		count(mapping.uploadedSize, "uploaded size");
	});
	if (bundle.runOptions != null) {
		const options = object(bundle.runOptions, "runOptions");
		bool(options.deferPhylogeny, "runOptions.deferPhylogeny");
		for (const key of [
			"deferContamination",
			"deferPostprocessing",
			"deferCollapse"
		]) if (options[key] != null) bool(options[key], `runOptions.${key}`);
		if (options.spoolStorage != null && !["automatic", "external-directory"].includes(text(options.spoolStorage, "runOptions.spoolStorage"))) throw new Error("runOptions.spoolStorage is not recognized.");
	}
	if (rawStageStatuses != null) {
		let prerequisiteIncomplete = false;
		for (const stage of [
			"contamination",
			"postprocessing",
			"collapse",
			"tree"
		]) {
			const status = object(rawStageStatuses[stage], `optionalStages.${stage}`), state = text(status.state, `optionalStages.${stage}.state`);
			if (![
				"completed",
				"deferred",
				"skipped"
			].includes(state)) throw new Error(`optionalStages.${stage}.state is not recognized.`);
			text(status.detail, `optionalStages.${stage}.detail`);
			text(status.updatedUtc, `optionalStages.${stage}.updatedUtc`);
			if (prerequisiteIncomplete && state === "completed") throw new Error(`optionalStages.${stage} cannot be completed before its prerequisite.`);
			if (state !== "completed") prerequisiteIncomplete = true;
		}
	}
	if (bundle.alignmentEdits != null) for (const [name, rawEdit] of Object.entries(object(bundle.alignmentEdits, "alignmentEdits"))) {
		const sample = name.split("/", 1)[0];
		knownSample(sample, `alignmentEdits.${name}`);
		if (name !== `${sample}/nucleotide` && name !== `${sample}/uncollapsed-nucleotide` && name !== `${sample}/functional-nucleotide`) throw new Error("Edited alignment keys must identify a stored nucleotide view.");
		const edit = object(rawEdit, `alignmentEdits.${name}`), fasta = text(edit.fasta, `alignmentEdits.${name}.fasta`);
		if (count(edit.frameOffset, `alignmentEdits.${name}.frameOffset`) > 2) throw new Error("Edited alignment frame offsets must be 0, 1, or 2.");
		const baselineFingerprint = text(edit.baselineFingerprint, `alignmentEdits.${name}.baselineFingerprint`);
		const editedFingerprint = text(edit.editedFingerprint, `alignmentEdits.${name}.editedFingerprint`);
		text(edit.source, `alignmentEdits.${name}.source`);
		text(edit.savedUtc, `alignmentEdits.${name}.savedUtc`);
		optionalText(edit.treeNewick, `alignmentEdits.${name}.treeNewick`);
		optionalText(edit.treeFingerprint, `alignmentEdits.${name}.treeFingerprint`);
		optionalBool(edit.treeStale, `alignmentEdits.${name}.treeStale`);
		if (edit.warnings != null) array(edit.warnings, `alignmentEdits.${name}.warnings`).forEach((warning) => text(warning, "alignment edit warning"));
		if (edit.changes != null) {
			const changes = object(edit.changes, `alignmentEdits.${name}.changes`);
			for (const field of [
				"rowsBefore",
				"rowsAfter",
				"columnsBefore",
				"columnsAfter",
				"removedNucleotides",
				"insertedNucleotides",
				"substitutedNucleotides"
			]) count(changes[field], `alignmentEdits.${name}.changes.${field}`);
			bool(changes.rowOrderChanged, `alignmentEdits.${name}.changes.rowOrderChanged`);
			for (const field of [
				"rowOrderBefore",
				"rowOrderAfter",
				"removedRows",
				"addedRows",
				"changedRows"
			]) array(changes[field], `alignmentEdits.${name}.changes.${field}`).forEach((value) => text(value, `alignmentEdits.${name}.changes.${field} entry`));
			array(changes.rowChanges, `alignmentEdits.${name}.changes.rowChanges`).forEach((rawRow, index) => {
				const row = object(rawRow, `alignmentEdits.${name}.changes.rowChanges[${index}]`);
				text(row.name, "alignment row change name");
				for (const field of [
					"removedNucleotides",
					"insertedNucleotides",
					"substitutedNucleotides"
				]) count(row[field], `alignment row change ${field}`);
				bool(row.gapPlacementChanged, "alignment row gap-placement flag");
			});
		}
		const storedAlignments = object(bundle.alignments, "alignments");
		const original = text(storedAlignments[name] ?? (name === `${sample}/uncollapsed-nucleotide` ? storedAlignments[`${sample}/nucleotide`] : void 0), `alignments.${name}`);
		if (inspectAlignment(original, 1).fingerprint !== baselineFingerprint) throw new Error(`alignmentEdits.${name} has an inconsistent baseline fingerprint.`);
		if (inspectAlignment(fasta, 1).fingerprint !== editedFingerprint) throw new Error(`alignmentEdits.${name} has an inconsistent fingerprint.`);
		validateCorrectedAlignment(original, fasta);
		if (edit.changes != null) {
			const expected = summarizeAlignmentChanges(original, fasta), stored = edit.changes;
			for (const field of [
				"rowsBefore",
				"rowsAfter",
				"columnsBefore",
				"columnsAfter",
				"removedNucleotides",
				"insertedNucleotides",
				"substitutedNucleotides"
			]) if (stored[field] !== expected[field]) throw new Error(`alignmentEdits.${name}.changes.${field} is inconsistent with the stored alignments.`);
			if (stored.rowOrderChanged !== expected.rowOrderChanged) throw new Error(`alignmentEdits.${name}.changes.rowOrderChanged is inconsistent with the stored alignments.`);
			for (const field of [
				"rowOrderBefore",
				"rowOrderAfter",
				"removedRows",
				"addedRows",
				"changedRows"
			]) if (JSON.stringify(stored[field]) !== JSON.stringify(expected[field])) throw new Error(`alignmentEdits.${name}.changes.${field} is inconsistent with the stored alignments.`);
			if (JSON.stringify(stored.rowChanges) !== JSON.stringify(expected.rowChanges)) throw new Error(`alignmentEdits.${name}.changes.rowChanges is inconsistent with the stored alignments.`);
		}
	}
	if (bundle.alignmentEditHistory != null) array(bundle.alignmentEditHistory, "alignmentEditHistory").forEach((rawEntry, index) => {
		const entry = object(rawEntry, `alignmentEditHistory[${index}]`), key = text(entry.alignmentKey, `alignmentEditHistory[${index}].alignmentKey`);
		const sample = key.split("/", 1)[0];
		knownSample(sample, `alignmentEditHistory[${index}]`);
		if (key !== `${sample}/nucleotide` && key !== `${sample}/uncollapsed-nucleotide` && key !== `${sample}/functional-nucleotide`) throw new Error("Alignment audit keys must identify a nucleotide view.");
		if (![
			"alignment-edit",
			"frame-change",
			"tree-recalculation",
			"edit-reset"
		].includes(text(entry.action, `alignmentEditHistory[${index}].action`))) throw new Error("Alignment audit action is not recognized.");
		text(entry.timestamp, `alignmentEditHistory[${index}].timestamp`);
		text(entry.source, `alignmentEditHistory[${index}].source`);
		array(entry.details, `alignmentEditHistory[${index}].details`).forEach((value) => text(value, "alignment audit detail"));
		optionalText(entry.beforeFingerprint, `alignmentEditHistory[${index}].beforeFingerprint`);
		optionalText(entry.afterFingerprint, `alignmentEditHistory[${index}].afterFingerprint`);
	});
	if (bundle.timings != null) array(bundle.timings, "timings").forEach((entry, index) => {
		const timing = object(entry, `timings[${index}]`);
		const stage = text(timing.stage, `timings[${index}].stage`);
		if (![
			"setup",
			"preprocessing",
			"umi",
			"consensus",
			"contamination",
			"postprocessing",
			"collapse",
			"tree",
			"analysis-total"
		].includes(stage)) throw new Error(`timings[${index}] has an unknown stage.`);
		if (numeric(timing.seconds, `timings[${index}].seconds`) < 0) throw new Error(`timings[${index}].seconds must be non-negative.`);
		if (timing.workItems != null) count(timing.workItems, `timings[${index}].workItems`);
	});
	array(bundle.log, "log").forEach((entry, index) => text(entry, `log[${index}]`));
}
function encodeResult(bundle) {
	validateResult(bundle);
	const body = gzipSync(encode(bundle), { level: 9 }), output = new Uint8Array(MAGIC.byteLength + body.byteLength);
	output.set(MAGIC);
	output.set(body, MAGIC.byteLength);
	return output;
}
function decodeResult(bytes) {
	if (bytes.byteLength < MAGIC.byteLength || MAGIC.some((value, index) => bytes[index] !== value)) throw new Error("This is not a webPORPID results file.");
	const compressed = bytes.subarray(MAGIC.byteLength);
	if (compressed.byteLength > MAX_COMPRESSED_BYTES) throw new Error("The webPORPID results file is too large to load safely.");
	if (compressed.byteLength < 18) throw new Error("The webPORPID results payload is truncated.");
	if (new DataView(compressed.buffer, compressed.byteOffset + compressed.byteLength - 4, 4).getUint32(0, true) > MAX_UNCOMPRESSED_BYTES) throw new Error("The uncompressed webPORPID results payload is too large to load safely.");
	let unpacked;
	try {
		unpacked = gunzipSync(compressed);
	} catch {
		throw new Error("The webPORPID results payload is corrupt or truncated.");
	}
	if (unpacked.byteLength > MAX_UNCOMPRESSED_BYTES) throw new Error("The uncompressed webPORPID results payload is too large to load safely.");
	let value;
	try {
		value = decode(unpacked);
	} catch {
		throw new Error("The webPORPID results payload contains invalid MessagePack data.");
	}
	validateResult(value);
	return value;
}
const quote = (value) => {
	const valueText = value == null ? "" : String(value);
	return /[",\r\n]/.test(valueText) ? `"${valueText.replaceAll("\"", "\"\"")}"` : valueText;
};
const csv = (headers, rows) => [headers, ...rows].map((row) => row.map(quote).join(",")).join("\n") + "\n";
const fasta = (rows) => rows.map((row) => `>${row.id}\n${row.sequence.match(/.{1,80}/g)?.join("\n") ?? ""}`).join("\n") + (rows.length ? "\n" : "");
const passed = (record) => record.artefactPass && record.agreementPass && record.contaminationPass && record.panelPass;
function alignmentSample(bundle, sample) {
	if (sample) return sample;
	if (bundle.config.samples.length === 1) return bundle.config.samples[0].name;
	throw new Error("Choose a sample when exporting a sample-specific alignment or tree.");
}
function exportComponent(bundle, kind, sample) {
	const consensuses = bundle.consensuses.filter((record) => !sample || record.sample === sample);
	const records = bundle.records.filter((record) => !sample || record.sample === sample);
	switch (kind) {
		case "consensus-fasta": return {
			extension: "consensus.fasta",
			mime: "text/x-fasta",
			text: fasta(consensuses.map((record) => ({
				id: record.id,
				sequence: record.sequence
			})))
		};
		case "passed-consensus-fasta": return {
			extension: "passed-consensus.fasta",
			mime: "text/x-fasta",
			text: fasta(records.filter(passed).map((record) => ({
				id: record.id,
				sequence: record.consensusNt
			})))
		};
		case "rejected-consensus-fasta": return {
			extension: "rejected-consensus.fasta",
			mime: "text/x-fasta",
			text: fasta(records.filter((record) => !passed(record)).map((record) => ({
				id: record.id,
				sequence: record.consensusNt
			})))
		};
		case "trimmed-nt-fasta": return {
			extension: "trimmed-nt.fasta",
			mime: "text/x-fasta",
			text: fasta(records.filter((record) => record.functionalPass && record.trimmedNt).map((record) => ({
				id: record.id,
				sequence: record.trimmedNt
			})))
		};
		case "trimmed-aa-fasta": return {
			extension: "trimmed-aa.fasta",
			mime: "text/x-fasta",
			text: fasta(records.filter((record) => record.functionalPass && record.trimmedAa).map((record) => ({
				id: record.id,
				sequence: record.trimmedAa
			})))
		};
		case "family-csv": return {
			extension: "families.csv",
			mime: "text/csv",
			text: csv([
				"sample",
				"UMI",
				"fs",
				"tags",
				"posterior_probability",
				"log_offspring_probability",
				"minag"
			], bundle.umiFamilies.filter((row) => !sample || row.sample === sample).map((row) => [
				row.sample,
				row.umi,
				row.familySize,
				row.disposition,
				row.posteriorProbability,
				row.logOffspringProbability,
				row.minimumAgreement
			]))
		};
		case "low-agreement-csv": return {
			extension: "low-agreement.csv",
			mime: "text/csv",
			text: csv([
				"sample",
				"sequence_id",
				"UMI",
				"position_from_3prime",
				"agreement",
				"modal_read_base",
				"modal_run_length"
			], consensuses.flatMap((record) => record.lowAgreementSites.map((site) => [
				record.sample,
				record.id,
				record.umi,
				site.position,
				site.agreement,
				site.modalReadBase,
				site.modalRunLength
			])))
		};
		case "contamination-csv": return {
			extension: "contamination.csv",
			mime: "text/csv",
			text: csv([
				"sample",
				"sequence_name",
				"nearest_nonself_variant",
				"nearest_nonself_distance",
				"flagged",
				"discarded",
				"suspect_only"
			], deduplicateContaminationCalls(bundle.contamination).filter((row) => !sample || row.sample === sample).map((row) => [
				row.sample,
				row.sequenceId,
				row.nearestNonselfVariant,
				row.nearestNonselfDistance,
				row.flagged,
				row.discarded,
				row.suspectOnly
			]))
		};
		case "postproc-csv": return {
			extension: "postproc.csv",
			mime: "text/csv",
			text: csv([
				"sample",
				"id",
				"UMI",
				"fs",
				"minag",
				"panel_score",
				"artefact_pass",
				"agreement_pass",
				"contamination_pass",
				"panel_pass",
				"functional_pass",
				"rejection_reasons"
			], records.map((row) => [
				row.sample,
				row.id,
				row.umi,
				row.familySize,
				row.minimumAgreement,
				row.panelScore,
				row.artefactPass,
				row.agreementPass,
				row.contaminationPass,
				row.panelPass,
				row.functionalPass,
				row.rejectionReasons.join(";")
			]))
		};
		case "apobec-csv": return {
			extension: "apobec.csv",
			mime: "text/csv",
			text: csv([
				"sample",
				"id",
				"posterior_mean_GA_multiplier",
				"posterior_probability_GA_inflated",
				"posterior_mean_mutation_rate",
				"GA_mutations",
				"total_mutations"
			], records.filter((row) => row.apobec).map((row) => [
				row.sample,
				row.id,
				row.apobec.posteriorMeanGaMultiplier,
				row.apobec.posteriorGaInflated,
				row.apobec.posteriorMeanMutationRate,
				row.apobec.gaMutations,
				row.apobec.totalMutations
			]))
		};
		case "collapse-csv": {
			const selected = alignmentSample(bundle, sample);
			return {
				extension: "collapsed-families.csv",
				mime: "text/csv",
				text: csv([
					"sample",
					"representative_id",
					"family_count",
					"member_ids"
				], (bundle.collapseGroups?.[selected] ?? []).map((group) => [
					selected,
					group.representativeId,
					group.familyCount,
					group.memberIds.join(";")
				]))
			};
		}
		case "nucleotide-alignment": {
			const key = `${alignmentSample(bundle, sample)}/nucleotide`;
			return {
				extension: "nucleotide-alignment.fasta",
				mime: "text/x-fasta",
				text: bundle.alignmentEdits?.[key]?.fasta ?? bundle.alignments[key] ?? ""
			};
		}
		case "protein-alignment": {
			const key = `${alignmentSample(bundle, sample)}/nucleotide`, edit = bundle.alignmentEdits?.[key];
			const nucleotide = edit?.fasta ?? bundle.alignments[key];
			return {
				extension: "protein-alignment.fasta",
				mime: "text/x-fasta",
				text: nucleotide ? translateAlignmentFasta(nucleotide, edit?.frameOffset ?? 0) : ""
			};
		}
		case "newick": {
			const key = `${alignmentSample(bundle, sample)}/nucleotide`, edit = bundle.alignmentEdits?.[key];
			return {
				extension: "tree.newick",
				mime: "text/plain",
				text: edit ? edit.treeNewick ?? "" : bundle.trees[key] ?? ""
			};
		}
		case "uncollapsed-nucleotide-alignment": {
			const selected = alignmentSample(bundle, sample), key = `${selected}/uncollapsed-nucleotide`;
			return {
				extension: "uncollapsed-nucleotide-alignment.fasta",
				mime: "text/x-fasta",
				text: (bundle.alignmentEdits?.[key])?.fasta ?? bundle.alignments[key] ?? bundle.alignments[`${selected}/nucleotide`] ?? ""
			};
		}
		case "uncollapsed-protein-alignment": {
			const selected = alignmentSample(bundle, sample), key = `${selected}/uncollapsed-nucleotide`, edit = bundle.alignmentEdits?.[key];
			const nucleotide = edit?.fasta ?? bundle.alignments[key] ?? bundle.alignments[`${selected}/nucleotide`];
			return {
				extension: "uncollapsed-protein-alignment.fasta",
				mime: "text/x-fasta",
				text: nucleotide ? translateAlignmentFasta(nucleotide, edit?.frameOffset ?? 0) : ""
			};
		}
		case "uncollapsed-newick": {
			const key = `${alignmentSample(bundle, sample)}/uncollapsed-nucleotide`;
			return {
				extension: "uncollapsed-tree.newick",
				mime: "text/plain",
				text: (bundle.alignmentEdits?.[key])?.treeNewick ?? bundle.trees[key] ?? ""
			};
		}
		case "functional-nucleotide-alignment": {
			const key = `${alignmentSample(bundle, sample)}/functional-nucleotide`;
			return {
				extension: "functional-nucleotide-alignment.fasta",
				mime: "text/x-fasta",
				text: (bundle.alignmentEdits?.[key])?.fasta ?? bundle.alignments[key] ?? ""
			};
		}
		case "functional-protein-alignment": {
			const selected = alignmentSample(bundle, sample), key = `${selected}/functional-nucleotide`, edit = bundle.alignmentEdits?.[key];
			const nucleotide = edit?.fasta ?? bundle.alignments[key];
			return {
				extension: "functional-protein-alignment.fasta",
				mime: "text/x-fasta",
				text: nucleotide ? translateAlignmentFasta(nucleotide, edit?.frameOffset ?? 0) : bundle.alignments[`${selected}/functional-protein`] ?? ""
			};
		}
		case "functional-newick": {
			const key = `${alignmentSample(bundle, sample)}/functional-nucleotide`;
			return {
				extension: "functional-tree.newick",
				mime: "text/plain",
				text: (bundle.alignmentEdits?.[key])?.treeNewick ?? bundle.trees[key] ?? ""
			};
		}
		case "log": return {
			extension: "log.txt",
			mime: "text/plain",
			text: bundle.log.join("\n") + "\n"
		};
	}
}
function safeDatasetName(value) {
	return value.replace(/[^A-Za-z0-9_.-]+/g, "_") || "webporpid";
}
const MAX_SPOOL_RECORD_BYTES = 256 * 1024 * 1024;
function parseSpoolRecordHeader(bytes) {
	if (bytes.byteLength < 24) throw new Error("A temporary spool record has a truncated header.");
	const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
	const bodyLength = view.getUint32(0, true), sample = view.getUint16(4, true);
	const umiLength = view.getUint16(6, true), nameLength = view.getUint32(8, true), sequenceLength = view.getUint32(12, true);
	const recordLength = bodyLength + 4;
	const expectedBody = 20 + umiLength + nameLength + sequenceLength * 2;
	if (bodyLength < 20 || recordLength > MAX_SPOOL_RECORD_BYTES || bodyLength !== expectedBody || !Number.isSafeInteger(recordLength)) throw new Error("A temporary spool record has inconsistent lengths.");
	return {
		recordLength,
		sample,
		samplingHash: view.getBigUint64(16, true)
	};
}
function selectedSpoolRecord(header, cutoffs) {
	return header.sample >= cutoffs.length || header.samplingHash <= cutoffs[header.sample];
}
function concatenateSpoolRecords(records) {
	const length = records.reduce((sum, record) => sum + record.byteLength, 0);
	const output = new Uint8Array(length);
	let offset = 0;
	for (const record of records) {
		output.set(record, offset);
		offset += record.byteLength;
	}
	return output;
}
//#endregion
//#region cli-src/porpid-cli.mjs
const VERSION = "0.3.6";
const UPSTREAM_COMMIT = "201af7942029cfb7974880e41674be9f0ddfaf3b";
const CLI_DIRECTORY = dirname(new URL(import.meta.url).pathname);
function defaultCliAssets() {
	const directory = join(CLI_DIRECTORY, "assets");
	return {
		wasmPath: join(directory, "webporpid.wasm"),
		msaPath: join(directory, "alivibe-msa.wasm"),
		fastTreeJavascriptPath: join(directory, "fasttree.cjs"),
		fastTreeWasmPath: join(directory, "fasttree.wasm"),
		msaWorkerPath: join(CLI_DIRECTORY, "porpid-msa-worker.mjs"),
		fastTreeWorkerPath: join(CLI_DIRECTORY, "porpid-fasttree-worker.mjs")
	};
}
function usage() {
	return `porpid-cli ${VERSION}\n\nRun the complete nanopore/PacBio pipeline:\n  porpid-cli run reads.fastq.gz --config config.yaml --output results.webporpid [--workers N] [--defer-phylogeny]\n\nInspect or export a saved analysis:\n  porpid-cli inspect results.webporpid\n  porpid-cli export results.webporpid --component consensus-fasta [--sample NAME] --output consensus.fasta\n\nWorkers default to all logical CPUs (${availableParallelism()}). Temporary read partitions are streamed to disk and removed after consensus.\nComponents: consensus-fasta, passed-consensus-fasta, rejected-consensus-fasta, trimmed-nt-fasta, trimmed-aa-fasta,\n            family-csv, low-agreement-csv, contamination-csv, postproc-csv, apobec-csv, collapse-csv,\n            nucleotide-alignment, protein-alignment, newick, uncollapsed-nucleotide-alignment,\n            uncollapsed-protein-alignment, uncollapsed-newick, functional-nucleotide-alignment,\n            functional-protein-alignment, functional-newick, log`;
}
function option(args, name) {
	const index = args.indexOf(name);
	if (index >= 0) return args[index + 1];
	return args.find((value) => value.startsWith(`${name}=`))?.slice(name.length + 1);
}
function integer(value, label, fallback) {
	if (value == null) return fallback;
	if (!/^\d+$/.test(value) || Number(value) < 1 || !Number.isSafeInteger(Number(value))) throw new Error(`${label} requires an integer of at least one.`);
	return Number(value);
}
function status(message) {
	process.stderr.write(`[webPORPID] ${message}\n`);
}
const now = () => (/* @__PURE__ */ new Date()).toISOString();
var WorkerClient = class {
	constructor(worker, webWorker) {
		this.worker = worker;
		this.webWorker = webWorker;
		this.pending = /* @__PURE__ */ new Map();
		this.nextId = 1;
		this.tail = Promise.resolve();
		const receive = (message) => {
			const pending = this.pending.get(message.id);
			if (!pending) return;
			this.pending.delete(message.id);
			message.error ? pending.reject(new Error(message.error)) : pending.resolve(message.result);
		};
		const fail = (cause) => {
			const error = cause instanceof Error ? cause : new Error(cause?.message ?? String(cause));
			for (const pending of this.pending.values()) pending.reject(error);
			this.pending.clear();
		};
		if (webWorker) {
			worker.addEventListener("message", (event) => receive(event.data));
			worker.addEventListener("error", fail);
		} else {
			worker.on("message", receive);
			worker.on("error", fail);
		}
	}
	raw(message, transfer = []) {
		return new Promise((resolvePromise, reject) => {
			const id = this.nextId++;
			this.pending.set(id, {
				resolve: resolvePromise,
				reject
			});
			if (this.webWorker) this.worker.postMessage({
				id,
				...message
			}, { transfer });
			else this.worker.postMessage({
				id,
				...message
			}, transfer);
		});
	}
	call(message, transfer = []) {
		const start = () => this.raw(message, transfer), result = this.tail.then(start, start);
		this.tail = result.catch(() => {});
		return result;
	}
	terminate() {
		return Promise.resolve(this.worker.terminate());
	}
};
var WorkerPool = class WorkerPool {
	constructor(clients) {
		this.clients = clients;
		this.cursor = 0;
	}
	static async create(size, wasmPath, compiledConfig) {
		const clients = [];
		for (let index = 0; index < size; index++) {
			const webWorker = Boolean(process.versions.bun && globalThis.Worker);
			const worker = webWorker ? new globalThis.Worker(new URL("./porpid-worker.mjs", import.meta.url)) : new Worker$1(new URL("./porpid-worker.mjs", import.meta.url));
			clients.push(new WorkerClient(worker, webWorker));
		}
		await Promise.all(clients.map((client) => {
			const copy = compiledConfig.slice().buffer;
			return client.call({
				type: "init",
				wasmPath,
				config: copy
			}, [copy]);
		}));
		return new WorkerPool(clients);
	}
	any(message, transfer = []) {
		return this.clients[this.cursor++ % this.clients.length].call(message, transfer);
	}
	at(index, message, transfer = []) {
		return this.clients[index].call(message, transfer);
	}
	async close() {
		await Promise.all(this.clients.map((client) => client.terminate()));
	}
};
var DiskPartitions = class DiskPartitions {
	static async create(count) {
		let base = process.env.WEBPORPID_TMPDIR || tmpdir();
		try {
			await mkdir(base, { recursive: true });
		} catch {
			base = process.cwd();
		}
		const directory = await mkdtemp(join(base, "webporpid-")), handles = [], paths = [];
		for (let index = 0; index < count; index++) {
			const path = join(directory, `partition-${index}.bin`);
			paths.push(path);
			handles.push(await open(path, "w+"));
		}
		return new DiskPartitions(directory, handles, paths, Array.from({ length: count }, () => Promise.resolve()), Array(count).fill(0));
	}
	constructor(directory, handles, paths, chains, lengths) {
		Object.assign(this, {
			directory,
			handles,
			paths,
			chains,
			lengths
		});
		this.closed = false;
	}
	async appendFrames(bytes) {
		const grouped = /* @__PURE__ */ new Map();
		const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
		let offset = 0;
		while (offset < bytes.byteLength) {
			if (offset + 5 > bytes.byteLength) throw new Error("Truncated spool routing frame.");
			const partition = bytes[offset++], length = view.getUint32(offset, true);
			offset += 4;
			if (partition >= this.handles.length || offset + length > bytes.byteLength) throw new Error("Invalid spool routing frame.");
			const parts = grouped.get(partition) ?? [];
			parts.push(bytes.slice(offset, offset + length));
			grouped.set(partition, parts);
			offset += length;
		}
		await Promise.all([...grouped].map(([partition, parts]) => {
			const length = parts.reduce((sum, part) => sum + part.length, 0), joined = new Uint8Array(length);
			let at = 0;
			for (const part of parts) {
				joined.set(part, at);
				at += part.length;
			}
			const operation = this.chains[partition].then(async () => {
				const result = await this.handles[partition].write(joined, 0, joined.length, this.lengths[partition]);
				if (result.bytesWritten !== joined.length) throw new Error("Temporary partition write was incomplete.");
				this.lengths[partition] += result.bytesWritten;
			});
			this.chains[partition] = operation;
			return operation;
		}));
	}
	async readFully(partition, output, position) {
		let offset = 0;
		while (offset < output.byteLength) {
			const result = await this.handles[partition].read(output, offset, output.byteLength - offset, position + offset);
			if (!result.bytesRead) throw new Error("Temporary partition read was incomplete.");
			offset += result.bytesRead;
		}
	}
	async readSelected(partition, cutoffs) {
		await this.chains[partition];
		const records = [], headerBytes = new Uint8Array(24);
		let offset = 0;
		while (offset < this.lengths[partition]) {
			await this.readFully(partition, headerBytes, offset);
			const header = parseSpoolRecordHeader(headerBytes);
			if (offset + header.recordLength > this.lengths[partition]) throw new Error("Temporary partition contains a truncated spool record.");
			if (selectedSpoolRecord(header, cutoffs)) {
				const record = new Uint8Array(header.recordLength);
				record.set(headerBytes);
				await this.readFully(partition, record.subarray(24), offset + 24);
				records.push(record);
			}
			offset += header.recordLength;
		}
		if (offset !== this.lengths[partition]) throw new Error("Temporary partition contains trailing spool bytes.");
		return concatenateSpoolRecords(records);
	}
	async close() {
		if (this.closed) return;
		this.closed = true;
		await Promise.all(this.chains);
		await Promise.all(this.handles.map((handle) => handle.close()));
		await rm(this.directory, {
			recursive: true,
			force: true
		});
	}
};
async function* fastqBatches(path, hash, batchRecords = 256, maximumBatchBytes = 4 * 1024 * 1024) {
	const raw = createReadStream(path, { highWaterMark: 1024 * 1024 });
	raw.on("data", (chunk) => hash.update(chunk));
	const lines = createInterface({
		input: path.toLowerCase().endsWith(".gz") ? raw.pipe(createGunzip({ chunkSize: 1024 * 1024 })) : raw,
		crlfDelay: Infinity
	});
	let record = [], batch = "", count = 0, records = 0, firstOrdinal = 1;
	for await (const line of lines) {
		record.push(line);
		if (record.length !== 4) continue;
		if (!record[0].startsWith("@") || !record[2].startsWith("+") || record[1].length !== record[3].length) throw new Error(`Malformed FASTQ record ${records + 1}.`);
		batch += `${record.join("\n")}\n`;
		record = [];
		count++;
		records++;
		if (count >= batchRecords || batch.length >= maximumBatchBytes) {
			yield {
				text: batch,
				count,
				firstOrdinal,
				records
			};
			firstOrdinal += count;
			batch = "";
			count = 0;
		}
	}
	if (record.length) throw new Error("The FASTQ stream is truncated at its final record.");
	if (count) yield {
		text: batch,
		count,
		firstOrdinal,
		records
	};
	if (!records) throw new Error("The input contains no FASTQ records.");
}
async function loadConfiguration(path) {
	const source = await readFile(path, "utf8"), config = parseConfigYaml(source), base = dirname(resolve(path));
	const requested = /* @__PURE__ */ new Map();
	const add = (name, role) => {
		if (!requested.has(name)) requested.set(name, role);
	};
	for (const sample of config.samples) {
		add(sample.panel, "panel");
		if (sample.functionalReference) add(sample.functionalReference, "functional-reference");
	}
	if (config.parameters.contaminationFilter) add(config.contaminationPanel, "contamination-panel");
	const resolved = new Map([...requested].map(([name]) => [name, isAbsolute(name) ? name : resolve(base, name)]));
	const loaded = await resolveReferenceFiles(config, new Map([...resolved].map(([name, filePath]) => [name, () => readFile(filePath, "utf8")])));
	const mappings = [{
		slot: "configuration",
		role: "configuration",
		uploadedName: basename(path),
		uploadedSize: Buffer.byteLength(source)
	}];
	for (const [name, role] of requested) {
		const filePath = resolved.get(name), information = await stat(filePath);
		mappings.push({
			slot: name,
			role,
			expectedName: name,
			uploadedName: basename(filePath),
			uploadedSize: information.size
		});
	}
	return {
		config: loaded,
		mappings
	};
}
async function runPipeline({ inputPath, configPath, outputPath, workers, assets, deferPhylogeny = false }) {
	const runStarted = performance.now(), timings = [];
	const input = resolve(inputPath), configuration = resolve(configPath), inputInformation = await stat(input);
	const loadedConfiguration = await loadConfiguration(configuration), config = loadedConfiguration.config, compiledConfig = compileConfig(config);
	const configHash = createHash("sha256").update(compiledConfig).digest("hex"), inputHash = createHash("sha256");
	status(`starting ${config.dataset} with ${workers} workers`);
	const pool = await WorkerPool.create(workers, assets.wasmPath, compiledConfig), store = await DiskPartitions.create(config.parameters.spoolPartitions);
	const log = [
		`${now()} webPORPID ${VERSION} started`,
		`${now()} execution: ${workers} WASM workers; disk-backed partition spool`,
		`${now()} parameters: error_rate=${config.parameters.errorRate}, lengths=(${config.parameters.minLength},${config.parameters.maxLength}), lda=${config.parameters.ldaThreshold}`
	];
	const inputMappings = [{
		slot: "reads",
		role: "reads",
		uploadedName: basename(input),
		uploadedSize: inputInformation.size
	}, ...loadedConfiguration.mappings];
	for (const mapping of inputMappings) log.push(`${now()} input mapping: ${mapping.role} slot ${mapping.slot}${mapping.expectedName ? ` (${mapping.expectedName})` : ""} <- ${mapping.uploadedName} (${mapping.uploadedSize} bytes)`);
	const storeTiming = (stage, seconds, workItems) => {
		const entry = {
			stage,
			seconds
		};
		if (workItems != null) entry.workItems = workItems;
		timings.push(entry);
		log.push(`${now()} timing ${stage}: ${entry.seconds.toFixed(6)} s${workItems == null ? "" : `; ${workItems} work items`}`);
	};
	const recordTiming = (stage, started, workItems) => {
		storeTiming(stage, (performance.now() - started) / 1e3, workItems);
		return performance.now();
	};
	recordTiming("setup", runStarted);
	try {
		let stageStarted = performance.now();
		const pending = /* @__PURE__ */ new Set();
		let batches = 0, streamed = 0;
		for await (const batch of fastqBatches(input, inputHash)) {
			const task = pool.any({
				type: "preprocess",
				text: batch.text,
				firstOrdinal: batch.firstOrdinal
			}).then((buffer) => store.appendFrames(new Uint8Array(buffer)));
			pending.add(task);
			task.then(() => pending.delete(task), () => pending.delete(task));
			batches++;
			streamed = batch.records;
			if (batches % 20 === 0) status(`preprocessing: ${streamed.toLocaleString()} reads streamed`);
			if (pending.size >= workers * 2) await Promise.race(pending);
		}
		await Promise.all(pending);
		const quality = mergeStats(await Promise.all(pool.clients.map((_, index) => pool.at(index, { type: "stats" }))), config.samples.length);
		log.push(`${now()} preprocessing: ${quality.totalReads} raw; ${quality.qualityReads} quality; ${quality.demultiplexedReads} demultiplexed; ${batches} bounded batches`);
		status(`preprocessing complete: ${quality.demultiplexedReads.toLocaleString()} demultiplexed reads`);
		stageStarted = recordTiming("preprocessing", stageStarted, quality.totalReads);
		const sampleCounts = quality.perSample.map(BigInt), cutoffValues = makeCutoffValues(sampleCounts, config.parameters.maxReadsPerSample);
		const cutoffs = makeCutoffs(sampleCounts, config.parameters.maxReadsPerSample), countParts = Array(config.parameters.spoolPartitions);
		await Promise.all(pool.clients.map(async (_, worker) => {
			for (let partition = worker; partition < countParts.length; partition += workers) {
				const bytes = await store.readSelected(partition, cutoffValues), data = bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength), cutoffCopy = cutoffs.slice().buffer;
				countParts[partition] = new Uint8Array(await pool.at(worker, {
					type: "countFamilies",
					bytes: data,
					cutoffs: cutoffCopy
				}, [data, cutoffCopy]));
			}
		}));
		const mergedCounts = mergeFamilyCounts(countParts), decodedCounts = decodeFamilyCounts(mergedCounts);
		const selectedReadsBySample = Array(config.samples.length).fill(0);
		for (const entry of decodedCounts) selectedReadsBySample[entry.sample] += entry.count;
		const selectedReads = selectedReadsBySample.reduce((sum, count) => sum + count, 0);
		quality.downsampledReads = Math.max(0, quality.demultiplexedReads - selectedReads);
		const modelData = mergedCounts.buffer.slice(mergedCounts.byteOffset, mergedCounts.byteOffset + mergedCounts.byteLength);
		const familyModel = new Uint8Array(await pool.at(0, {
			type: "buildModel",
			bytes: modelData
		}, [modelData]));
		const umiFamilies = decodeFamilyModel(familyModel, config);
		quality.bpbRejects = umiFamilies.filter((row) => row.disposition === "BPB-rejects").reduce((sum, row) => sum + row.familySize, 0);
		await Promise.all(pool.clients.map((_, index) => {
			const copy = familyModel.slice().buffer;
			return pool.at(index, {
				type: "initModel",
				bytes: copy
			}, [copy]);
		}));
		log.push(`${now()} UMI model: ${umiFamilies.filter((row) => row.disposition !== "BPB-rejects").length} observed families; ${quality.bpbRejects} BPB rejects; ${umiFamilies.filter((row) => row.disposition === "likely_real").length} initially likely real`);
		status(`UMI model complete: ${umiFamilies.filter((row) => row.disposition !== "BPB-rejects").length.toLocaleString()} observed families`);
		stageStarted = recordTiming("umi", stageStarted, quality.demultiplexedReads - quality.downsampledReads);
		const consensusParts = Array(countParts.length);
		await Promise.all(pool.clients.map(async (_, worker) => {
			for (let partition = worker; partition < consensusParts.length; partition += workers) {
				const bytes = await store.readSelected(partition, cutoffValues), data = bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength), cutoffCopy = cutoffs.slice().buffer;
				if (process.env.WEBPORPID_DEBUG) status(`debug: worker ${worker} consensus partition ${partition} bytes=${bytes.byteLength}`);
				const response = await pool.at(worker, {
					type: "consensus",
					bytes: data,
					cutoffs: cutoffCopy
				}, [data, cutoffCopy]);
				consensusParts[partition] = decodeConsensusOutput(new Uint8Array(response), config);
				if (process.env.WEBPORPID_DEBUG) status(`debug: worker ${worker} finished partition ${partition}`);
			}
		}));
		const consensuses = consensusParts.flatMap((part) => part.consensuses).sort((a, b) => a.sampleIndex - b.sampleIndex || a.umi.localeCompare(b.umi));
		const heteroduplexes = new Set(consensusParts.flatMap((part) => part.heteroduplexes));
		const consensusByFamily = new Map(consensuses.map((record) => [`${record.sampleIndex}\0${record.umi}`, record]));
		for (const family of umiFamilies) {
			const key = `${family.sampleIndex}\0${family.umi}`;
			if (heteroduplexes.has(key)) family.disposition = "heteroduplex";
			const consensus = consensusByFamily.get(key);
			if (consensus) family.minimumAgreement = consensus.minimumAgreement;
		}
		log.push(`${now()} consensus: ${consensuses.length} sequences; ${heteroduplexes.size} heteroduplex families`);
		status(`consensus complete: ${consensuses.length.toLocaleString()} sequences`);
		await store.close();
		stageStarted = recordTiming("consensus", stageStarted, consensuses.length);
		const contamination = classifyContamination(consensuses, config);
		log.push(`${now()} contamination: ${contamination.filter((call) => call.discarded).length} discarded; ${contamination.filter((call) => call.suspectOnly).length} suspect calls`);
		stageStarted = recordTiming("contamination", stageStarted, consensuses.length);
		const msaRunner = createMsaRunner(assets.msaPath, Math.min(workers, config.samples.length), assets.msaWorkerPath);
		const downstreamStarted = performance.now();
		let downstream;
		try {
			downstream = await postprocess(consensuses, contamination, config, void 0, msaRunner, workers);
		} finally {
			await msaRunner.close?.();
		}
		downstream.summaries.forEach((summary, index) => {
			summary.demultiplexedReads = quality.perSample[index] ?? 0;
			summary.selectedReads = selectedReadsBySample[index] ?? 0;
			summary.downsampledReads = Math.max(0, summary.demultiplexedReads - summary.selectedReads);
			summary.observedUmis = umiFamilies.filter((family) => family.sampleIndex === index && family.disposition !== "BPB-rejects").length;
			summary.likelyRealUmis = umiFamilies.filter((family) => family.sampleIndex === index && family.disposition === "likely_real").length;
		});
		const downstreamFinished = performance.now(), downstreamSeconds = (downstreamFinished - downstreamStarted) / 1e3;
		const collapseSeconds = Math.max(0, Math.min(downstreamSeconds, downstream.collapseSeconds));
		storeTiming("postprocessing", downstreamSeconds - collapseSeconds, downstream.records.length);
		const collapsedHaplotypes = Object.values(downstream.collapseGroups).reduce((sum, groups) => sum + groups.length, 0);
		storeTiming("collapse", collapseSeconds, collapsedHaplotypes);
		stageStarted = downstreamFinished;
		log.push(`${now()} collapse: ${collapsedHaplotypes} distinct haplotypes from ${downstream.records.filter((record) => record.alignedNt).length} retained UMI families; multiplicities count families, not reads`);
		const treeInputs = Object.entries(downstream.alignments).filter(([name]) => name.endsWith("/nucleotide"));
		let treeEntries = [];
		if (!deferPhylogeny) {
			const fastTree = createFastTreeRunner(assets.fastTreeJavascriptPath, assets.fastTreeWasmPath, Math.min(workers, Math.max(1, treeInputs.length)), assets.fastTreeWorkerPath);
			try {
				treeEntries = await Promise.all(treeInputs.map(async ([name, alignment]) => {
					status(`FastTree: ${name.split("/")[0]}`);
					return [name, await fastTree(alignment)];
				}));
			} finally {
				await fastTree.close?.();
			}
		} else log.push(`${now()} phylogeny: deferred by user; collapsed alignments are stored and trees can be inferred in the results explorer`);
		const trees = Object.fromEntries(treeEntries);
		recordTiming("tree", stageStarted, Object.keys(trees).length);
		timings.push({
			stage: "analysis-total",
			seconds: (performance.now() - runStarted) / 1e3
		});
		log.push(`${now()} postprocessing: ${downstream.records.filter((record) => record.artefactPass && record.agreementPass && record.contaminationPass && record.panelPass).length} sequences passed all non-functional filters`);
		const result = {
			schema: "webporpid-results/1",
			provenance: {
				webporpidVersion: VERSION,
				createdUtc: now(),
				engine: "C++20 WASM/WASI SIMD",
				workers,
				inputName: basename(input),
				inputSha256: inputHash.digest("hex"),
				configSha256: configHash,
				deterministicSeed: config.parameters.deterministicSeed.toString(),
				upstreamBranch: "nanopore",
				upstreamCommit: UPSTREAM_COMMIT
			},
			config: resultConfig(config),
			quality,
			summaries: downstream.summaries,
			umiFamilies,
			consensuses,
			contamination,
			contaminationReferences: config.contaminationPanelSequences,
			downstreamResources: downstreamResources(config),
			records: downstream.records,
			alignments: downstream.alignments,
			trees,
			referenceAlignments: downstream.referenceAlignments,
			collapseGroups: downstream.collapseGroups,
			inputMappings,
			runOptions: {
				deferPhylogeny,
				deferContamination: false,
				deferPostprocessing: false,
				deferCollapse: false
			},
			optionalStages: {
				contamination: statusRecord("completed", `${contamination.filter((call) => call.discarded).length} consensus sequences excluded.`),
				postprocessing: statusRecord("completed", `${downstream.records.length} consensus-family records evaluated.`),
				collapse: statusRecord("completed", `${collapsedHaplotypes} haplotypes; multiplicities count UMI families.`),
				tree: statusRecord(deferPhylogeny ? "deferred" : "completed", deferPhylogeny ? "Deferred by user; collapsed alignments are stored for on-demand inference." : `${Object.keys(trees).length} collapsed phylogenies inferred.`)
			},
			timings,
			log
		};
		await mkdir(dirname(resolve(outputPath)), { recursive: true });
		await writeFile(outputPath, encodeResult(result));
		status(`wrote ${outputPath}`);
		return result;
	} finally {
		await pool.close();
		await store.close();
	}
}
async function inspect(path) {
	const result = decodeResult(new Uint8Array(await readFile(path)));
	process.stdout.write(JSON.stringify({
		schema: result.schema,
		provenance: result.provenance,
		quality: result.quality,
		summaries: result.summaries,
		timings: result.timings ?? [],
		components: {
			consensuses: result.consensuses.length,
			families: result.umiFamilies.length,
			contaminationCalls: result.contamination.length,
			contaminationReferences: result.contaminationReferences?.length ?? 0,
			records: result.records.length,
			collapsedHaplotypes: Object.values(result.collapseGroups ?? {}).reduce((sum, groups) => sum + groups.length, 0),
			alignments: Object.keys(result.alignments),
			trees: Object.keys(result.trees)
		}
	}, null, 2) + "\n");
}
async function exportResult(args) {
	const path = args[0], component = option(args, "--component"), sample = option(args, "--sample");
	if (!path || !component) throw new Error("export requires a results file and --component.");
	const result = decodeResult(new Uint8Array(await readFile(path))), exported = exportComponent(result, component, sample);
	const output = option(args, "--output") ?? `${safeDatasetName(result.config.dataset)}${sample ? `-${safeDatasetName(sample)}` : ""}.${exported.extension}`;
	await writeFile(output, exported.text);
	status(`wrote ${output}`);
}
async function runCli(overrideAssets) {
	const args = process.argv.slice(2), command = args.shift();
	if (!command || command === "--help" || command === "-h" || command === "help") {
		process.stdout.write(`${usage()}\n`);
		return;
	}
	if (command === "--version" || command === "version") {
		process.stdout.write(`${VERSION}\n`);
		return;
	}
	if (command === "inspect") {
		if (!args[0]) throw new Error("inspect requires a .webporpid file.");
		await inspect(args[0]);
		return;
	}
	if (command === "export") {
		await exportResult(args);
		return;
	}
	if (command !== "run") throw new Error(`Unknown command ${command}.\n\n${usage()}`);
	const inputPath = args[0], configPath = option(args, "--config");
	if (!inputPath || inputPath.startsWith("--") || !configPath) throw new Error("run requires an input FASTQ and --config config.yaml.");
	await runPipeline({
		inputPath,
		configPath,
		outputPath: option(args, "--output") ?? option(args, "--out") ?? `${basename(inputPath).replace(/\.(fastq|fq)(\.gz)?$/i, "")}.webporpid`,
		workers: integer(option(args, "--workers"), "--workers", availableParallelism()),
		assets: overrideAssets ?? defaultCliAssets(),
		deferPhylogeny: args.includes("--defer-phylogeny")
	});
}
//#endregion
//#region cli-src/porpid-cli-node.mjs
runCli().catch((cause) => {
	process.stderr.write(`porpid-cli: ${cause instanceof Error ? cause.message : String(cause)}\n`);
	process.exitCode = 1;
});
//#endregion
