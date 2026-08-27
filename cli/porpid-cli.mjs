#!/usr/bin/env node
import { a as makeCutoffValues, c as mergeStats, i as decodeFamilyModel, l as WASI, n as decodeConsensusOutput, o as makeCutoffs, r as decodeFamilyCounts, s as mergeFamilyCounts, u as BinaryWriter } from "./chunks/wasm-runtime-CjLR_5MC.mjs";
import { createHash } from "node:crypto";
import { createReadStream } from "node:fs";
import { mkdir, mkdtemp, open, readFile, rm, stat, writeFile } from "node:fs/promises";
import { availableParallelism, tmpdir } from "node:os";
import { basename, dirname, isAbsolute, join, resolve } from "node:path";
import { createInterface } from "node:readline";
import { createGunzip } from "node:zlib";
import { Worker as Worker$1 } from "node:worker_threads";
import { createRequire } from "node:module";
//#region \0rolldown/runtime.js
var __defProp = Object.defineProperty;
var __exportAll = (all, no_symbols) => {
	let target = {};
	for (var name in all) __defProp(target, name, {
		get: all[name],
		enumerable: true
	});
	if (!no_symbols) __defProp(target, Symbol.toStringTag, { value: "Module" });
	return target;
};
//#endregion
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
		if (possible) output.push({
			sample: record.sample,
			sequenceId: record.id,
			nearestNonselfVariant: possible.label,
			nearestNonselfDistance: possible.distance,
			flagged: true,
			discarded: false,
			suspectOnly: true
		});
	}
	return output.sort((a, b) => a.nearestNonselfDistance - b.nearestNonselfDistance);
}
//#endregion
//#region node_modules/yaml/browser/dist/nodes/identity.js
const ALIAS = Symbol.for("yaml.alias");
const DOC = Symbol.for("yaml.document");
const MAP = Symbol.for("yaml.map");
const PAIR = Symbol.for("yaml.pair");
const SCALAR$1 = Symbol.for("yaml.scalar");
const SEQ = Symbol.for("yaml.seq");
const NODE_TYPE = Symbol.for("yaml.node.type");
const isAlias = (node) => !!node && typeof node === "object" && node[NODE_TYPE] === ALIAS;
const isDocument = (node) => !!node && typeof node === "object" && node[NODE_TYPE] === DOC;
const isMap = (node) => !!node && typeof node === "object" && node[NODE_TYPE] === MAP;
const isPair = (node) => !!node && typeof node === "object" && node[NODE_TYPE] === PAIR;
const isScalar$1 = (node) => !!node && typeof node === "object" && node[NODE_TYPE] === SCALAR$1;
const isSeq = (node) => !!node && typeof node === "object" && node[NODE_TYPE] === SEQ;
function isCollection$1(node) {
	if (node && typeof node === "object") switch (node[NODE_TYPE]) {
		case MAP:
		case SEQ: return true;
	}
	return false;
}
function isNode(node) {
	if (node && typeof node === "object") switch (node[NODE_TYPE]) {
		case ALIAS:
		case MAP:
		case SCALAR$1:
		case SEQ: return true;
	}
	return false;
}
const hasAnchor = (node) => (isScalar$1(node) || isCollection$1(node)) && !!node.anchor;
//#endregion
//#region node_modules/yaml/browser/dist/visit.js
const BREAK$1 = Symbol("break visit");
const SKIP$1 = Symbol("skip children");
const REMOVE$1 = Symbol("remove node");
/**
* Apply a visitor to an AST node or document.
*
* Walks through the tree (depth-first) starting from `node`, calling a
* `visitor` function with three arguments:
*   - `key`: For sequence values and map `Pair`, the node's index in the
*     collection. Within a `Pair`, `'key'` or `'value'`, correspondingly.
*     `null` for the root node.
*   - `node`: The current node.
*   - `path`: The ancestry of the current node.
*
* The return value of the visitor may be used to control the traversal:
*   - `undefined` (default): Do nothing and continue
*   - `visit.SKIP`: Do not visit the children of this node, continue with next
*     sibling
*   - `visit.BREAK`: Terminate traversal completely
*   - `visit.REMOVE`: Remove the current node, then continue with the next one
*   - `Node`: Replace the current node, then continue by visiting it
*   - `number`: While iterating the items of a sequence or map, set the index
*     of the next step. This is useful especially if the index of the current
*     node has changed.
*
* If `visitor` is a single function, it will be called with all values
* encountered in the tree, including e.g. `null` values. Alternatively,
* separate visitor functions may be defined for each `Map`, `Pair`, `Seq`,
* `Alias` and `Scalar` node. To define the same visitor function for more than
* one node type, use the `Collection` (map and seq), `Value` (map, seq & scalar)
* and `Node` (alias, map, seq & scalar) targets. Of all these, only the most
* specific defined one will be used for each node.
*/
function visit$1(node, visitor) {
	const visitor_ = initVisitor(visitor);
	if (isDocument(node)) {
		if (visit_(null, node.contents, visitor_, Object.freeze([node])) === REMOVE$1) node.contents = null;
	} else visit_(null, node, visitor_, Object.freeze([]));
}
/** Terminate visit traversal completely */
visit$1.BREAK = BREAK$1;
/** Do not visit the children of the current node */
visit$1.SKIP = SKIP$1;
/** Remove the current node */
visit$1.REMOVE = REMOVE$1;
function visit_(key, node, visitor, path) {
	const ctrl = callVisitor(key, node, visitor, path);
	if (isNode(ctrl) || isPair(ctrl)) {
		replaceNode(key, path, ctrl);
		return visit_(key, ctrl, visitor, path);
	}
	if (typeof ctrl !== "symbol") {
		if (isCollection$1(node)) {
			path = Object.freeze(path.concat(node));
			for (let i = 0; i < node.items.length; ++i) {
				const ci = visit_(i, node.items[i], visitor, path);
				if (typeof ci === "number") i = ci - 1;
				else if (ci === BREAK$1) return BREAK$1;
				else if (ci === REMOVE$1) {
					node.items.splice(i, 1);
					i -= 1;
				}
			}
		} else if (isPair(node)) {
			path = Object.freeze(path.concat(node));
			const ck = visit_("key", node.key, visitor, path);
			if (ck === BREAK$1) return BREAK$1;
			else if (ck === REMOVE$1) node.key = null;
			const cv = visit_("value", node.value, visitor, path);
			if (cv === BREAK$1) return BREAK$1;
			else if (cv === REMOVE$1) node.value = null;
		}
	}
	return ctrl;
}
/**
* Apply an async visitor to an AST node or document.
*
* Walks through the tree (depth-first) starting from `node`, calling a
* `visitor` function with three arguments:
*   - `key`: For sequence values and map `Pair`, the node's index in the
*     collection. Within a `Pair`, `'key'` or `'value'`, correspondingly.
*     `null` for the root node.
*   - `node`: The current node.
*   - `path`: The ancestry of the current node.
*
* The return value of the visitor may be used to control the traversal:
*   - `Promise`: Must resolve to one of the following values
*   - `undefined` (default): Do nothing and continue
*   - `visit.SKIP`: Do not visit the children of this node, continue with next
*     sibling
*   - `visit.BREAK`: Terminate traversal completely
*   - `visit.REMOVE`: Remove the current node, then continue with the next one
*   - `Node`: Replace the current node, then continue by visiting it
*   - `number`: While iterating the items of a sequence or map, set the index
*     of the next step. This is useful especially if the index of the current
*     node has changed.
*
* If `visitor` is a single function, it will be called with all values
* encountered in the tree, including e.g. `null` values. Alternatively,
* separate visitor functions may be defined for each `Map`, `Pair`, `Seq`,
* `Alias` and `Scalar` node. To define the same visitor function for more than
* one node type, use the `Collection` (map and seq), `Value` (map, seq & scalar)
* and `Node` (alias, map, seq & scalar) targets. Of all these, only the most
* specific defined one will be used for each node.
*/
async function visitAsync(node, visitor) {
	const visitor_ = initVisitor(visitor);
	if (isDocument(node)) {
		if (await visitAsync_(null, node.contents, visitor_, Object.freeze([node])) === REMOVE$1) node.contents = null;
	} else await visitAsync_(null, node, visitor_, Object.freeze([]));
}
/** Terminate visit traversal completely */
visitAsync.BREAK = BREAK$1;
/** Do not visit the children of the current node */
visitAsync.SKIP = SKIP$1;
/** Remove the current node */
visitAsync.REMOVE = REMOVE$1;
async function visitAsync_(key, node, visitor, path) {
	const ctrl = await callVisitor(key, node, visitor, path);
	if (isNode(ctrl) || isPair(ctrl)) {
		replaceNode(key, path, ctrl);
		return visitAsync_(key, ctrl, visitor, path);
	}
	if (typeof ctrl !== "symbol") {
		if (isCollection$1(node)) {
			path = Object.freeze(path.concat(node));
			for (let i = 0; i < node.items.length; ++i) {
				const ci = await visitAsync_(i, node.items[i], visitor, path);
				if (typeof ci === "number") i = ci - 1;
				else if (ci === BREAK$1) return BREAK$1;
				else if (ci === REMOVE$1) {
					node.items.splice(i, 1);
					i -= 1;
				}
			}
		} else if (isPair(node)) {
			path = Object.freeze(path.concat(node));
			const ck = await visitAsync_("key", node.key, visitor, path);
			if (ck === BREAK$1) return BREAK$1;
			else if (ck === REMOVE$1) node.key = null;
			const cv = await visitAsync_("value", node.value, visitor, path);
			if (cv === BREAK$1) return BREAK$1;
			else if (cv === REMOVE$1) node.value = null;
		}
	}
	return ctrl;
}
function initVisitor(visitor) {
	if (typeof visitor === "object" && (visitor.Collection || visitor.Node || visitor.Value)) return Object.assign({
		Alias: visitor.Node,
		Map: visitor.Node,
		Scalar: visitor.Node,
		Seq: visitor.Node
	}, visitor.Value && {
		Map: visitor.Value,
		Scalar: visitor.Value,
		Seq: visitor.Value
	}, visitor.Collection && {
		Map: visitor.Collection,
		Seq: visitor.Collection
	}, visitor);
	return visitor;
}
function callVisitor(key, node, visitor, path) {
	if (typeof visitor === "function") return visitor(key, node, path);
	if (isMap(node)) return visitor.Map?.(key, node, path);
	if (isSeq(node)) return visitor.Seq?.(key, node, path);
	if (isPair(node)) return visitor.Pair?.(key, node, path);
	if (isScalar$1(node)) return visitor.Scalar?.(key, node, path);
	if (isAlias(node)) return visitor.Alias?.(key, node, path);
}
function replaceNode(key, path, node) {
	const parent = path[path.length - 1];
	if (isCollection$1(parent)) parent.items[key] = node;
	else if (isPair(parent)) if (key === "key") parent.key = node;
	else parent.value = node;
	else if (isDocument(parent)) parent.contents = node;
	else {
		const pt = isAlias(parent) ? "alias" : "scalar";
		throw new Error(`Cannot replace node with ${pt} parent`);
	}
}
//#endregion
//#region node_modules/yaml/browser/dist/doc/directives.js
const escapeChars = {
	"!": "%21",
	",": "%2C",
	"[": "%5B",
	"]": "%5D",
	"{": "%7B",
	"}": "%7D"
};
const escapeTagName = (tn) => tn.replace(/[!,[\]{}]/g, (ch) => escapeChars[ch]);
var Directives = class Directives {
	constructor(yaml, tags) {
		/**
		* The directives-end/doc-start marker `---`. If `null`, a marker may still be
		* included in the document's stringified representation.
		*/
		this.docStart = null;
		/** The doc-end marker `...`.  */
		this.docEnd = false;
		this.yaml = Object.assign({}, Directives.defaultYaml, yaml);
		this.tags = Object.assign({}, Directives.defaultTags, tags);
	}
	clone() {
		const copy = new Directives(this.yaml, this.tags);
		copy.docStart = this.docStart;
		return copy;
	}
	/**
	* During parsing, get a Directives instance for the current document and
	* update the stream state according to the current version's spec.
	*/
	atDocument() {
		const res = new Directives(this.yaml, this.tags);
		switch (this.yaml.version) {
			case "1.1":
				this.atNextDocument = true;
				break;
			case "1.2":
				this.atNextDocument = false;
				this.yaml = {
					explicit: Directives.defaultYaml.explicit,
					version: "1.2"
				};
				this.tags = Object.assign({}, Directives.defaultTags);
				break;
		}
		return res;
	}
	/**
	* @param onError - May be called even if the action was successful
	* @returns `true` on success
	*/
	add(line, onError) {
		if (this.atNextDocument) {
			this.yaml = {
				explicit: Directives.defaultYaml.explicit,
				version: "1.1"
			};
			this.tags = Object.assign({}, Directives.defaultTags);
			this.atNextDocument = false;
		}
		const parts = line.trim().split(/[ \t]+/);
		const name = parts.shift();
		switch (name) {
			case "%TAG": {
				if (parts.length !== 2) {
					onError(0, "%TAG directive should contain exactly two parts");
					if (parts.length < 2) return false;
				}
				const [handle, prefix] = parts;
				this.tags[handle] = prefix;
				return true;
			}
			case "%YAML": {
				this.yaml.explicit = true;
				if (parts.length !== 1) {
					onError(0, "%YAML directive should contain exactly one part");
					return false;
				}
				const [version] = parts;
				if (version === "1.1" || version === "1.2") {
					this.yaml.version = version;
					return true;
				} else {
					const isValid = /^\d+\.\d+$/.test(version);
					onError(6, `Unsupported YAML version ${version}`, isValid);
					return false;
				}
			}
			default:
				onError(0, `Unknown directive ${name}`, true);
				return false;
		}
	}
	/**
	* Resolves a tag, matching handles to those defined in %TAG directives.
	*
	* @returns Resolved tag, which may also be the non-specific tag `'!'` or a
	*   `'!local'` tag, or `null` if unresolvable.
	*/
	tagName(source, onError) {
		if (source === "!") return "!";
		if (source[0] !== "!") {
			onError(`Not a valid tag: ${source}`);
			return null;
		}
		if (source[1] === "<") {
			const verbatim = source.slice(2, -1);
			if (verbatim === "!" || verbatim === "!!") {
				onError(`Verbatim tags aren't resolved, so ${source} is invalid.`);
				return null;
			}
			if (source[source.length - 1] !== ">") onError("Verbatim tags must end with a >");
			return verbatim;
		}
		const [, handle, suffix] = source.match(/^(.*!)([^!]*)$/s);
		if (!suffix) onError(`The ${source} tag has no suffix`);
		const prefix = this.tags[handle];
		if (prefix) try {
			return prefix + decodeURIComponent(suffix);
		} catch (error) {
			onError(String(error));
			return null;
		}
		if (handle === "!") return source;
		onError(`Could not resolve tag: ${source}`);
		return null;
	}
	/**
	* Given a fully resolved tag, returns its printable string form,
	* taking into account current tag prefixes and defaults.
	*/
	tagString(tag) {
		for (const [handle, prefix] of Object.entries(this.tags)) if (tag.startsWith(prefix)) return handle + escapeTagName(tag.substring(prefix.length));
		return tag[0] === "!" ? tag : `!<${tag}>`;
	}
	toString(doc) {
		const lines = this.yaml.explicit ? [`%YAML ${this.yaml.version || "1.2"}`] : [];
		const tagEntries = Object.entries(this.tags);
		let tagNames;
		if (doc && tagEntries.length > 0 && isNode(doc.contents)) {
			const tags = {};
			visit$1(doc.contents, (_key, node) => {
				if (isNode(node) && node.tag) tags[node.tag] = true;
			});
			tagNames = Object.keys(tags);
		} else tagNames = [];
		for (const [handle, prefix] of tagEntries) {
			if (handle === "!!" && prefix === "tag:yaml.org,2002:") continue;
			if (!doc || tagNames.some((tn) => tn.startsWith(prefix))) lines.push(`%TAG ${handle} ${prefix}`);
		}
		return lines.join("\n");
	}
};
Directives.defaultYaml = {
	explicit: false,
	version: "1.2"
};
Directives.defaultTags = { "!!": "tag:yaml.org,2002:" };
//#endregion
//#region node_modules/yaml/browser/dist/doc/anchors.js
/**
* Verify that the input string is a valid anchor.
*
* Will throw on errors.
*/
function anchorIsValid(anchor) {
	if (/[\x00-\x19\s,[\]{}]/.test(anchor)) {
		const msg = `Anchor must not contain whitespace or control characters: ${JSON.stringify(anchor)}`;
		throw new Error(msg);
	}
	return true;
}
function anchorNames(root) {
	const anchors = /* @__PURE__ */ new Set();
	visit$1(root, { Value(_key, node) {
		if (node.anchor) anchors.add(node.anchor);
	} });
	return anchors;
}
/** Find a new anchor name with the given `prefix` and a one-indexed suffix. */
function findNewAnchor(prefix, exclude) {
	for (let i = 1;; ++i) {
		const name = `${prefix}${i}`;
		if (!exclude.has(name)) return name;
	}
}
function createNodeAnchors(doc, prefix) {
	const aliasObjects = [];
	const sourceObjects = /* @__PURE__ */ new Map();
	let prevAnchors = null;
	return {
		onAnchor: (source) => {
			aliasObjects.push(source);
			prevAnchors ?? (prevAnchors = anchorNames(doc));
			const anchor = findNewAnchor(prefix, prevAnchors);
			prevAnchors.add(anchor);
			return anchor;
		},
		/**
		* With circular references, the source node is only resolved after all
		* of its child nodes are. This is why anchors are set only after all of
		* the nodes have been created.
		*/
		setAnchors: () => {
			for (const source of aliasObjects) {
				const ref = sourceObjects.get(source);
				if (typeof ref === "object" && ref.anchor && (isScalar$1(ref.node) || isCollection$1(ref.node))) ref.node.anchor = ref.anchor;
				else {
					const error = /* @__PURE__ */ new Error("Failed to resolve repeated object (this should not happen)");
					error.source = source;
					throw error;
				}
			}
		},
		sourceObjects
	};
}
//#endregion
//#region node_modules/yaml/browser/dist/doc/applyReviver.js
/**
* Applies the JSON.parse reviver algorithm as defined in the ECMA-262 spec,
* in section 24.5.1.1 "Runtime Semantics: InternalizeJSONProperty" of the
* 2021 edition: https://tc39.es/ecma262/#sec-json.parse
*
* Includes extensions for handling Map and Set objects.
*/
function applyReviver(reviver, obj, key, val) {
	if (val && typeof val === "object") if (Array.isArray(val)) for (let i = 0, len = val.length; i < len; ++i) {
		const v0 = val[i];
		const v1 = applyReviver(reviver, val, String(i), v0);
		if (v1 === void 0) delete val[i];
		else if (v1 !== v0) val[i] = v1;
	}
	else if (val instanceof Map) for (const k of Array.from(val.keys())) {
		const v0 = val.get(k);
		const v1 = applyReviver(reviver, val, k, v0);
		if (v1 === void 0) val.delete(k);
		else if (v1 !== v0) val.set(k, v1);
	}
	else if (val instanceof Set) for (const v0 of Array.from(val)) {
		const v1 = applyReviver(reviver, val, v0, v0);
		if (v1 === void 0) val.delete(v0);
		else if (v1 !== v0) {
			val.delete(v0);
			val.add(v1);
		}
	}
	else for (const [k, v0] of Object.entries(val)) {
		const v1 = applyReviver(reviver, val, k, v0);
		if (v1 === void 0) delete val[k];
		else if (v1 !== v0) val[k] = v1;
	}
	return reviver.call(obj, key, val);
}
//#endregion
//#region node_modules/yaml/browser/dist/nodes/toJS.js
/**
* Recursively convert any node or its contents to native JavaScript
*
* @param value - The input value
* @param arg - If `value` defines a `toJSON()` method, use this
*   as its first argument
* @param ctx - Conversion context, originally set in Document#toJS(). If
*   `{ keep: true }` is not set, output should be suitable for JSON
*   stringification.
*/
function toJS(value, arg, ctx) {
	if (Array.isArray(value)) return value.map((v, i) => toJS(v, String(i), ctx));
	if (value && typeof value.toJSON === "function") {
		if (!ctx || !hasAnchor(value)) return value.toJSON(arg, ctx);
		const data = {
			aliasCount: 0,
			count: 1,
			res: void 0
		};
		ctx.anchors.set(value, data);
		ctx.onCreate = (res) => {
			data.res = res;
			delete ctx.onCreate;
		};
		const res = value.toJSON(arg, ctx);
		if (ctx.onCreate) ctx.onCreate(res);
		return res;
	}
	if (typeof value === "bigint" && !ctx?.keep) return Number(value);
	return value;
}
//#endregion
//#region node_modules/yaml/browser/dist/nodes/Node.js
var NodeBase = class {
	constructor(type) {
		Object.defineProperty(this, NODE_TYPE, { value: type });
	}
	/** Create a copy of this node.  */
	clone() {
		const copy = Object.create(Object.getPrototypeOf(this), Object.getOwnPropertyDescriptors(this));
		if (this.range) copy.range = this.range.slice();
		return copy;
	}
	/** A plain JavaScript representation of this node. */
	toJS(doc, { mapAsMap, maxAliasCount, onAnchor, reviver } = {}) {
		if (!isDocument(doc)) throw new TypeError("A document argument is required");
		const ctx = {
			anchors: /* @__PURE__ */ new Map(),
			doc,
			keep: true,
			mapAsMap: mapAsMap === true,
			mapKeyWarned: false,
			maxAliasCount: typeof maxAliasCount === "number" ? maxAliasCount : 100
		};
		const res = toJS(this, "", ctx);
		if (typeof onAnchor === "function") for (const { count, res } of ctx.anchors.values()) onAnchor(res, count);
		return typeof reviver === "function" ? applyReviver(reviver, { "": res }, "", res) : res;
	}
};
//#endregion
//#region node_modules/yaml/browser/dist/nodes/Alias.js
var Alias = class extends NodeBase {
	constructor(source) {
		super(ALIAS);
		this.source = source;
		Object.defineProperty(this, "tag", { set() {
			throw new Error("Alias nodes cannot have tags");
		} });
	}
	/**
	* Resolve the value of this alias within `doc`, finding the last
	* instance of the `source` anchor before this node.
	*/
	resolve(doc, ctx) {
		if (ctx?.maxAliasCount === 0) throw new ReferenceError("Alias resolution is disabled");
		let nodes;
		if (ctx?.aliasResolveCache) nodes = ctx.aliasResolveCache;
		else {
			nodes = [];
			visit$1(doc, { Node: (_key, node) => {
				if (isAlias(node) || hasAnchor(node)) nodes.push(node);
			} });
			if (ctx) ctx.aliasResolveCache = nodes;
		}
		let found = void 0;
		for (const node of nodes) {
			if (node === this) break;
			if (node.anchor === this.source) found = node;
		}
		return found;
	}
	toJSON(_arg, ctx) {
		if (!ctx) return { source: this.source };
		const { anchors, doc, maxAliasCount } = ctx;
		const source = this.resolve(doc, ctx);
		if (!source) {
			const msg = `Unresolved alias (the anchor must be set before the alias): ${this.source}`;
			throw new ReferenceError(msg);
		}
		let data = anchors.get(source);
		if (!data) {
			toJS(source, null, ctx);
			data = anchors.get(source);
		}
		/* istanbul ignore if */
		if (data?.res === void 0) throw new ReferenceError("This should not happen: Alias anchor was not resolved?");
		if (maxAliasCount >= 0) {
			data.count += 1;
			if (data.aliasCount === 0) data.aliasCount = getAliasCount(doc, source, anchors);
			if (data.count * data.aliasCount > maxAliasCount) throw new ReferenceError("Excessive alias count indicates a resource exhaustion attack");
		}
		return data.res;
	}
	toString(ctx, _onComment, _onChompKeep) {
		const src = `*${this.source}`;
		if (ctx) {
			anchorIsValid(this.source);
			if (ctx.options.verifyAliasOrder && !ctx.anchors.has(this.source)) {
				const msg = `Unresolved alias (the anchor must be set before the alias): ${this.source}`;
				throw new Error(msg);
			}
			if (ctx.implicitKey) return `${src} `;
		}
		return src;
	}
};
function getAliasCount(doc, node, anchors) {
	if (isAlias(node)) {
		const source = node.resolve(doc);
		const anchor = anchors && source && anchors.get(source);
		return anchor ? anchor.count * anchor.aliasCount : 0;
	} else if (isCollection$1(node)) {
		let count = 0;
		for (const item of node.items) {
			const c = getAliasCount(doc, item, anchors);
			if (c > count) count = c;
		}
		return count;
	} else if (isPair(node)) {
		const kc = getAliasCount(doc, node.key, anchors);
		const vc = getAliasCount(doc, node.value, anchors);
		return Math.max(kc, vc);
	}
	return 1;
}
//#endregion
//#region node_modules/yaml/browser/dist/nodes/Scalar.js
const isScalarValue = (value) => !value || typeof value !== "function" && typeof value !== "object";
var Scalar = class extends NodeBase {
	constructor(value) {
		super(SCALAR$1);
		this.value = value;
	}
	toJSON(arg, ctx) {
		return ctx?.keep ? this.value : toJS(this.value, arg, ctx);
	}
	toString() {
		return String(this.value);
	}
};
Scalar.BLOCK_FOLDED = "BLOCK_FOLDED";
Scalar.BLOCK_LITERAL = "BLOCK_LITERAL";
Scalar.PLAIN = "PLAIN";
Scalar.QUOTE_DOUBLE = "QUOTE_DOUBLE";
Scalar.QUOTE_SINGLE = "QUOTE_SINGLE";
//#endregion
//#region node_modules/yaml/browser/dist/doc/createNode.js
const defaultTagPrefix = "tag:yaml.org,2002:";
function findTagObject(value, tagName, tags) {
	if (tagName) {
		const match = tags.filter((t) => t.tag === tagName);
		const tagObj = match.find((t) => !t.format) ?? match[0];
		if (!tagObj) throw new Error(`Tag ${tagName} not found`);
		return tagObj;
	}
	return tags.find((t) => t.identify?.(value) && !t.format);
}
function createNode(value, tagName, ctx) {
	if (isDocument(value)) value = value.contents;
	if (isNode(value)) return value;
	if (isPair(value)) {
		const map = ctx.schema[MAP].createNode?.(ctx.schema, null, ctx);
		map.items.push(value);
		return map;
	}
	if (value instanceof String || value instanceof Number || value instanceof Boolean || typeof BigInt !== "undefined" && value instanceof BigInt) value = value.valueOf();
	const { aliasDuplicateObjects, onAnchor, onTagObj, schema, sourceObjects } = ctx;
	let ref = void 0;
	if (aliasDuplicateObjects && value && typeof value === "object") {
		ref = sourceObjects.get(value);
		if (ref) {
			ref.anchor ?? (ref.anchor = onAnchor(value));
			return new Alias(ref.anchor);
		} else {
			ref = {
				anchor: null,
				node: null
			};
			sourceObjects.set(value, ref);
		}
	}
	if (tagName?.startsWith("!!")) tagName = defaultTagPrefix + tagName.slice(2);
	let tagObj = findTagObject(value, tagName, schema.tags);
	if (!tagObj) {
		if (value && typeof value.toJSON === "function") value = value.toJSON();
		if (!value || typeof value !== "object") {
			const node = new Scalar(value);
			if (ref) ref.node = node;
			return node;
		}
		tagObj = value instanceof Map ? schema[MAP] : Symbol.iterator in Object(value) ? schema[SEQ] : schema[MAP];
	}
	if (onTagObj) {
		onTagObj(tagObj);
		delete ctx.onTagObj;
	}
	const node = tagObj?.createNode ? tagObj.createNode(ctx.schema, value, ctx) : typeof tagObj?.nodeClass?.from === "function" ? tagObj.nodeClass.from(ctx.schema, value, ctx) : new Scalar(value);
	if (tagName) node.tag = tagName;
	else if (!tagObj.default) node.tag = tagObj.tag;
	if (ref) ref.node = node;
	return node;
}
//#endregion
//#region node_modules/yaml/browser/dist/nodes/Collection.js
function collectionFromPath(schema, path, value) {
	let v = value;
	for (let i = path.length - 1; i >= 0; --i) {
		const k = path[i];
		if (typeof k === "number" && Number.isInteger(k) && k >= 0) {
			const a = [];
			a[k] = v;
			v = a;
		} else v = new Map([[k, v]]);
	}
	return createNode(v, void 0, {
		aliasDuplicateObjects: false,
		keepUndefined: false,
		onAnchor: () => {
			throw new Error("This should not happen, please report a bug.");
		},
		schema,
		sourceObjects: /* @__PURE__ */ new Map()
	});
}
const isEmptyPath = (path) => path == null || typeof path === "object" && !!path[Symbol.iterator]().next().done;
var Collection = class extends NodeBase {
	constructor(type, schema) {
		super(type);
		Object.defineProperty(this, "schema", {
			value: schema,
			configurable: true,
			enumerable: false,
			writable: true
		});
	}
	/**
	* Create a copy of this collection.
	*
	* @param schema - If defined, overwrites the original's schema
	*/
	clone(schema) {
		const copy = Object.create(Object.getPrototypeOf(this), Object.getOwnPropertyDescriptors(this));
		if (schema) copy.schema = schema;
		copy.items = copy.items.map((it) => isNode(it) || isPair(it) ? it.clone(schema) : it);
		if (this.range) copy.range = this.range.slice();
		return copy;
	}
	/**
	* Adds a value to the collection. For `!!map` and `!!omap` the value must
	* be a Pair instance or a `{ key, value }` object, which may not have a key
	* that already exists in the map.
	*/
	addIn(path, value) {
		if (isEmptyPath(path)) this.add(value);
		else {
			const [key, ...rest] = path;
			const node = this.get(key, true);
			if (isCollection$1(node)) node.addIn(rest, value);
			else if (node === void 0 && this.schema) this.set(key, collectionFromPath(this.schema, rest, value));
			else throw new Error(`Expected YAML collection at ${key}. Remaining path: ${rest}`);
		}
	}
	/**
	* Removes a value from the collection.
	* @returns `true` if the item was found and removed.
	*/
	deleteIn(path) {
		const [key, ...rest] = path;
		if (rest.length === 0) return this.delete(key);
		const node = this.get(key, true);
		if (isCollection$1(node)) return node.deleteIn(rest);
		else throw new Error(`Expected YAML collection at ${key}. Remaining path: ${rest}`);
	}
	/**
	* Returns item at `key`, or `undefined` if not found. By default unwraps
	* scalar values from their surrounding node; to disable set `keepScalar` to
	* `true` (collections are always returned intact).
	*/
	getIn(path, keepScalar) {
		const [key, ...rest] = path;
		const node = this.get(key, true);
		if (rest.length === 0) return !keepScalar && isScalar$1(node) ? node.value : node;
		else return isCollection$1(node) ? node.getIn(rest, keepScalar) : void 0;
	}
	hasAllNullValues(allowScalar) {
		return this.items.every((node) => {
			if (!isPair(node)) return false;
			const n = node.value;
			return n == null || allowScalar && isScalar$1(n) && n.value == null && !n.commentBefore && !n.comment && !n.tag;
		});
	}
	/**
	* Checks if the collection includes a value with the key `key`.
	*/
	hasIn(path) {
		const [key, ...rest] = path;
		if (rest.length === 0) return this.has(key);
		const node = this.get(key, true);
		return isCollection$1(node) ? node.hasIn(rest) : false;
	}
	/**
	* Sets a value in this collection. For `!!set`, `value` needs to be a
	* boolean to add/remove the item from the set.
	*/
	setIn(path, value) {
		const [key, ...rest] = path;
		if (rest.length === 0) this.set(key, value);
		else {
			const node = this.get(key, true);
			if (isCollection$1(node)) node.setIn(rest, value);
			else if (node === void 0 && this.schema) this.set(key, collectionFromPath(this.schema, rest, value));
			else throw new Error(`Expected YAML collection at ${key}. Remaining path: ${rest}`);
		}
	}
};
//#endregion
//#region node_modules/yaml/browser/dist/stringify/stringifyComment.js
/**
* Stringifies a comment.
*
* Empty comment lines are left empty,
* lines consisting of a single space are replaced by `#`,
* and all other lines are prefixed with a `#`.
*/
const stringifyComment = (str) => str.replace(/^(?!$)(?: $)?/gm, "#");
function indentComment(comment, indent) {
	if (/^\n+$/.test(comment)) return comment.substring(1);
	return indent ? comment.replace(/^(?! *$)/gm, indent) : comment;
}
const lineComment = (str, indent, comment) => str.endsWith("\n") ? indentComment(comment, indent) : comment.includes("\n") ? "\n" + indentComment(comment, indent) : (str.endsWith(" ") ? "" : " ") + comment;
//#endregion
//#region node_modules/yaml/browser/dist/stringify/foldFlowLines.js
const FOLD_FLOW = "flow";
const FOLD_BLOCK = "block";
const FOLD_QUOTED = "quoted";
/**
* Tries to keep input at up to `lineWidth` characters, splitting only on spaces
* not followed by newlines or spaces unless `mode` is `'quoted'`. Lines are
* terminated with `\n` and started with `indent`.
*/
function foldFlowLines(text, indent, mode = "flow", { indentAtStart, lineWidth = 80, minContentWidth = 20, onFold, onOverflow } = {}) {
	if (!lineWidth || lineWidth < 0) return text;
	if (lineWidth < minContentWidth) minContentWidth = 0;
	const endStep = Math.max(1 + minContentWidth, 1 + lineWidth - indent.length);
	if (text.length <= endStep) return text;
	const folds = [];
	const escapedFolds = {};
	let end = lineWidth - indent.length;
	if (typeof indentAtStart === "number") if (indentAtStart > lineWidth - Math.max(2, minContentWidth)) folds.push(0);
	else end = lineWidth - indentAtStart;
	let split = void 0;
	let prev = void 0;
	let overflow = false;
	let i = -1;
	let escStart = -1;
	let escEnd = -1;
	if (mode === "block") {
		i = consumeMoreIndentedLines(text, i, indent.length);
		if (i !== -1) end = i + endStep;
	}
	for (let ch; ch = text[i += 1];) {
		if (mode === "quoted" && ch === "\\") {
			escStart = i;
			switch (text[i + 1]) {
				case "x":
					i += 3;
					break;
				case "u":
					i += 5;
					break;
				case "U":
					i += 9;
					break;
				default: i += 1;
			}
			escEnd = i;
		}
		if (ch === "\n") {
			if (mode === "block") i = consumeMoreIndentedLines(text, i, indent.length);
			end = i + indent.length + endStep;
			split = void 0;
		} else {
			if (ch === " " && prev && prev !== " " && prev !== "\n" && prev !== "	") {
				const next = text[i + 1];
				if (next && next !== " " && next !== "\n" && next !== "	") split = i;
			}
			if (i >= end) if (split) {
				folds.push(split);
				end = split + endStep;
				split = void 0;
			} else if (mode === "quoted") {
				while (prev === " " || prev === "	") {
					prev = ch;
					ch = text[i += 1];
					overflow = true;
				}
				const j = i > escEnd + 1 ? i - 2 : escStart - 1;
				if (escapedFolds[j]) return text;
				folds.push(j);
				escapedFolds[j] = true;
				end = j + endStep;
				split = void 0;
			} else overflow = true;
		}
		prev = ch;
	}
	if (overflow && onOverflow) onOverflow();
	if (folds.length === 0) return text;
	if (onFold) onFold();
	let res = text.slice(0, folds[0]);
	for (let i = 0; i < folds.length; ++i) {
		const fold = folds[i];
		const end = folds[i + 1] || text.length;
		if (fold === 0) res = `\n${indent}${text.slice(0, end)}`;
		else {
			if (mode === "quoted" && escapedFolds[fold]) res += `${text[fold]}\\`;
			res += `\n${indent}${text.slice(fold + 1, end)}`;
		}
	}
	return res;
}
/**
* Presumes `i + 1` is at the start of a line
* @returns index of last newline in more-indented block
*/
function consumeMoreIndentedLines(text, i, indent) {
	let end = i;
	let start = i + 1;
	let ch = text[start];
	while (ch === " " || ch === "	") if (i < start + indent) ch = text[++i];
	else {
		do
			ch = text[++i];
		while (ch && ch !== "\n");
		end = i;
		start = i + 1;
		ch = text[start];
	}
	return end;
}
//#endregion
//#region node_modules/yaml/browser/dist/stringify/stringifyString.js
const getFoldOptions = (ctx, isBlock) => ({
	indentAtStart: isBlock ? ctx.indent.length : ctx.indentAtStart,
	lineWidth: ctx.options.lineWidth,
	minContentWidth: ctx.options.minContentWidth
});
const containsDocumentMarker = (str) => /^(%|---|\.\.\.)/m.test(str);
function lineLengthOverLimit(str, lineWidth, indentLength) {
	if (!lineWidth || lineWidth < 0) return false;
	const limit = lineWidth - indentLength;
	const strLen = str.length;
	if (strLen <= limit) return false;
	for (let i = 0, start = 0; i < strLen; ++i) if (str[i] === "\n") {
		if (i - start > limit) return true;
		start = i + 1;
		if (strLen - start <= limit) return false;
	}
	return true;
}
function doubleQuotedString(value, ctx) {
	const json = JSON.stringify(value);
	if (ctx.options.doubleQuotedAsJSON) return json;
	const { implicitKey } = ctx;
	const minMultiLineLength = ctx.options.doubleQuotedMinMultiLineLength;
	const indent = ctx.indent || (containsDocumentMarker(value) ? "  " : "");
	let str = "";
	let start = 0;
	for (let i = 0, ch = json[i]; ch; ch = json[++i]) {
		if (ch === " " && json[i + 1] === "\\" && json[i + 2] === "n") {
			str += json.slice(start, i) + "\\ ";
			i += 1;
			start = i;
			ch = "\\";
		}
		if (ch === "\\") switch (json[i + 1]) {
			case "u":
				{
					str += json.slice(start, i);
					const code = json.substr(i + 2, 4);
					switch (code) {
						case "0000":
							str += "\\0";
							break;
						case "0007":
							str += "\\a";
							break;
						case "000b":
							str += "\\v";
							break;
						case "001b":
							str += "\\e";
							break;
						case "0085":
							str += "\\N";
							break;
						case "00a0":
							str += "\\_";
							break;
						case "2028":
							str += "\\L";
							break;
						case "2029":
							str += "\\P";
							break;
						default: if (code.substr(0, 2) === "00") str += "\\x" + code.substr(2);
						else str += json.substr(i, 6);
					}
					i += 5;
					start = i + 1;
				}
				break;
			case "n":
				if (implicitKey || json[i + 2] === "\"" || json.length < minMultiLineLength) i += 1;
				else {
					str += json.slice(start, i) + "\n\n";
					while (json[i + 2] === "\\" && json[i + 3] === "n" && json[i + 4] !== "\"") {
						str += "\n";
						i += 2;
					}
					str += indent;
					if (json[i + 2] === " ") str += "\\";
					i += 1;
					start = i + 1;
				}
				break;
			default: i += 1;
		}
	}
	str = start ? str + json.slice(start) : json;
	return implicitKey ? str : foldFlowLines(str, indent, FOLD_QUOTED, getFoldOptions(ctx, false));
}
function singleQuotedString(value, ctx) {
	if (ctx.options.singleQuote === false || ctx.implicitKey && value.includes("\n") || /[ \t]\n|\n[ \t]/.test(value)) return doubleQuotedString(value, ctx);
	const indent = ctx.indent || (containsDocumentMarker(value) ? "  " : "");
	const res = "'" + value.replace(/'/g, "''").replace(/\n+/g, `$&\n${indent}`) + "'";
	return ctx.implicitKey ? res : foldFlowLines(res, indent, FOLD_FLOW, getFoldOptions(ctx, false));
}
function quotedString(value, ctx) {
	const { singleQuote } = ctx.options;
	let qs;
	if (singleQuote === false) qs = doubleQuotedString;
	else {
		const hasDouble = value.includes("\"");
		const hasSingle = value.includes("'");
		if (hasDouble && !hasSingle) qs = singleQuotedString;
		else if (hasSingle && !hasDouble) qs = doubleQuotedString;
		else qs = singleQuote ? singleQuotedString : doubleQuotedString;
	}
	return qs(value, ctx);
}
let blockEndNewlines;
try {
	blockEndNewlines = /* @__PURE__ */ new RegExp("(^|(?<!\n))\n+(?!\n|$)", "g");
} catch {
	blockEndNewlines = /\n+(?!\n|$)/g;
}
function blockString({ comment, type, value }, ctx, onComment, onChompKeep) {
	const { blockQuote, commentString, lineWidth } = ctx.options;
	if (!blockQuote || /\n[\t ]+$/.test(value)) return quotedString(value, ctx);
	const indent = ctx.indent || (ctx.forceBlockIndent || containsDocumentMarker(value) ? "  " : "");
	const literal = blockQuote === "literal" ? true : blockQuote === "folded" || type === Scalar.BLOCK_FOLDED ? false : type === Scalar.BLOCK_LITERAL ? true : !lineLengthOverLimit(value, lineWidth, indent.length);
	if (!value) return literal ? "|\n" : ">\n";
	let chomp;
	let endStart;
	for (endStart = value.length; endStart > 0; --endStart) {
		const ch = value[endStart - 1];
		if (ch !== "\n" && ch !== "	" && ch !== " ") break;
	}
	let end = value.substring(endStart);
	const endNlPos = end.indexOf("\n");
	if (endNlPos === -1) chomp = "-";
	else if (value === end || endNlPos !== end.length - 1) {
		chomp = "+";
		if (onChompKeep) onChompKeep();
	} else chomp = "";
	if (end) {
		value = value.slice(0, -end.length);
		if (end[end.length - 1] === "\n") end = end.slice(0, -1);
		end = end.replace(blockEndNewlines, `$&${indent}`);
	}
	let startWithSpace = false;
	let startEnd;
	let startNlPos = -1;
	for (startEnd = 0; startEnd < value.length; ++startEnd) {
		const ch = value[startEnd];
		if (ch === " ") startWithSpace = true;
		else if (ch === "\n") startNlPos = startEnd;
		else break;
	}
	let start = value.substring(0, startNlPos < startEnd ? startNlPos + 1 : startEnd);
	if (start) {
		value = value.substring(start.length);
		start = start.replace(/\n+/g, `$&${indent}`);
	}
	let header = (startWithSpace ? indent ? "2" : "1" : "") + chomp;
	if (comment) {
		header += " " + commentString(comment.replace(/ ?[\r\n]+/g, " "));
		if (onComment) onComment();
	}
	if (!literal) {
		const foldedValue = value.replace(/\n+/g, "\n$&").replace(/(?:^|\n)([\t ].*)(?:([\n\t ]*)\n(?![\n\t ]))?/g, "$1$2").replace(/\n+/g, `$&${indent}`);
		let literalFallback = false;
		const foldOptions = getFoldOptions(ctx, true);
		if (blockQuote !== "folded" && type !== Scalar.BLOCK_FOLDED) foldOptions.onOverflow = () => {
			literalFallback = true;
		};
		const body = foldFlowLines(`${start}${foldedValue}${end}`, indent, FOLD_BLOCK, foldOptions);
		if (!literalFallback) return `>${header}\n${indent}${body}`;
	}
	value = value.replace(/\n+/g, `$&${indent}`);
	return `|${header}\n${indent}${start}${value}${end}`;
}
function plainString(item, ctx, onComment, onChompKeep) {
	const { type, value } = item;
	const { actualString, implicitKey, indent, indentStep, inFlow } = ctx;
	if (implicitKey && value.includes("\n") || inFlow && /[[\]{},]/.test(value)) return quotedString(value, ctx);
	if (/^[\n\t ,[\]{}#&*!|>'"%@`]|^[?-]$|^[?-][ \t]|[\n:][ \t]|[ \t]\n|[\n\t ]#|[\n\t :]$/.test(value)) return implicitKey || inFlow || !value.includes("\n") ? quotedString(value, ctx) : blockString(item, ctx, onComment, onChompKeep);
	if (!implicitKey && !inFlow && type !== Scalar.PLAIN && value.includes("\n")) return blockString(item, ctx, onComment, onChompKeep);
	if (containsDocumentMarker(value)) {
		if (indent === "") {
			ctx.forceBlockIndent = true;
			return blockString(item, ctx, onComment, onChompKeep);
		} else if (implicitKey && indent === indentStep) return quotedString(value, ctx);
	}
	const str = value.replace(/\n+/g, `$&\n${indent}`);
	if (actualString) {
		const test = (tag) => tag.default && tag.tag !== "tag:yaml.org,2002:str" && tag.test?.test(str);
		const { compat, tags } = ctx.doc.schema;
		if (tags.some(test) || compat?.some(test)) return quotedString(value, ctx);
	}
	return implicitKey ? str : foldFlowLines(str, indent, FOLD_FLOW, getFoldOptions(ctx, false));
}
function stringifyString(item, ctx, onComment, onChompKeep) {
	const { implicitKey, inFlow } = ctx;
	const ss = typeof item.value === "string" ? item : Object.assign({}, item, { value: String(item.value) });
	let { type } = item;
	if (type !== Scalar.QUOTE_DOUBLE) {
		if (/[\x00-\x08\x0b-\x1f\x7f-\x9f\u{D800}-\u{DFFF}]/u.test(ss.value)) type = Scalar.QUOTE_DOUBLE;
	}
	const _stringify = (_type) => {
		switch (_type) {
			case Scalar.BLOCK_FOLDED:
			case Scalar.BLOCK_LITERAL: return implicitKey || inFlow ? quotedString(ss.value, ctx) : blockString(ss, ctx, onComment, onChompKeep);
			case Scalar.QUOTE_DOUBLE: return doubleQuotedString(ss.value, ctx);
			case Scalar.QUOTE_SINGLE: return singleQuotedString(ss.value, ctx);
			case Scalar.PLAIN: return plainString(ss, ctx, onComment, onChompKeep);
			default: return null;
		}
	};
	let res = _stringify(type);
	if (res === null) {
		const { defaultKeyType, defaultStringType } = ctx.options;
		const t = implicitKey && defaultKeyType || defaultStringType;
		res = _stringify(t);
		if (res === null) throw new Error(`Unsupported default string type ${t}`);
	}
	return res;
}
//#endregion
//#region node_modules/yaml/browser/dist/stringify/stringify.js
function createStringifyContext(doc, options) {
	const opt = Object.assign({
		blockQuote: true,
		commentString: stringifyComment,
		defaultKeyType: null,
		defaultStringType: "PLAIN",
		directives: null,
		doubleQuotedAsJSON: false,
		doubleQuotedMinMultiLineLength: 40,
		falseStr: "false",
		flowCollectionPadding: true,
		indentSeq: true,
		lineWidth: 80,
		minContentWidth: 20,
		nullStr: "null",
		simpleKeys: false,
		singleQuote: null,
		trailingComma: false,
		trueStr: "true",
		verifyAliasOrder: true
	}, doc.schema.toStringOptions, options);
	let inFlow;
	switch (opt.collectionStyle) {
		case "block":
			inFlow = false;
			break;
		case "flow":
			inFlow = true;
			break;
		default: inFlow = null;
	}
	return {
		anchors: /* @__PURE__ */ new Set(),
		doc,
		flowCollectionPadding: opt.flowCollectionPadding ? " " : "",
		indent: "",
		indentStep: typeof opt.indent === "number" ? " ".repeat(opt.indent) : "  ",
		inFlow,
		options: opt
	};
}
function getTagObject(tags, item) {
	if (item.tag) {
		const match = tags.filter((t) => t.tag === item.tag);
		if (match.length > 0) return match.find((t) => t.format === item.format) ?? match[0];
	}
	let tagObj = void 0;
	let obj;
	if (isScalar$1(item)) {
		obj = item.value;
		let match = tags.filter((t) => t.identify?.(obj));
		if (match.length > 1) {
			const testMatch = match.filter((t) => t.test);
			if (testMatch.length > 0) match = testMatch;
		}
		tagObj = match.find((t) => t.format === item.format) ?? match.find((t) => !t.format);
	} else {
		obj = item;
		tagObj = tags.find((t) => t.nodeClass && obj instanceof t.nodeClass);
	}
	if (!tagObj) {
		const name = obj?.constructor?.name ?? (obj === null ? "null" : typeof obj);
		throw new Error(`Tag not resolved for ${name} value`);
	}
	return tagObj;
}
function stringifyProps(node, tagObj, { anchors, doc }) {
	if (!doc.directives) return "";
	const props = [];
	const anchor = (isScalar$1(node) || isCollection$1(node)) && node.anchor;
	if (anchor && anchorIsValid(anchor)) {
		anchors.add(anchor);
		props.push(`&${anchor}`);
	}
	const tag = node.tag ?? (tagObj.default ? null : tagObj.tag);
	if (tag) props.push(doc.directives.tagString(tag));
	return props.join(" ");
}
function stringify$2(item, ctx, onComment, onChompKeep) {
	if (isPair(item)) return item.toString(ctx, onComment, onChompKeep);
	if (isAlias(item)) {
		if (ctx.doc.directives) return item.toString(ctx);
		if (ctx.resolvedAliases?.has(item)) throw new TypeError(`Cannot stringify circular structure without alias nodes`);
		else {
			if (ctx.resolvedAliases) ctx.resolvedAliases.add(item);
			else ctx.resolvedAliases = new Set([item]);
			item = item.resolve(ctx.doc);
		}
	}
	let tagObj = void 0;
	const node = isNode(item) ? item : ctx.doc.createNode(item, { onTagObj: (o) => tagObj = o });
	tagObj ?? (tagObj = getTagObject(ctx.doc.schema.tags, node));
	const props = stringifyProps(node, tagObj, ctx);
	if (props.length > 0) ctx.indentAtStart = (ctx.indentAtStart ?? 0) + props.length + 1;
	const str = typeof tagObj.stringify === "function" ? tagObj.stringify(node, ctx, onComment, onChompKeep) : isScalar$1(node) ? stringifyString(node, ctx, onComment, onChompKeep) : node.toString(ctx, onComment, onChompKeep);
	if (!props) return str;
	return isScalar$1(node) || str[0] === "{" || str[0] === "[" ? `${props} ${str}` : `${props}\n${ctx.indent}${str}`;
}
//#endregion
//#region node_modules/yaml/browser/dist/stringify/stringifyPair.js
function stringifyPair({ key, value }, ctx, onComment, onChompKeep) {
	const { allNullValues, doc, indent, indentStep, options: { commentString, indentSeq, simpleKeys } } = ctx;
	let keyComment = isNode(key) && key.comment || null;
	if (simpleKeys) {
		if (keyComment) throw new Error("With simple keys, key nodes cannot have comments");
		if (isCollection$1(key) || !isNode(key) && typeof key === "object") throw new Error("With simple keys, collection cannot be used as a key value");
	}
	let explicitKey = !simpleKeys && (!key || keyComment && value == null && !ctx.inFlow || isCollection$1(key) || (isScalar$1(key) ? key.type === Scalar.BLOCK_FOLDED || key.type === Scalar.BLOCK_LITERAL : typeof key === "object"));
	ctx = Object.assign({}, ctx, {
		allNullValues: false,
		implicitKey: !explicitKey && (simpleKeys || !allNullValues),
		indent: indent + indentStep
	});
	let keyCommentDone = false;
	let chompKeep = false;
	let str = stringify$2(key, ctx, () => keyCommentDone = true, () => chompKeep = true);
	if (!explicitKey && !ctx.inFlow && str.length > 1024) {
		if (simpleKeys) throw new Error("With simple keys, single line scalar must not span more than 1024 characters");
		explicitKey = true;
	}
	if (ctx.inFlow) {
		if (allNullValues || value == null) {
			if (keyCommentDone && onComment) onComment();
			return str === "" ? "?" : explicitKey ? `? ${str}` : str;
		}
	} else if (allNullValues && !simpleKeys || value == null && explicitKey) {
		str = `? ${str}`;
		if (keyComment && !keyCommentDone) str += lineComment(str, ctx.indent, commentString(keyComment));
		else if (chompKeep && onChompKeep) onChompKeep();
		return str;
	}
	if (keyCommentDone) keyComment = null;
	if (explicitKey) {
		if (keyComment) str += lineComment(str, ctx.indent, commentString(keyComment));
		str = `? ${str}\n${indent}:`;
	} else {
		str = `${str}:`;
		if (keyComment) str += lineComment(str, ctx.indent, commentString(keyComment));
	}
	let vsb, vcb, valueComment;
	if (isNode(value)) {
		vsb = !!value.spaceBefore;
		vcb = value.commentBefore;
		valueComment = value.comment;
	} else {
		vsb = false;
		vcb = null;
		valueComment = null;
		if (value && typeof value === "object") value = doc.createNode(value);
	}
	ctx.implicitKey = false;
	if (!explicitKey && !keyComment && isScalar$1(value)) ctx.indentAtStart = str.length + 1;
	chompKeep = false;
	if (!indentSeq && indentStep.length >= 2 && !ctx.inFlow && !explicitKey && isSeq(value) && !value.flow && !value.tag && !value.anchor) ctx.indent = ctx.indent.substring(2);
	let valueCommentDone = false;
	const valueStr = stringify$2(value, ctx, () => valueCommentDone = true, () => chompKeep = true);
	let ws = " ";
	if (keyComment || vsb || vcb) {
		ws = vsb ? "\n" : "";
		if (vcb) {
			const cs = commentString(vcb);
			ws += `\n${indentComment(cs, ctx.indent)}`;
		}
		if (valueStr === "" && !ctx.inFlow) {
			if (ws === "\n" && valueComment) ws = "\n\n";
		} else ws += `\n${ctx.indent}`;
	} else if (!explicitKey && isCollection$1(value)) {
		const vs0 = valueStr[0];
		const nl0 = valueStr.indexOf("\n");
		const hasNewline = nl0 !== -1;
		const flow = ctx.inFlow ?? value.flow ?? value.items.length === 0;
		if (hasNewline || !flow) {
			let hasPropsLine = false;
			if (hasNewline && (vs0 === "&" || vs0 === "!")) {
				let sp0 = valueStr.indexOf(" ");
				if (vs0 === "&" && sp0 !== -1 && sp0 < nl0 && valueStr[sp0 + 1] === "!") sp0 = valueStr.indexOf(" ", sp0 + 1);
				if (sp0 === -1 || nl0 < sp0) hasPropsLine = true;
			}
			if (!hasPropsLine) ws = `\n${ctx.indent}`;
		}
	} else if (valueStr === "" || valueStr[0] === "\n") ws = "";
	str += ws + valueStr;
	if (ctx.inFlow) {
		if (valueCommentDone && onComment) onComment();
	} else if (valueComment && !valueCommentDone) str += lineComment(str, ctx.indent, commentString(valueComment));
	else if (chompKeep && onChompKeep) onChompKeep();
	return str;
}
//#endregion
//#region node_modules/yaml/browser/dist/log.js
function warn(logLevel, warning) {
	if (logLevel === "debug" || logLevel === "warn") console.warn(warning);
}
//#endregion
//#region node_modules/yaml/browser/dist/schema/yaml-1.1/merge.js
const MERGE_KEY = "<<";
const merge = {
	identify: (value) => value === MERGE_KEY || typeof value === "symbol" && value.description === MERGE_KEY,
	default: "key",
	tag: "tag:yaml.org,2002:merge",
	test: /^<<$/,
	resolve: () => Object.assign(new Scalar(Symbol(MERGE_KEY)), { addToJSMap: addMergeToJSMap }),
	stringify: () => MERGE_KEY
};
const isMergeKey = (ctx, key) => (merge.identify(key) || isScalar$1(key) && (!key.type || key.type === Scalar.PLAIN) && merge.identify(key.value)) && ctx?.doc.schema.tags.some((tag) => tag.tag === merge.tag && tag.default);
function addMergeToJSMap(ctx, map, value) {
	const source = resolveAliasValue(ctx, value);
	if (isSeq(source)) for (const it of source.items) mergeValue(ctx, map, it);
	else if (Array.isArray(source)) for (const it of source) mergeValue(ctx, map, it);
	else mergeValue(ctx, map, source);
}
function mergeValue(ctx, map, value) {
	const source = resolveAliasValue(ctx, value);
	if (!isMap(source)) throw new Error("Merge sources must be maps or map aliases");
	const srcMap = source.toJSON(null, ctx, Map);
	for (const [key, value] of srcMap) if (map instanceof Map) {
		if (!map.has(key)) map.set(key, value);
	} else if (map instanceof Set) map.add(key);
	else if (!Object.prototype.hasOwnProperty.call(map, key)) Object.defineProperty(map, key, {
		value,
		writable: true,
		enumerable: true,
		configurable: true
	});
	return map;
}
function resolveAliasValue(ctx, value) {
	return ctx && isAlias(value) ? value.resolve(ctx.doc, ctx) : value;
}
//#endregion
//#region node_modules/yaml/browser/dist/nodes/addPairToJSMap.js
function addPairToJSMap(ctx, map, { key, value }) {
	if (isNode(key) && key.addToJSMap) key.addToJSMap(ctx, map, value);
	else if (isMergeKey(ctx, key)) addMergeToJSMap(ctx, map, value);
	else {
		const jsKey = toJS(key, "", ctx);
		if (map instanceof Map) map.set(jsKey, toJS(value, jsKey, ctx));
		else if (map instanceof Set) map.add(jsKey);
		else {
			const stringKey = stringifyKey(key, jsKey, ctx);
			const jsValue = toJS(value, stringKey, ctx);
			if (stringKey in map) Object.defineProperty(map, stringKey, {
				value: jsValue,
				writable: true,
				enumerable: true,
				configurable: true
			});
			else map[stringKey] = jsValue;
		}
	}
	return map;
}
function stringifyKey(key, jsKey, ctx) {
	if (jsKey === null) return "";
	if (typeof jsKey !== "object") return String(jsKey);
	if (isNode(key) && ctx?.doc) {
		const strCtx = createStringifyContext(ctx.doc, {});
		strCtx.anchors = /* @__PURE__ */ new Set();
		for (const node of ctx.anchors.keys()) strCtx.anchors.add(node.anchor);
		strCtx.inFlow = true;
		strCtx.inStringifyKey = true;
		const strKey = key.toString(strCtx);
		if (!ctx.mapKeyWarned) {
			let jsonStr = JSON.stringify(strKey);
			if (jsonStr.length > 40) jsonStr = jsonStr.substring(0, 36) + "...\"";
			warn(ctx.doc.options.logLevel, `Keys with collection values will be stringified due to JS Object restrictions: ${jsonStr}. Set mapAsMap: true to use object keys.`);
			ctx.mapKeyWarned = true;
		}
		return strKey;
	}
	return JSON.stringify(jsKey);
}
//#endregion
//#region node_modules/yaml/browser/dist/nodes/Pair.js
function createPair(key, value, ctx) {
	return new Pair(createNode(key, void 0, ctx), createNode(value, void 0, ctx));
}
var Pair = class Pair {
	constructor(key, value = null) {
		Object.defineProperty(this, NODE_TYPE, { value: PAIR });
		this.key = key;
		this.value = value;
	}
	clone(schema) {
		let { key, value } = this;
		if (isNode(key)) key = key.clone(schema);
		if (isNode(value)) value = value.clone(schema);
		return new Pair(key, value);
	}
	toJSON(_, ctx) {
		return addPairToJSMap(ctx, ctx?.mapAsMap ? /* @__PURE__ */ new Map() : {}, this);
	}
	toString(ctx, onComment, onChompKeep) {
		return ctx?.doc ? stringifyPair(this, ctx, onComment, onChompKeep) : JSON.stringify(this);
	}
};
//#endregion
//#region node_modules/yaml/browser/dist/stringify/stringifyCollection.js
function stringifyCollection(collection, ctx, options) {
	return (ctx.inFlow ?? collection.flow ? stringifyFlowCollection : stringifyBlockCollection)(collection, ctx, options);
}
function stringifyBlockCollection({ comment, items }, ctx, { blockItemPrefix, flowChars, itemIndent, onChompKeep, onComment }) {
	const { indent, options: { commentString } } = ctx;
	const itemCtx = Object.assign({}, ctx, {
		indent: itemIndent,
		type: null
	});
	let chompKeep = false;
	const lines = [];
	for (let i = 0; i < items.length; ++i) {
		const item = items[i];
		let comment = null;
		if (isNode(item)) {
			if (!chompKeep && item.spaceBefore) lines.push("");
			addCommentBefore(ctx, lines, item.commentBefore, chompKeep);
			if (item.comment) comment = item.comment;
		} else if (isPair(item)) {
			const ik = isNode(item.key) ? item.key : null;
			if (ik) {
				if (!chompKeep && ik.spaceBefore) lines.push("");
				addCommentBefore(ctx, lines, ik.commentBefore, chompKeep);
			}
		}
		chompKeep = false;
		let str = stringify$2(item, itemCtx, () => comment = null, () => chompKeep = true);
		if (comment) str += lineComment(str, itemIndent, commentString(comment));
		if (chompKeep && comment) chompKeep = false;
		lines.push(blockItemPrefix + str);
	}
	let str;
	if (lines.length === 0) str = flowChars.start + flowChars.end;
	else {
		str = lines[0];
		for (let i = 1; i < lines.length; ++i) {
			const line = lines[i];
			str += line ? `\n${indent}${line}` : "\n";
		}
	}
	if (comment) {
		str += "\n" + indentComment(commentString(comment), indent);
		if (onComment) onComment();
	} else if (chompKeep && onChompKeep) onChompKeep();
	return str;
}
function stringifyFlowCollection({ items }, ctx, { flowChars, itemIndent }) {
	const { indent, indentStep, flowCollectionPadding: fcPadding, options: { commentString } } = ctx;
	itemIndent += indentStep;
	const itemCtx = Object.assign({}, ctx, {
		indent: itemIndent,
		inFlow: true,
		type: null
	});
	let reqNewline = false;
	let linesAtValue = 0;
	const lines = [];
	for (let i = 0; i < items.length; ++i) {
		const item = items[i];
		let comment = null;
		if (isNode(item)) {
			if (item.spaceBefore) lines.push("");
			addCommentBefore(ctx, lines, item.commentBefore, false);
			if (item.comment) comment = item.comment;
		} else if (isPair(item)) {
			const ik = isNode(item.key) ? item.key : null;
			if (ik) {
				if (ik.spaceBefore) lines.push("");
				addCommentBefore(ctx, lines, ik.commentBefore, false);
				if (ik.comment) reqNewline = true;
			}
			const iv = isNode(item.value) ? item.value : null;
			if (iv) {
				if (iv.comment) comment = iv.comment;
				if (iv.commentBefore) reqNewline = true;
			} else if (item.value == null && ik?.comment) comment = ik.comment;
		}
		if (comment) reqNewline = true;
		let str = stringify$2(item, itemCtx, () => comment = null);
		reqNewline || (reqNewline = lines.length > linesAtValue || str.includes("\n"));
		if (i < items.length - 1) str += ",";
		else if (ctx.options.trailingComma) {
			if (ctx.options.lineWidth > 0) reqNewline || (reqNewline = lines.reduce((sum, line) => sum + line.length + 2, 2) + (str.length + 2) > ctx.options.lineWidth);
			if (reqNewline) str += ",";
		}
		if (comment) str += lineComment(str, itemIndent, commentString(comment));
		lines.push(str);
		linesAtValue = lines.length;
	}
	const { start, end } = flowChars;
	if (lines.length === 0) return start + end;
	else {
		if (!reqNewline) {
			const len = lines.reduce((sum, line) => sum + line.length + 2, 2);
			reqNewline = ctx.options.lineWidth > 0 && len > ctx.options.lineWidth;
		}
		if (reqNewline) {
			let str = start;
			for (const line of lines) str += line ? `\n${indentStep}${indent}${line}` : "\n";
			return `${str}\n${indent}${end}`;
		} else return `${start}${fcPadding}${lines.join(" ")}${fcPadding}${end}`;
	}
}
function addCommentBefore({ indent, options: { commentString } }, lines, comment, chompKeep) {
	if (comment && chompKeep) comment = comment.replace(/^\n+/, "");
	if (comment) {
		const ic = indentComment(commentString(comment), indent);
		lines.push(ic.trimStart());
	}
}
//#endregion
//#region node_modules/yaml/browser/dist/nodes/YAMLMap.js
function findPair(items, key) {
	const k = isScalar$1(key) ? key.value : key;
	for (const it of items) if (isPair(it)) {
		if (it.key === key || it.key === k) return it;
		if (isScalar$1(it.key) && it.key.value === k) return it;
	}
}
var YAMLMap = class extends Collection {
	static get tagName() {
		return "tag:yaml.org,2002:map";
	}
	constructor(schema) {
		super(MAP, schema);
		this.items = [];
	}
	/**
	* A generic collection parsing method that can be extended
	* to other node classes that inherit from YAMLMap
	*/
	static from(schema, obj, ctx) {
		const { keepUndefined, replacer } = ctx;
		const map = new this(schema);
		const add = (key, value) => {
			if (typeof replacer === "function") value = replacer.call(obj, key, value);
			else if (Array.isArray(replacer) && !replacer.includes(key)) return;
			if (value !== void 0 || keepUndefined) map.items.push(createPair(key, value, ctx));
		};
		if (obj instanceof Map) for (const [key, value] of obj) add(key, value);
		else if (obj && typeof obj === "object") for (const key of Object.keys(obj)) add(key, obj[key]);
		if (typeof schema.sortMapEntries === "function") map.items.sort(schema.sortMapEntries);
		return map;
	}
	/**
	* Adds a value to the collection.
	*
	* @param overwrite - If not set `true`, using a key that is already in the
	*   collection will throw. Otherwise, overwrites the previous value.
	*/
	add(pair, overwrite) {
		let _pair;
		if (isPair(pair)) _pair = pair;
		else if (!pair || typeof pair !== "object" || !("key" in pair)) _pair = new Pair(pair, pair?.value);
		else _pair = new Pair(pair.key, pair.value);
		const prev = findPair(this.items, _pair.key);
		const sortEntries = this.schema?.sortMapEntries;
		if (prev) {
			if (!overwrite) throw new Error(`Key ${_pair.key} already set`);
			if (isScalar$1(prev.value) && isScalarValue(_pair.value)) prev.value.value = _pair.value;
			else prev.value = _pair.value;
		} else if (sortEntries) {
			const i = this.items.findIndex((item) => sortEntries(_pair, item) < 0);
			if (i === -1) this.items.push(_pair);
			else this.items.splice(i, 0, _pair);
		} else this.items.push(_pair);
	}
	delete(key) {
		const it = findPair(this.items, key);
		if (!it) return false;
		return this.items.splice(this.items.indexOf(it), 1).length > 0;
	}
	get(key, keepScalar) {
		const node = findPair(this.items, key)?.value;
		return (!keepScalar && isScalar$1(node) ? node.value : node) ?? void 0;
	}
	has(key) {
		return !!findPair(this.items, key);
	}
	set(key, value) {
		this.add(new Pair(key, value), true);
	}
	/**
	* @param ctx - Conversion context, originally set in Document#toJS()
	* @param {Class} Type - If set, forces the returned collection type
	* @returns Instance of Type, Map, or Object
	*/
	toJSON(_, ctx, Type) {
		const map = Type ? new Type() : ctx?.mapAsMap ? /* @__PURE__ */ new Map() : {};
		if (ctx?.onCreate) ctx.onCreate(map);
		for (const item of this.items) addPairToJSMap(ctx, map, item);
		return map;
	}
	toString(ctx, onComment, onChompKeep) {
		if (!ctx) return JSON.stringify(this);
		for (const item of this.items) if (!isPair(item)) throw new Error(`Map items must all be pairs; found ${JSON.stringify(item)} instead`);
		if (!ctx.allNullValues && this.hasAllNullValues(false)) ctx = Object.assign({}, ctx, { allNullValues: true });
		return stringifyCollection(this, ctx, {
			blockItemPrefix: "",
			flowChars: {
				start: "{",
				end: "}"
			},
			itemIndent: ctx.indent || "",
			onChompKeep,
			onComment
		});
	}
};
//#endregion
//#region node_modules/yaml/browser/dist/schema/common/map.js
const map = {
	collection: "map",
	default: true,
	nodeClass: YAMLMap,
	tag: "tag:yaml.org,2002:map",
	resolve(map, onError) {
		if (!isMap(map)) onError("Expected a mapping for this tag");
		return map;
	},
	createNode: (schema, obj, ctx) => YAMLMap.from(schema, obj, ctx)
};
//#endregion
//#region node_modules/yaml/browser/dist/nodes/YAMLSeq.js
var YAMLSeq = class extends Collection {
	static get tagName() {
		return "tag:yaml.org,2002:seq";
	}
	constructor(schema) {
		super(SEQ, schema);
		this.items = [];
	}
	add(value) {
		this.items.push(value);
	}
	/**
	* Removes a value from the collection.
	*
	* `key` must contain a representation of an integer for this to succeed.
	* It may be wrapped in a `Scalar`.
	*
	* @returns `true` if the item was found and removed.
	*/
	delete(key) {
		const idx = asItemIndex(key);
		if (typeof idx !== "number") return false;
		return this.items.splice(idx, 1).length > 0;
	}
	get(key, keepScalar) {
		const idx = asItemIndex(key);
		if (typeof idx !== "number") return void 0;
		const it = this.items[idx];
		return !keepScalar && isScalar$1(it) ? it.value : it;
	}
	/**
	* Checks if the collection includes a value with the key `key`.
	*
	* `key` must contain a representation of an integer for this to succeed.
	* It may be wrapped in a `Scalar`.
	*/
	has(key) {
		const idx = asItemIndex(key);
		return typeof idx === "number" && idx < this.items.length;
	}
	/**
	* Sets a value in this collection. For `!!set`, `value` needs to be a
	* boolean to add/remove the item from the set.
	*
	* If `key` does not contain a representation of an integer, this will throw.
	* It may be wrapped in a `Scalar`.
	*/
	set(key, value) {
		const idx = asItemIndex(key);
		if (typeof idx !== "number") throw new Error(`Expected a valid index, not ${key}.`);
		const prev = this.items[idx];
		if (isScalar$1(prev) && isScalarValue(value)) prev.value = value;
		else this.items[idx] = value;
	}
	toJSON(_, ctx) {
		const seq = [];
		if (ctx?.onCreate) ctx.onCreate(seq);
		let i = 0;
		for (const item of this.items) seq.push(toJS(item, String(i++), ctx));
		return seq;
	}
	toString(ctx, onComment, onChompKeep) {
		if (!ctx) return JSON.stringify(this);
		return stringifyCollection(this, ctx, {
			blockItemPrefix: "- ",
			flowChars: {
				start: "[",
				end: "]"
			},
			itemIndent: (ctx.indent || "") + "  ",
			onChompKeep,
			onComment
		});
	}
	static from(schema, obj, ctx) {
		const { replacer } = ctx;
		const seq = new this(schema);
		if (obj && Symbol.iterator in Object(obj)) {
			let i = 0;
			for (let it of obj) {
				if (typeof replacer === "function") {
					const key = obj instanceof Set ? it : String(i++);
					it = replacer.call(obj, key, it);
				}
				seq.items.push(createNode(it, void 0, ctx));
			}
		}
		return seq;
	}
};
function asItemIndex(key) {
	let idx = isScalar$1(key) ? key.value : key;
	if (idx && typeof idx === "string") idx = Number(idx);
	return typeof idx === "number" && Number.isInteger(idx) && idx >= 0 ? idx : null;
}
//#endregion
//#region node_modules/yaml/browser/dist/schema/common/seq.js
const seq = {
	collection: "seq",
	default: true,
	nodeClass: YAMLSeq,
	tag: "tag:yaml.org,2002:seq",
	resolve(seq, onError) {
		if (!isSeq(seq)) onError("Expected a sequence for this tag");
		return seq;
	},
	createNode: (schema, obj, ctx) => YAMLSeq.from(schema, obj, ctx)
};
//#endregion
//#region node_modules/yaml/browser/dist/schema/common/string.js
const string = {
	identify: (value) => typeof value === "string",
	default: true,
	tag: "tag:yaml.org,2002:str",
	resolve: (str) => str,
	stringify(item, ctx, onComment, onChompKeep) {
		ctx = Object.assign({ actualString: true }, ctx);
		return stringifyString(item, ctx, onComment, onChompKeep);
	}
};
//#endregion
//#region node_modules/yaml/browser/dist/schema/common/null.js
const nullTag = {
	identify: (value) => value == null,
	createNode: () => new Scalar(null),
	default: true,
	tag: "tag:yaml.org,2002:null",
	test: /^(?:~|[Nn]ull|NULL)?$/,
	resolve: () => new Scalar(null),
	stringify: ({ source }, ctx) => typeof source === "string" && nullTag.test.test(source) ? source : ctx.options.nullStr
};
//#endregion
//#region node_modules/yaml/browser/dist/schema/core/bool.js
const boolTag = {
	identify: (value) => typeof value === "boolean",
	default: true,
	tag: "tag:yaml.org,2002:bool",
	test: /^(?:[Tt]rue|TRUE|[Ff]alse|FALSE)$/,
	resolve: (str) => new Scalar(str[0] === "t" || str[0] === "T"),
	stringify({ source, value }, ctx) {
		if (source && boolTag.test.test(source)) {
			if (value === (source[0] === "t" || source[0] === "T")) return source;
		}
		return value ? ctx.options.trueStr : ctx.options.falseStr;
	}
};
//#endregion
//#region node_modules/yaml/browser/dist/stringify/stringifyNumber.js
function stringifyNumber({ format, minFractionDigits, tag, value }) {
	if (typeof value === "bigint") return String(value);
	const num = typeof value === "number" ? value : Number(value);
	if (!isFinite(num)) return isNaN(num) ? ".nan" : num < 0 ? "-.inf" : ".inf";
	let n = Object.is(value, -0) ? "-0" : JSON.stringify(value);
	if (!format && minFractionDigits && (!tag || tag === "tag:yaml.org,2002:float") && /^-?\d/.test(n) && !n.includes("e")) {
		let i = n.indexOf(".");
		if (i < 0) {
			i = n.length;
			n += ".";
		}
		let d = minFractionDigits - (n.length - i - 1);
		while (d-- > 0) n += "0";
	}
	return n;
}
//#endregion
//#region node_modules/yaml/browser/dist/schema/core/float.js
const floatNaN$1 = {
	identify: (value) => typeof value === "number",
	default: true,
	tag: "tag:yaml.org,2002:float",
	test: /^(?:[-+]?\.(?:inf|Inf|INF)|\.nan|\.NaN|\.NAN)$/,
	resolve: (str) => str.slice(-3).toLowerCase() === "nan" ? NaN : str[0] === "-" ? Number.NEGATIVE_INFINITY : Number.POSITIVE_INFINITY,
	stringify: stringifyNumber
};
const floatExp$1 = {
	identify: (value) => typeof value === "number",
	default: true,
	tag: "tag:yaml.org,2002:float",
	format: "EXP",
	test: /^[-+]?(?:\.[0-9]+|[0-9]+(?:\.[0-9]*)?)[eE][-+]?[0-9]+$/,
	resolve: (str) => parseFloat(str),
	stringify(node) {
		const num = Number(node.value);
		return isFinite(num) ? num.toExponential() : stringifyNumber(node);
	}
};
const float$1 = {
	identify: (value) => typeof value === "number",
	default: true,
	tag: "tag:yaml.org,2002:float",
	test: /^[-+]?(?:\.[0-9]+|[0-9]+\.[0-9]*)$/,
	resolve(str) {
		const node = new Scalar(parseFloat(str));
		const dot = str.indexOf(".");
		if (dot !== -1 && str[str.length - 1] === "0") node.minFractionDigits = str.length - dot - 1;
		return node;
	},
	stringify: stringifyNumber
};
//#endregion
//#region node_modules/yaml/browser/dist/schema/core/int.js
const intIdentify$2 = (value) => typeof value === "bigint" || Number.isInteger(value);
const intResolve$1 = (str, offset, radix, { intAsBigInt }) => intAsBigInt ? BigInt(str) : parseInt(str.substring(offset), radix);
function intStringify$1(node, radix, prefix) {
	const { value } = node;
	if (intIdentify$2(value) && value >= 0) return prefix + value.toString(radix);
	return stringifyNumber(node);
}
const intOct$1 = {
	identify: (value) => intIdentify$2(value) && value >= 0,
	default: true,
	tag: "tag:yaml.org,2002:int",
	format: "OCT",
	test: /^0o[0-7]+$/,
	resolve: (str, _onError, opt) => intResolve$1(str, 2, 8, opt),
	stringify: (node) => intStringify$1(node, 8, "0o")
};
const int$1 = {
	identify: intIdentify$2,
	default: true,
	tag: "tag:yaml.org,2002:int",
	test: /^[-+]?[0-9]+$/,
	resolve: (str, _onError, opt) => intResolve$1(str, 0, 10, opt),
	stringify: stringifyNumber
};
const intHex$1 = {
	identify: (value) => intIdentify$2(value) && value >= 0,
	default: true,
	tag: "tag:yaml.org,2002:int",
	format: "HEX",
	test: /^0x[0-9a-fA-F]+$/,
	resolve: (str, _onError, opt) => intResolve$1(str, 2, 16, opt),
	stringify: (node) => intStringify$1(node, 16, "0x")
};
//#endregion
//#region node_modules/yaml/browser/dist/schema/core/schema.js
const schema$2 = [
	map,
	seq,
	string,
	nullTag,
	boolTag,
	intOct$1,
	int$1,
	intHex$1,
	floatNaN$1,
	floatExp$1,
	float$1
];
//#endregion
//#region node_modules/yaml/browser/dist/schema/json/schema.js
function intIdentify$1(value) {
	return typeof value === "bigint" || Number.isInteger(value);
}
const stringifyJSON = ({ value }) => JSON.stringify(value);
const jsonScalars = [
	{
		identify: (value) => typeof value === "string",
		default: true,
		tag: "tag:yaml.org,2002:str",
		resolve: (str) => str,
		stringify: stringifyJSON
	},
	{
		identify: (value) => value == null,
		createNode: () => new Scalar(null),
		default: true,
		tag: "tag:yaml.org,2002:null",
		test: /^null$/,
		resolve: () => null,
		stringify: stringifyJSON
	},
	{
		identify: (value) => typeof value === "boolean",
		default: true,
		tag: "tag:yaml.org,2002:bool",
		test: /^true$|^false$/,
		resolve: (str) => str === "true",
		stringify: stringifyJSON
	},
	{
		identify: intIdentify$1,
		default: true,
		tag: "tag:yaml.org,2002:int",
		test: /^-?(?:0|[1-9][0-9]*)$/,
		resolve: (str, _onError, { intAsBigInt }) => intAsBigInt ? BigInt(str) : parseInt(str, 10),
		stringify: ({ value }) => intIdentify$1(value) ? value.toString() : JSON.stringify(value)
	},
	{
		identify: (value) => typeof value === "number",
		default: true,
		tag: "tag:yaml.org,2002:float",
		test: /^-?(?:0|[1-9][0-9]*)(?:\.[0-9]*)?(?:[eE][-+]?[0-9]+)?$/,
		resolve: (str) => parseFloat(str),
		stringify: stringifyJSON
	}
];
const schema$1 = [map, seq].concat(jsonScalars, {
	default: true,
	tag: "",
	test: /^/,
	resolve(str, onError) {
		onError(`Unresolved plain scalar ${JSON.stringify(str)}`);
		return str;
	}
});
//#endregion
//#region node_modules/yaml/browser/dist/schema/yaml-1.1/binary.js
const binary = {
	identify: (value) => value instanceof Uint8Array,
	default: false,
	tag: "tag:yaml.org,2002:binary",
	/**
	* Returns a Buffer in node and an Uint8Array in browsers
	*
	* To use the resulting buffer as an image, you'll want to do something like:
	*
	*   const blob = new Blob([buffer], { type: 'image/jpeg' })
	*   document.querySelector('#photo').src = URL.createObjectURL(blob)
	*/
	resolve(src, onError) {
		if (typeof atob === "function") {
			const str = atob(src.replace(/[\n\r]/g, ""));
			const buffer = new Uint8Array(str.length);
			for (let i = 0; i < str.length; ++i) buffer[i] = str.charCodeAt(i);
			return buffer;
		} else {
			onError("This environment does not support reading binary tags; either Buffer or atob is required");
			return src;
		}
	},
	stringify({ comment, type, value }, ctx, onComment, onChompKeep) {
		if (!value) return "";
		const buf = value;
		let str;
		if (typeof btoa === "function") {
			let s = "";
			for (let i = 0; i < buf.length; ++i) s += String.fromCharCode(buf[i]);
			str = btoa(s);
		} else throw new Error("This environment does not support writing binary tags; either Buffer or btoa is required");
		type ?? (type = Scalar.BLOCK_LITERAL);
		if (type !== Scalar.QUOTE_DOUBLE) {
			const lineWidth = Math.max(ctx.options.lineWidth - ctx.indent.length, ctx.options.minContentWidth);
			const n = Math.ceil(str.length / lineWidth);
			const lines = new Array(n);
			for (let i = 0, o = 0; i < n; ++i, o += lineWidth) lines[i] = str.substr(o, lineWidth);
			str = lines.join(type === Scalar.BLOCK_LITERAL ? "\n" : " ");
		}
		return stringifyString({
			comment,
			type,
			value: str
		}, ctx, onComment, onChompKeep);
	}
};
//#endregion
//#region node_modules/yaml/browser/dist/schema/yaml-1.1/pairs.js
function resolvePairs(seq, onError) {
	if (isSeq(seq)) for (let i = 0; i < seq.items.length; ++i) {
		let item = seq.items[i];
		if (isPair(item)) continue;
		else if (isMap(item)) {
			if (item.items.length > 1) onError("Each pair must have its own sequence indicator");
			const pair = item.items[0] || new Pair(new Scalar(null));
			if (item.commentBefore) pair.key.commentBefore = pair.key.commentBefore ? `${item.commentBefore}\n${pair.key.commentBefore}` : item.commentBefore;
			if (item.comment) {
				const cn = pair.value ?? pair.key;
				cn.comment = cn.comment ? `${item.comment}\n${cn.comment}` : item.comment;
			}
			item = pair;
		}
		seq.items[i] = isPair(item) ? item : new Pair(item);
	}
	else onError("Expected a sequence for this tag");
	return seq;
}
function createPairs(schema, iterable, ctx) {
	const { replacer } = ctx;
	const pairs = new YAMLSeq(schema);
	pairs.tag = "tag:yaml.org,2002:pairs";
	let i = 0;
	if (iterable && Symbol.iterator in Object(iterable)) for (let it of iterable) {
		if (typeof replacer === "function") it = replacer.call(iterable, String(i++), it);
		let key, value;
		if (Array.isArray(it)) if (it.length === 2) {
			key = it[0];
			value = it[1];
		} else throw new TypeError(`Expected [key, value] tuple: ${it}`);
		else if (it && it instanceof Object) {
			const keys = Object.keys(it);
			if (keys.length === 1) {
				key = keys[0];
				value = it[key];
			} else throw new TypeError(`Expected tuple with one key, not ${keys.length} keys`);
		} else key = it;
		pairs.items.push(createPair(key, value, ctx));
	}
	return pairs;
}
const pairs = {
	collection: "seq",
	default: false,
	tag: "tag:yaml.org,2002:pairs",
	resolve: resolvePairs,
	createNode: createPairs
};
//#endregion
//#region node_modules/yaml/browser/dist/schema/yaml-1.1/omap.js
var YAMLOMap = class YAMLOMap extends YAMLSeq {
	constructor() {
		super();
		this.add = YAMLMap.prototype.add.bind(this);
		this.delete = YAMLMap.prototype.delete.bind(this);
		this.get = YAMLMap.prototype.get.bind(this);
		this.has = YAMLMap.prototype.has.bind(this);
		this.set = YAMLMap.prototype.set.bind(this);
		this.tag = YAMLOMap.tag;
	}
	/**
	* If `ctx` is given, the return type is actually `Map<unknown, unknown>`,
	* but TypeScript won't allow widening the signature of a child method.
	*/
	toJSON(_, ctx) {
		if (!ctx) return super.toJSON(_);
		const map = /* @__PURE__ */ new Map();
		if (ctx?.onCreate) ctx.onCreate(map);
		for (const pair of this.items) {
			let key, value;
			if (isPair(pair)) {
				key = toJS(pair.key, "", ctx);
				value = toJS(pair.value, key, ctx);
			} else key = toJS(pair, "", ctx);
			if (map.has(key)) throw new Error("Ordered maps must not include duplicate keys");
			map.set(key, value);
		}
		return map;
	}
	static from(schema, iterable, ctx) {
		const pairs = createPairs(schema, iterable, ctx);
		const omap = new this();
		omap.items = pairs.items;
		return omap;
	}
};
YAMLOMap.tag = "tag:yaml.org,2002:omap";
const omap = {
	collection: "seq",
	identify: (value) => value instanceof Map,
	nodeClass: YAMLOMap,
	default: false,
	tag: "tag:yaml.org,2002:omap",
	resolve(seq, onError) {
		const pairs = resolvePairs(seq, onError);
		const seenKeys = [];
		for (const { key } of pairs.items) if (isScalar$1(key)) if (seenKeys.includes(key.value)) onError(`Ordered maps must not include duplicate keys: ${key.value}`);
		else seenKeys.push(key.value);
		return Object.assign(new YAMLOMap(), pairs);
	},
	createNode: (schema, iterable, ctx) => YAMLOMap.from(schema, iterable, ctx)
};
//#endregion
//#region node_modules/yaml/browser/dist/schema/yaml-1.1/bool.js
function boolStringify({ value, source }, ctx) {
	if (source && (value ? trueTag : falseTag).test.test(source)) return source;
	return value ? ctx.options.trueStr : ctx.options.falseStr;
}
const trueTag = {
	identify: (value) => value === true,
	default: true,
	tag: "tag:yaml.org,2002:bool",
	test: /^(?:Y|y|[Yy]es|YES|[Tt]rue|TRUE|[Oo]n|ON)$/,
	resolve: () => new Scalar(true),
	stringify: boolStringify
};
const falseTag = {
	identify: (value) => value === false,
	default: true,
	tag: "tag:yaml.org,2002:bool",
	test: /^(?:N|n|[Nn]o|NO|[Ff]alse|FALSE|[Oo]ff|OFF)$/,
	resolve: () => new Scalar(false),
	stringify: boolStringify
};
//#endregion
//#region node_modules/yaml/browser/dist/schema/yaml-1.1/float.js
const floatNaN = {
	identify: (value) => typeof value === "number",
	default: true,
	tag: "tag:yaml.org,2002:float",
	test: /^(?:[-+]?\.(?:inf|Inf|INF)|\.nan|\.NaN|\.NAN)$/,
	resolve: (str) => str.slice(-3).toLowerCase() === "nan" ? NaN : str[0] === "-" ? Number.NEGATIVE_INFINITY : Number.POSITIVE_INFINITY,
	stringify: stringifyNumber
};
const floatExp = {
	identify: (value) => typeof value === "number",
	default: true,
	tag: "tag:yaml.org,2002:float",
	format: "EXP",
	test: /^[-+]?(?:[0-9][0-9_]*)?(?:\.[0-9_]*)?[eE][-+]?[0-9]+$/,
	resolve: (str) => parseFloat(str.replace(/_/g, "")),
	stringify(node) {
		const num = Number(node.value);
		return isFinite(num) ? num.toExponential() : stringifyNumber(node);
	}
};
const float = {
	identify: (value) => typeof value === "number",
	default: true,
	tag: "tag:yaml.org,2002:float",
	test: /^[-+]?(?:[0-9][0-9_]*)?\.[0-9_]*$/,
	resolve(str) {
		const node = new Scalar(parseFloat(str.replace(/_/g, "")));
		const dot = str.indexOf(".");
		if (dot !== -1) {
			const f = str.substring(dot + 1).replace(/_/g, "");
			if (f[f.length - 1] === "0") node.minFractionDigits = f.length;
		}
		return node;
	},
	stringify: stringifyNumber
};
//#endregion
//#region node_modules/yaml/browser/dist/schema/yaml-1.1/int.js
const intIdentify = (value) => typeof value === "bigint" || Number.isInteger(value);
function intResolve(str, offset, radix, { intAsBigInt }) {
	const sign = str[0];
	if (sign === "-" || sign === "+") offset += 1;
	str = str.substring(offset).replace(/_/g, "");
	if (intAsBigInt) {
		switch (radix) {
			case 2:
				str = `0b${str}`;
				break;
			case 8:
				str = `0o${str}`;
				break;
			case 16:
				str = `0x${str}`;
				break;
		}
		const n = BigInt(str);
		return sign === "-" ? BigInt(-1) * n : n;
	}
	const n = parseInt(str, radix);
	return sign === "-" ? -1 * n : n;
}
function intStringify(node, radix, prefix) {
	const { value } = node;
	if (intIdentify(value)) {
		const str = value.toString(radix);
		return value < 0 ? "-" + prefix + str.substr(1) : prefix + str;
	}
	return stringifyNumber(node);
}
const intBin = {
	identify: intIdentify,
	default: true,
	tag: "tag:yaml.org,2002:int",
	format: "BIN",
	test: /^[-+]?0b[0-1_]+$/,
	resolve: (str, _onError, opt) => intResolve(str, 2, 2, opt),
	stringify: (node) => intStringify(node, 2, "0b")
};
const intOct = {
	identify: intIdentify,
	default: true,
	tag: "tag:yaml.org,2002:int",
	format: "OCT",
	test: /^[-+]?0[0-7_]+$/,
	resolve: (str, _onError, opt) => intResolve(str, 1, 8, opt),
	stringify: (node) => intStringify(node, 8, "0")
};
const int = {
	identify: intIdentify,
	default: true,
	tag: "tag:yaml.org,2002:int",
	test: /^[-+]?[0-9][0-9_]*$/,
	resolve: (str, _onError, opt) => intResolve(str, 0, 10, opt),
	stringify: stringifyNumber
};
const intHex = {
	identify: intIdentify,
	default: true,
	tag: "tag:yaml.org,2002:int",
	format: "HEX",
	test: /^[-+]?0x[0-9a-fA-F_]+$/,
	resolve: (str, _onError, opt) => intResolve(str, 2, 16, opt),
	stringify: (node) => intStringify(node, 16, "0x")
};
//#endregion
//#region node_modules/yaml/browser/dist/schema/yaml-1.1/set.js
var YAMLSet = class YAMLSet extends YAMLMap {
	constructor(schema) {
		super(schema);
		this.tag = YAMLSet.tag;
	}
	add(key) {
		let pair;
		if (isPair(key)) pair = key;
		else if (key && typeof key === "object" && "key" in key && "value" in key && key.value === null) pair = new Pair(key.key, null);
		else pair = new Pair(key, null);
		if (!findPair(this.items, pair.key)) this.items.push(pair);
	}
	/**
	* If `keepPair` is `true`, returns the Pair matching `key`.
	* Otherwise, returns the value of that Pair's key.
	*/
	get(key, keepPair) {
		const pair = findPair(this.items, key);
		return !keepPair && isPair(pair) ? isScalar$1(pair.key) ? pair.key.value : pair.key : pair;
	}
	set(key, value) {
		if (typeof value !== "boolean") throw new Error(`Expected boolean value for set(key, value) in a YAML set, not ${typeof value}`);
		const prev = findPair(this.items, key);
		if (prev && !value) this.items.splice(this.items.indexOf(prev), 1);
		else if (!prev && value) this.items.push(new Pair(key));
	}
	toJSON(_, ctx) {
		return super.toJSON(_, ctx, Set);
	}
	toString(ctx, onComment, onChompKeep) {
		if (!ctx) return JSON.stringify(this);
		if (this.hasAllNullValues(true)) return super.toString(Object.assign({}, ctx, { allNullValues: true }), onComment, onChompKeep);
		else throw new Error("Set items must all have null values");
	}
	static from(schema, iterable, ctx) {
		const { replacer } = ctx;
		const set = new this(schema);
		if (iterable && Symbol.iterator in Object(iterable)) for (let value of iterable) {
			if (typeof replacer === "function") value = replacer.call(iterable, value, value);
			set.items.push(createPair(value, null, ctx));
		}
		return set;
	}
};
YAMLSet.tag = "tag:yaml.org,2002:set";
const set = {
	collection: "map",
	identify: (value) => value instanceof Set,
	nodeClass: YAMLSet,
	default: false,
	tag: "tag:yaml.org,2002:set",
	createNode: (schema, iterable, ctx) => YAMLSet.from(schema, iterable, ctx),
	resolve(map, onError) {
		if (isMap(map)) if (map.hasAllNullValues(true)) return Object.assign(new YAMLSet(), map);
		else onError("Set items must all have null values");
		else onError("Expected a mapping for this tag");
		return map;
	}
};
//#endregion
//#region node_modules/yaml/browser/dist/schema/yaml-1.1/timestamp.js
/** Internal types handle bigint as number, because TS can't figure it out. */
function parseSexagesimal(str, asBigInt) {
	const sign = str[0];
	const parts = sign === "-" || sign === "+" ? str.substring(1) : str;
	const num = (n) => asBigInt ? BigInt(n) : Number(n);
	const res = parts.replace(/_/g, "").split(":").reduce((res, p) => res * num(60) + num(p), num(0));
	return sign === "-" ? num(-1) * res : res;
}
/**
* hhhh:mm:ss.sss
*
* Internal types handle bigint as number, because TS can't figure it out.
*/
function stringifySexagesimal(node) {
	let { value } = node;
	let num = (n) => n;
	if (typeof value === "bigint") num = (n) => BigInt(n);
	else if (isNaN(value) || !isFinite(value)) return stringifyNumber(node);
	let sign = "";
	if (value < 0) {
		sign = "-";
		value *= num(-1);
	}
	const _60 = num(60);
	const parts = [value % _60];
	if (value < 60) parts.unshift(0);
	else {
		value = (value - parts[0]) / _60;
		parts.unshift(value % _60);
		if (value >= 60) {
			value = (value - parts[0]) / _60;
			parts.unshift(value);
		}
	}
	return sign + parts.map((n) => String(n).padStart(2, "0")).join(":").replace(/000000\d*$/, "");
}
const intTime = {
	identify: (value) => typeof value === "bigint" || Number.isInteger(value),
	default: true,
	tag: "tag:yaml.org,2002:int",
	format: "TIME",
	test: /^[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+$/,
	resolve: (str, _onError, { intAsBigInt }) => parseSexagesimal(str, intAsBigInt),
	stringify: stringifySexagesimal
};
const floatTime = {
	identify: (value) => typeof value === "number",
	default: true,
	tag: "tag:yaml.org,2002:float",
	format: "TIME",
	test: /^[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\.[0-9_]*$/,
	resolve: (str) => parseSexagesimal(str, false),
	stringify: stringifySexagesimal
};
const timestamp = {
	identify: (value) => value instanceof Date,
	default: true,
	tag: "tag:yaml.org,2002:timestamp",
	test: RegExp("^([0-9]{4})-([0-9]{1,2})-([0-9]{1,2})(?:(?:t|T|[ \\t]+)([0-9]{1,2}):([0-9]{1,2}):([0-9]{1,2}(\\.[0-9]+)?)(?:[ \\t]*(Z|[-+][012]?[0-9](?::[0-9]{2})?))?)?$"),
	resolve(str) {
		const match = str.match(timestamp.test);
		if (!match) throw new Error("!!timestamp expects a date, starting with yyyy-mm-dd");
		const [, year, month, day, hour, minute, second] = match.map(Number);
		const millisec = match[7] ? Number((match[7] + "00").substr(1, 3)) : 0;
		let date = Date.UTC(year, month - 1, day, hour || 0, minute || 0, second || 0, millisec);
		const tz = match[8];
		if (tz && tz !== "Z") {
			let d = parseSexagesimal(tz, false);
			if (Math.abs(d) < 30) d *= 60;
			date -= 6e4 * d;
		}
		return new Date(date);
	},
	stringify: ({ value }) => value?.toISOString().replace(/(T00:00:00)?\.000Z$/, "") ?? ""
};
//#endregion
//#region node_modules/yaml/browser/dist/schema/yaml-1.1/schema.js
const schema = [
	map,
	seq,
	string,
	nullTag,
	trueTag,
	falseTag,
	intBin,
	intOct,
	int,
	intHex,
	floatNaN,
	floatExp,
	float,
	binary,
	merge,
	omap,
	pairs,
	set,
	intTime,
	floatTime,
	timestamp
];
//#endregion
//#region node_modules/yaml/browser/dist/schema/tags.js
const schemas = new Map([
	["core", schema$2],
	["failsafe", [
		map,
		seq,
		string
	]],
	["json", schema$1],
	["yaml11", schema],
	["yaml-1.1", schema]
]);
const tagsByName = {
	binary,
	bool: boolTag,
	float: float$1,
	floatExp: floatExp$1,
	floatNaN: floatNaN$1,
	floatTime,
	int: int$1,
	intHex: intHex$1,
	intOct: intOct$1,
	intTime,
	map,
	merge,
	null: nullTag,
	omap,
	pairs,
	seq,
	set,
	timestamp
};
const coreKnownTags = {
	"tag:yaml.org,2002:binary": binary,
	"tag:yaml.org,2002:merge": merge,
	"tag:yaml.org,2002:omap": omap,
	"tag:yaml.org,2002:pairs": pairs,
	"tag:yaml.org,2002:set": set,
	"tag:yaml.org,2002:timestamp": timestamp
};
function getTags(customTags, schemaName, addMergeTag) {
	const schemaTags = schemas.get(schemaName);
	if (schemaTags && !customTags) return addMergeTag && !schemaTags.includes(merge) ? schemaTags.concat(merge) : schemaTags.slice();
	let tags = schemaTags;
	if (!tags) if (Array.isArray(customTags)) tags = [];
	else {
		const keys = Array.from(schemas.keys()).filter((key) => key !== "yaml11").map((key) => JSON.stringify(key)).join(", ");
		throw new Error(`Unknown schema "${schemaName}"; use one of ${keys} or define customTags array`);
	}
	if (Array.isArray(customTags)) for (const tag of customTags) tags = tags.concat(tag);
	else if (typeof customTags === "function") tags = customTags(tags.slice());
	if (addMergeTag) tags = tags.concat(merge);
	return tags.reduce((tags, tag) => {
		const tagObj = typeof tag === "string" ? tagsByName[tag] : tag;
		if (!tagObj) {
			const tagName = JSON.stringify(tag);
			const keys = Object.keys(tagsByName).map((key) => JSON.stringify(key)).join(", ");
			throw new Error(`Unknown custom tag ${tagName}; use one of ${keys}`);
		}
		if (!tags.includes(tagObj)) tags.push(tagObj);
		return tags;
	}, []);
}
//#endregion
//#region node_modules/yaml/browser/dist/schema/Schema.js
const sortMapEntriesByKey = (a, b) => a.key < b.key ? -1 : a.key > b.key ? 1 : 0;
var Schema = class Schema {
	constructor({ compat, customTags, merge, resolveKnownTags, schema, sortMapEntries, toStringDefaults }) {
		this.compat = Array.isArray(compat) ? getTags(compat, "compat") : compat ? getTags(null, compat) : null;
		this.name = typeof schema === "string" && schema || "core";
		this.knownTags = resolveKnownTags ? coreKnownTags : {};
		this.tags = getTags(customTags, this.name, merge);
		this.toStringOptions = toStringDefaults ?? null;
		Object.defineProperty(this, MAP, { value: map });
		Object.defineProperty(this, SCALAR$1, { value: string });
		Object.defineProperty(this, SEQ, { value: seq });
		this.sortMapEntries = typeof sortMapEntries === "function" ? sortMapEntries : sortMapEntries === true ? sortMapEntriesByKey : null;
	}
	clone() {
		const copy = Object.create(Schema.prototype, Object.getOwnPropertyDescriptors(this));
		copy.tags = this.tags.slice();
		return copy;
	}
};
//#endregion
//#region node_modules/yaml/browser/dist/stringify/stringifyDocument.js
function stringifyDocument(doc, options) {
	const lines = [];
	let hasDirectives = options.directives === true;
	if (options.directives !== false && doc.directives) {
		const dir = doc.directives.toString(doc);
		if (dir) {
			lines.push(dir);
			hasDirectives = true;
		} else if (doc.directives.docStart) hasDirectives = true;
	}
	if (hasDirectives) lines.push("---");
	const ctx = createStringifyContext(doc, options);
	const { commentString } = ctx.options;
	if (doc.commentBefore) {
		if (lines.length !== 1) lines.unshift("");
		const cs = commentString(doc.commentBefore);
		lines.unshift(indentComment(cs, ""));
	}
	let chompKeep = false;
	let contentComment = null;
	if (doc.contents) {
		if (isNode(doc.contents)) {
			if (doc.contents.spaceBefore && hasDirectives) lines.push("");
			if (doc.contents.commentBefore) {
				const cs = commentString(doc.contents.commentBefore);
				lines.push(indentComment(cs, ""));
			}
			ctx.forceBlockIndent = !!doc.comment;
			contentComment = doc.contents.comment;
		}
		const onChompKeep = contentComment ? void 0 : () => chompKeep = true;
		let body = stringify$2(doc.contents, ctx, () => contentComment = null, onChompKeep);
		if (contentComment) body += lineComment(body, "", commentString(contentComment));
		if ((body[0] === "|" || body[0] === ">") && lines[lines.length - 1] === "---") lines[lines.length - 1] = `--- ${body}`;
		else lines.push(body);
	} else lines.push(stringify$2(doc.contents, ctx));
	if (doc.directives?.docEnd) if (doc.comment) {
		const cs = commentString(doc.comment);
		if (cs.includes("\n")) {
			lines.push("...");
			lines.push(indentComment(cs, ""));
		} else lines.push(`... ${cs}`);
	} else lines.push("...");
	else {
		let dc = doc.comment;
		if (dc && chompKeep) dc = dc.replace(/^\n+/, "");
		if (dc) {
			if ((!chompKeep || contentComment) && lines[lines.length - 1] !== "") lines.push("");
			lines.push(indentComment(commentString(dc), ""));
		}
	}
	return lines.join("\n") + "\n";
}
//#endregion
//#region node_modules/yaml/browser/dist/doc/Document.js
var Document = class Document {
	constructor(value, replacer, options) {
		/** A comment before this Document */
		this.commentBefore = null;
		/** A comment immediately after this Document */
		this.comment = null;
		/** Errors encountered during parsing. */
		this.errors = [];
		/** Warnings encountered during parsing. */
		this.warnings = [];
		Object.defineProperty(this, NODE_TYPE, { value: DOC });
		let _replacer = null;
		if (typeof replacer === "function" || Array.isArray(replacer)) _replacer = replacer;
		else if (options === void 0 && replacer) {
			options = replacer;
			replacer = void 0;
		}
		const opt = Object.assign({
			intAsBigInt: false,
			keepSourceTokens: false,
			logLevel: "warn",
			prettyErrors: true,
			strict: true,
			stringKeys: false,
			uniqueKeys: true,
			version: "1.2"
		}, options);
		this.options = opt;
		let { version } = opt;
		if (options?._directives) {
			this.directives = options._directives.atDocument();
			if (this.directives.yaml.explicit) version = this.directives.yaml.version;
		} else this.directives = new Directives({ version });
		this.setSchema(version, options);
		this.contents = value === void 0 ? null : this.createNode(value, _replacer, options);
	}
	/**
	* Create a deep copy of this Document and its contents.
	*
	* Custom Node values that inherit from `Object` still refer to their original instances.
	*/
	clone() {
		const copy = Object.create(Document.prototype, { [NODE_TYPE]: { value: DOC } });
		copy.commentBefore = this.commentBefore;
		copy.comment = this.comment;
		copy.errors = this.errors.slice();
		copy.warnings = this.warnings.slice();
		copy.options = Object.assign({}, this.options);
		if (this.directives) copy.directives = this.directives.clone();
		copy.schema = this.schema.clone();
		copy.contents = isNode(this.contents) ? this.contents.clone(copy.schema) : this.contents;
		if (this.range) copy.range = this.range.slice();
		return copy;
	}
	/** Adds a value to the document. */
	add(value) {
		if (assertCollection(this.contents)) this.contents.add(value);
	}
	/** Adds a value to the document. */
	addIn(path, value) {
		if (assertCollection(this.contents)) this.contents.addIn(path, value);
	}
	/**
	* Create a new `Alias` node, ensuring that the target `node` has the required anchor.
	*
	* If `node` already has an anchor, `name` is ignored.
	* Otherwise, the `node.anchor` value will be set to `name`,
	* or if an anchor with that name is already present in the document,
	* `name` will be used as a prefix for a new unique anchor.
	* If `name` is undefined, the generated anchor will use 'a' as a prefix.
	*/
	createAlias(node, name) {
		if (!node.anchor) {
			const prev = anchorNames(this);
			node.anchor = !name || prev.has(name) ? findNewAnchor(name || "a", prev) : name;
		}
		return new Alias(node.anchor);
	}
	createNode(value, replacer, options) {
		let _replacer = void 0;
		if (typeof replacer === "function") {
			value = replacer.call({ "": value }, "", value);
			_replacer = replacer;
		} else if (Array.isArray(replacer)) {
			const keyToStr = (v) => typeof v === "number" || v instanceof String || v instanceof Number;
			const asStr = replacer.filter(keyToStr).map(String);
			if (asStr.length > 0) replacer = replacer.concat(asStr);
			_replacer = replacer;
		} else if (options === void 0 && replacer) {
			options = replacer;
			replacer = void 0;
		}
		const { aliasDuplicateObjects, anchorPrefix, flow, keepUndefined, onTagObj, tag } = options ?? {};
		const { onAnchor, setAnchors, sourceObjects } = createNodeAnchors(this, anchorPrefix || "a");
		const ctx = {
			aliasDuplicateObjects: aliasDuplicateObjects ?? true,
			keepUndefined: keepUndefined ?? false,
			onAnchor,
			onTagObj,
			replacer: _replacer,
			schema: this.schema,
			sourceObjects
		};
		const node = createNode(value, tag, ctx);
		if (flow && isCollection$1(node)) node.flow = true;
		setAnchors();
		return node;
	}
	/**
	* Convert a key and a value into a `Pair` using the current schema,
	* recursively wrapping all values as `Scalar` or `Collection` nodes.
	*/
	createPair(key, value, options = {}) {
		return new Pair(this.createNode(key, null, options), this.createNode(value, null, options));
	}
	/**
	* Removes a value from the document.
	* @returns `true` if the item was found and removed.
	*/
	delete(key) {
		return assertCollection(this.contents) ? this.contents.delete(key) : false;
	}
	/**
	* Removes a value from the document.
	* @returns `true` if the item was found and removed.
	*/
	deleteIn(path) {
		if (isEmptyPath(path)) {
			if (this.contents == null) return false;
			this.contents = null;
			return true;
		}
		return assertCollection(this.contents) ? this.contents.deleteIn(path) : false;
	}
	/**
	* Returns item at `key`, or `undefined` if not found. By default unwraps
	* scalar values from their surrounding node; to disable set `keepScalar` to
	* `true` (collections are always returned intact).
	*/
	get(key, keepScalar) {
		return isCollection$1(this.contents) ? this.contents.get(key, keepScalar) : void 0;
	}
	/**
	* Returns item at `path`, or `undefined` if not found. By default unwraps
	* scalar values from their surrounding node; to disable set `keepScalar` to
	* `true` (collections are always returned intact).
	*/
	getIn(path, keepScalar) {
		if (isEmptyPath(path)) return !keepScalar && isScalar$1(this.contents) ? this.contents.value : this.contents;
		return isCollection$1(this.contents) ? this.contents.getIn(path, keepScalar) : void 0;
	}
	/**
	* Checks if the document includes a value with the key `key`.
	*/
	has(key) {
		return isCollection$1(this.contents) ? this.contents.has(key) : false;
	}
	/**
	* Checks if the document includes a value at `path`.
	*/
	hasIn(path) {
		if (isEmptyPath(path)) return this.contents !== void 0;
		return isCollection$1(this.contents) ? this.contents.hasIn(path) : false;
	}
	/**
	* Sets a value in this document. For `!!set`, `value` needs to be a
	* boolean to add/remove the item from the set.
	*/
	set(key, value) {
		if (this.contents == null) this.contents = collectionFromPath(this.schema, [key], value);
		else if (assertCollection(this.contents)) this.contents.set(key, value);
	}
	/**
	* Sets a value in this document. For `!!set`, `value` needs to be a
	* boolean to add/remove the item from the set.
	*/
	setIn(path, value) {
		if (isEmptyPath(path)) this.contents = value;
		else if (this.contents == null) this.contents = collectionFromPath(this.schema, Array.from(path), value);
		else if (assertCollection(this.contents)) this.contents.setIn(path, value);
	}
	/**
	* Change the YAML version and schema used by the document.
	* A `null` version disables support for directives, explicit tags, anchors, and aliases.
	* It also requires the `schema` option to be given as a `Schema` instance value.
	*
	* Overrides all previously set schema options.
	*/
	setSchema(version, options = {}) {
		if (typeof version === "number") version = String(version);
		let opt;
		switch (version) {
			case "1.1":
				if (this.directives) this.directives.yaml.version = "1.1";
				else this.directives = new Directives({ version: "1.1" });
				opt = {
					resolveKnownTags: false,
					schema: "yaml-1.1"
				};
				break;
			case "1.2":
			case "next":
				if (this.directives) this.directives.yaml.version = version;
				else this.directives = new Directives({ version });
				opt = {
					resolveKnownTags: true,
					schema: "core"
				};
				break;
			case null:
				if (this.directives) delete this.directives;
				opt = null;
				break;
			default: {
				const sv = JSON.stringify(version);
				throw new Error(`Expected '1.1', '1.2' or null as first argument, but found: ${sv}`);
			}
		}
		if (options.schema instanceof Object) this.schema = options.schema;
		else if (opt) this.schema = new Schema(Object.assign(opt, options));
		else throw new Error(`With a null YAML version, the { schema: Schema } option is required`);
	}
	toJS({ json, jsonArg, mapAsMap, maxAliasCount, onAnchor, reviver } = {}) {
		const ctx = {
			anchors: /* @__PURE__ */ new Map(),
			doc: this,
			keep: !json,
			mapAsMap: mapAsMap === true,
			mapKeyWarned: false,
			maxAliasCount: typeof maxAliasCount === "number" ? maxAliasCount : 100
		};
		const res = toJS(this.contents, jsonArg ?? "", ctx);
		if (typeof onAnchor === "function") for (const { count, res } of ctx.anchors.values()) onAnchor(res, count);
		return typeof reviver === "function" ? applyReviver(reviver, { "": res }, "", res) : res;
	}
	/**
	* A JSON representation of the document `contents`.
	*
	* @param jsonArg Used by `JSON.stringify` to indicate the array index or
	*   property name.
	*/
	toJSON(jsonArg, onAnchor) {
		return this.toJS({
			json: true,
			jsonArg,
			mapAsMap: false,
			onAnchor
		});
	}
	/** A YAML representation of the document. */
	toString(options = {}) {
		if (this.errors.length > 0) throw new Error("Document with errors cannot be stringified");
		if ("indent" in options && (!Number.isInteger(options.indent) || Number(options.indent) <= 0)) {
			const s = JSON.stringify(options.indent);
			throw new Error(`"indent" option must be a positive integer, not ${s}`);
		}
		return stringifyDocument(this, options);
	}
};
function assertCollection(contents) {
	if (isCollection$1(contents)) return true;
	throw new Error("Expected a YAML collection as document contents");
}
//#endregion
//#region node_modules/yaml/browser/dist/errors.js
var YAMLError = class extends Error {
	constructor(name, pos, code, message) {
		super();
		this.name = name;
		this.code = code;
		this.message = message;
		this.pos = pos;
	}
};
var YAMLParseError = class extends YAMLError {
	constructor(pos, code, message) {
		super("YAMLParseError", pos, code, message);
	}
};
var YAMLWarning = class extends YAMLError {
	constructor(pos, code, message) {
		super("YAMLWarning", pos, code, message);
	}
};
const prettifyError = (src, lc) => (error) => {
	if (error.pos[0] === -1) return;
	error.linePos = error.pos.map((pos) => lc.linePos(pos));
	const { line, col } = error.linePos[0];
	error.message += ` at line ${line}, column ${col}`;
	let ci = col - 1;
	let lineStr = src.substring(lc.lineStarts[line - 1], lc.lineStarts[line]).replace(/[\n\r]+$/, "");
	if (ci >= 60 && lineStr.length > 80) {
		const trimStart = Math.min(ci - 39, lineStr.length - 79);
		lineStr = "…" + lineStr.substring(trimStart);
		ci -= trimStart - 1;
	}
	if (lineStr.length > 80) lineStr = lineStr.substring(0, 79) + "…";
	if (line > 1 && /^ *$/.test(lineStr.substring(0, ci))) {
		let prev = src.substring(lc.lineStarts[line - 2], lc.lineStarts[line - 1]);
		if (prev.length > 80) prev = prev.substring(0, 79) + "…\n";
		lineStr = prev + lineStr;
	}
	if (/[^ ]/.test(lineStr)) {
		let count = 1;
		const end = error.linePos[1];
		if (end?.line === line && end.col > col) count = Math.max(1, Math.min(end.col - col, 80 - ci));
		const pointer = " ".repeat(ci) + "^".repeat(count);
		error.message += `:\n\n${lineStr}\n${pointer}\n`;
	}
};
//#endregion
//#region node_modules/yaml/browser/dist/compose/resolve-props.js
function resolveProps(tokens, { flow, indicator, next, offset, onError, parentIndent, startOnNewline }) {
	let spaceBefore = false;
	let atNewline = startOnNewline;
	let hasSpace = startOnNewline;
	let comment = "";
	let commentSep = "";
	let hasNewline = false;
	let reqSpace = false;
	let tab = null;
	let anchor = null;
	let tag = null;
	let newlineAfterProp = null;
	let comma = null;
	let found = null;
	let start = null;
	for (const token of tokens) {
		if (reqSpace) {
			if (token.type !== "space" && token.type !== "newline" && token.type !== "comma") onError(token.offset, "MISSING_CHAR", "Tags and anchors must be separated from the next token by white space");
			reqSpace = false;
		}
		if (tab) {
			if (atNewline && token.type !== "comment" && token.type !== "newline") onError(tab, "TAB_AS_INDENT", "Tabs are not allowed as indentation");
			tab = null;
		}
		switch (token.type) {
			case "space":
				if (!flow && (indicator !== "doc-start" || next?.type !== "flow-collection") && token.source.includes("	")) tab = token;
				hasSpace = true;
				break;
			case "comment": {
				if (!hasSpace) onError(token, "MISSING_CHAR", "Comments must be separated from other tokens by white space characters");
				const cb = token.source.substring(1) || " ";
				if (!comment) comment = cb;
				else comment += commentSep + cb;
				commentSep = "";
				atNewline = false;
				break;
			}
			case "newline":
				if (atNewline) {
					if (comment) comment += token.source;
					else if (!found || indicator !== "seq-item-ind") spaceBefore = true;
				} else commentSep += token.source;
				atNewline = true;
				hasNewline = true;
				if (anchor || tag) newlineAfterProp = token;
				hasSpace = true;
				break;
			case "anchor":
				if (anchor) onError(token, "MULTIPLE_ANCHORS", "A node can have at most one anchor");
				if (token.source.endsWith(":")) onError(token.offset + token.source.length - 1, "BAD_ALIAS", "Anchor ending in : is ambiguous", true);
				anchor = token;
				start ?? (start = token.offset);
				atNewline = false;
				hasSpace = false;
				reqSpace = true;
				break;
			case "tag":
				if (tag) onError(token, "MULTIPLE_TAGS", "A node can have at most one tag");
				tag = token;
				start ?? (start = token.offset);
				atNewline = false;
				hasSpace = false;
				reqSpace = true;
				break;
			case indicator:
				if (anchor || tag) onError(token, "BAD_PROP_ORDER", `Anchors and tags must be after the ${token.source} indicator`);
				if (found) onError(token, "UNEXPECTED_TOKEN", `Unexpected ${token.source} in ${flow ?? "collection"}`);
				found = token;
				atNewline = indicator === "seq-item-ind" || indicator === "explicit-key-ind";
				hasSpace = false;
				break;
			case "comma": if (flow) {
				if (comma) onError(token, "UNEXPECTED_TOKEN", `Unexpected , in ${flow}`);
				comma = token;
				atNewline = false;
				hasSpace = false;
				break;
			}
			default:
				onError(token, "UNEXPECTED_TOKEN", `Unexpected ${token.type} token`);
				atNewline = false;
				hasSpace = false;
		}
	}
	const last = tokens[tokens.length - 1];
	const end = last ? last.offset + last.source.length : offset;
	if (reqSpace && next && next.type !== "space" && next.type !== "newline" && next.type !== "comma" && (next.type !== "scalar" || next.source !== "")) onError(next.offset, "MISSING_CHAR", "Tags and anchors must be separated from the next token by white space");
	if (tab && (atNewline && tab.indent <= parentIndent || next?.type === "block-map" || next?.type === "block-seq")) onError(tab, "TAB_AS_INDENT", "Tabs are not allowed as indentation");
	return {
		comma,
		found,
		spaceBefore,
		comment,
		hasNewline,
		anchor,
		tag,
		newlineAfterProp,
		end,
		start: start ?? end
	};
}
//#endregion
//#region node_modules/yaml/browser/dist/compose/util-contains-newline.js
function containsNewline(key) {
	if (!key) return null;
	switch (key.type) {
		case "alias":
		case "scalar":
		case "double-quoted-scalar":
		case "single-quoted-scalar":
			if (key.source.includes("\n")) return true;
			if (key.end) {
				for (const st of key.end) if (st.type === "newline") return true;
			}
			return false;
		case "flow-collection":
			for (const it of key.items) {
				for (const st of it.start) if (st.type === "newline") return true;
				if (it.sep) {
					for (const st of it.sep) if (st.type === "newline") return true;
				}
				if (containsNewline(it.key) || containsNewline(it.value)) return true;
			}
			return false;
		default: return true;
	}
}
//#endregion
//#region node_modules/yaml/browser/dist/compose/util-flow-indent-check.js
function flowIndentCheck(indent, fc, onError) {
	if (fc?.type === "flow-collection") {
		const end = fc.end[0];
		if (end.indent === indent && (end.source === "]" || end.source === "}") && containsNewline(fc)) onError(end, "BAD_INDENT", "Flow end indicator should be more indented than parent", true);
	}
}
//#endregion
//#region node_modules/yaml/browser/dist/compose/util-map-includes.js
function mapIncludes(ctx, items, search) {
	const { uniqueKeys } = ctx.options;
	if (uniqueKeys === false) return false;
	const isEqual = typeof uniqueKeys === "function" ? uniqueKeys : (a, b) => a === b || isScalar$1(a) && isScalar$1(b) && a.value === b.value;
	return items.some((pair) => isEqual(pair.key, search));
}
//#endregion
//#region node_modules/yaml/browser/dist/compose/resolve-block-map.js
const startColMsg = "All mapping items must start at the same column";
function resolveBlockMap({ composeNode, composeEmptyNode }, ctx, bm, onError, tag) {
	const map = new (tag?.nodeClass ?? YAMLMap)(ctx.schema);
	if (ctx.atRoot) ctx.atRoot = false;
	let offset = bm.offset;
	let commentEnd = null;
	for (const collItem of bm.items) {
		const { start, key, sep, value } = collItem;
		const keyProps = resolveProps(start, {
			indicator: "explicit-key-ind",
			next: key ?? sep?.[0],
			offset,
			onError,
			parentIndent: bm.indent,
			startOnNewline: true
		});
		const implicitKey = !keyProps.found;
		if (implicitKey) {
			if (key) {
				if (key.type === "block-seq") onError(offset, "BLOCK_AS_IMPLICIT_KEY", "A block sequence may not be used as an implicit map key");
				else if ("indent" in key && key.indent !== bm.indent) onError(offset, "BAD_INDENT", startColMsg);
			}
			if (!keyProps.anchor && !keyProps.tag && !sep) {
				commentEnd = keyProps.end;
				if (keyProps.comment) if (map.comment) map.comment += "\n" + keyProps.comment;
				else map.comment = keyProps.comment;
				continue;
			}
			if (keyProps.newlineAfterProp || containsNewline(key)) onError(key ?? start[start.length - 1], "MULTILINE_IMPLICIT_KEY", "Implicit keys need to be on a single line");
		} else if (keyProps.found?.indent !== bm.indent) onError(offset, "BAD_INDENT", startColMsg);
		ctx.atKey = true;
		const keyStart = keyProps.end;
		const keyNode = key ? composeNode(ctx, key, keyProps, onError) : composeEmptyNode(ctx, keyStart, start, null, keyProps, onError);
		if (ctx.schema.compat) flowIndentCheck(bm.indent, key, onError);
		ctx.atKey = false;
		if (mapIncludes(ctx, map.items, keyNode)) onError(keyStart, "DUPLICATE_KEY", "Map keys must be unique");
		const valueProps = resolveProps(sep ?? [], {
			indicator: "map-value-ind",
			next: value,
			offset: keyNode.range[2],
			onError,
			parentIndent: bm.indent,
			startOnNewline: !key || key.type === "block-scalar"
		});
		offset = valueProps.end;
		if (valueProps.found) {
			if (implicitKey) {
				if (value?.type === "block-map" && !valueProps.hasNewline) onError(offset, "BLOCK_AS_IMPLICIT_KEY", "Nested mappings are not allowed in compact mappings");
				if (ctx.options.strict && keyProps.start < valueProps.found.offset - 1024) onError(keyNode.range, "KEY_OVER_1024_CHARS", "The : indicator must be at most 1024 chars after the start of an implicit block mapping key");
			}
			const valueNode = value ? composeNode(ctx, value, valueProps, onError) : composeEmptyNode(ctx, offset, sep, null, valueProps, onError);
			if (ctx.schema.compat) flowIndentCheck(bm.indent, value, onError);
			offset = valueNode.range[2];
			const pair = new Pair(keyNode, valueNode);
			if (ctx.options.keepSourceTokens) pair.srcToken = collItem;
			map.items.push(pair);
		} else {
			if (implicitKey) onError(keyNode.range, "MISSING_CHAR", "Implicit map keys need to be followed by map values");
			if (valueProps.comment) if (keyNode.comment) keyNode.comment += "\n" + valueProps.comment;
			else keyNode.comment = valueProps.comment;
			const pair = new Pair(keyNode);
			if (ctx.options.keepSourceTokens) pair.srcToken = collItem;
			map.items.push(pair);
		}
	}
	if (commentEnd && commentEnd < offset) onError(commentEnd, "IMPOSSIBLE", "Map comment with trailing content");
	map.range = [
		bm.offset,
		offset,
		commentEnd ?? offset
	];
	return map;
}
//#endregion
//#region node_modules/yaml/browser/dist/compose/resolve-block-seq.js
function resolveBlockSeq({ composeNode, composeEmptyNode }, ctx, bs, onError, tag) {
	const seq = new (tag?.nodeClass ?? YAMLSeq)(ctx.schema);
	if (ctx.atRoot) ctx.atRoot = false;
	if (ctx.atKey) ctx.atKey = false;
	let offset = bs.offset;
	let commentEnd = null;
	for (const { start, value } of bs.items) {
		const props = resolveProps(start, {
			indicator: "seq-item-ind",
			next: value,
			offset,
			onError,
			parentIndent: bs.indent,
			startOnNewline: true
		});
		if (!props.found) if (props.anchor || props.tag || value) if (value?.type === "block-seq") onError(props.end, "BAD_INDENT", "All sequence items must start at the same column");
		else onError(offset, "MISSING_CHAR", "Sequence item without - indicator");
		else {
			commentEnd = props.end;
			if (props.comment) seq.comment = props.comment;
			continue;
		}
		const node = value ? composeNode(ctx, value, props, onError) : composeEmptyNode(ctx, props.end, start, null, props, onError);
		if (ctx.schema.compat) flowIndentCheck(bs.indent, value, onError);
		offset = node.range[2];
		seq.items.push(node);
	}
	seq.range = [
		bs.offset,
		offset,
		commentEnd ?? offset
	];
	return seq;
}
//#endregion
//#region node_modules/yaml/browser/dist/compose/resolve-end.js
function resolveEnd(end, offset, reqSpace, onError) {
	let comment = "";
	if (end) {
		let hasSpace = false;
		let sep = "";
		for (const token of end) {
			const { source, type } = token;
			switch (type) {
				case "space":
					hasSpace = true;
					break;
				case "comment": {
					if (reqSpace && !hasSpace) onError(token, "MISSING_CHAR", "Comments must be separated from other tokens by white space characters");
					const cb = source.substring(1) || " ";
					if (!comment) comment = cb;
					else comment += sep + cb;
					sep = "";
					break;
				}
				case "newline":
					if (comment) sep += source;
					hasSpace = true;
					break;
				default: onError(token, "UNEXPECTED_TOKEN", `Unexpected ${type} at node end`);
			}
			offset += source.length;
		}
	}
	return {
		comment,
		offset
	};
}
//#endregion
//#region node_modules/yaml/browser/dist/compose/resolve-flow-collection.js
const blockMsg = "Block collections are not allowed within flow collections";
const isBlock = (token) => token && (token.type === "block-map" || token.type === "block-seq");
function resolveFlowCollection({ composeNode, composeEmptyNode }, ctx, fc, onError, tag) {
	const isMap = fc.start.source === "{";
	const fcName = isMap ? "flow map" : "flow sequence";
	const coll = new (tag?.nodeClass ?? (isMap ? YAMLMap : YAMLSeq))(ctx.schema);
	coll.flow = true;
	const atRoot = ctx.atRoot;
	if (atRoot) ctx.atRoot = false;
	if (ctx.atKey) ctx.atKey = false;
	let offset = fc.offset + fc.start.source.length;
	for (let i = 0; i < fc.items.length; ++i) {
		const collItem = fc.items[i];
		const { start, key, sep, value } = collItem;
		const props = resolveProps(start, {
			flow: fcName,
			indicator: "explicit-key-ind",
			next: key ?? sep?.[0],
			offset,
			onError,
			parentIndent: fc.indent,
			startOnNewline: false
		});
		if (!props.found) {
			if (!props.anchor && !props.tag && !sep && !value) {
				if (i === 0 && props.comma) onError(props.comma, "UNEXPECTED_TOKEN", `Unexpected , in ${fcName}`);
				else if (i < fc.items.length - 1) onError(props.start, "UNEXPECTED_TOKEN", `Unexpected empty item in ${fcName}`);
				if (props.comment) if (coll.comment) coll.comment += "\n" + props.comment;
				else coll.comment = props.comment;
				offset = props.end;
				continue;
			}
			if (!isMap && ctx.options.strict && containsNewline(key)) onError(key, "MULTILINE_IMPLICIT_KEY", "Implicit keys of flow sequence pairs need to be on a single line");
		}
		if (i === 0) {
			if (props.comma) onError(props.comma, "UNEXPECTED_TOKEN", `Unexpected , in ${fcName}`);
		} else {
			if (!props.comma) onError(props.start, "MISSING_CHAR", `Missing , between ${fcName} items`);
			if (props.comment) {
				let prevItemComment = "";
				loop: for (const st of start) switch (st.type) {
					case "comma":
					case "space": break;
					case "comment":
						prevItemComment = st.source.substring(1);
						break loop;
					default: break loop;
				}
				if (prevItemComment) {
					let prev = coll.items[coll.items.length - 1];
					if (isPair(prev)) prev = prev.value ?? prev.key;
					if (prev.comment) prev.comment += "\n" + prevItemComment;
					else prev.comment = prevItemComment;
					props.comment = props.comment.substring(prevItemComment.length + 1);
				}
			}
		}
		if (!isMap && !sep && !props.found) {
			const valueNode = value ? composeNode(ctx, value, props, onError) : composeEmptyNode(ctx, props.end, sep, null, props, onError);
			coll.items.push(valueNode);
			offset = valueNode.range[2];
			if (isBlock(value)) onError(valueNode.range, "BLOCK_IN_FLOW", blockMsg);
		} else {
			ctx.atKey = true;
			const keyStart = props.end;
			const keyNode = key ? composeNode(ctx, key, props, onError) : composeEmptyNode(ctx, keyStart, start, null, props, onError);
			if (isBlock(key)) onError(keyNode.range, "BLOCK_IN_FLOW", blockMsg);
			ctx.atKey = false;
			const valueProps = resolveProps(sep ?? [], {
				flow: fcName,
				indicator: "map-value-ind",
				next: value,
				offset: keyNode.range[2],
				onError,
				parentIndent: fc.indent,
				startOnNewline: false
			});
			if (valueProps.found) {
				if (!isMap && !props.found && ctx.options.strict) {
					if (sep) for (const st of sep) {
						if (st === valueProps.found) break;
						if (st.type === "newline") {
							onError(st, "MULTILINE_IMPLICIT_KEY", "Implicit keys of flow sequence pairs need to be on a single line");
							break;
						}
					}
					if (props.start < valueProps.found.offset - 1024) onError(valueProps.found, "KEY_OVER_1024_CHARS", "The : indicator must be at most 1024 chars after the start of an implicit flow sequence key");
				}
			} else if (value) if ("source" in value && value.source?.[0] === ":") onError(value, "MISSING_CHAR", `Missing space after : in ${fcName}`);
			else onError(valueProps.start, "MISSING_CHAR", `Missing , or : between ${fcName} items`);
			const valueNode = value ? composeNode(ctx, value, valueProps, onError) : valueProps.found ? composeEmptyNode(ctx, valueProps.end, sep, null, valueProps, onError) : null;
			if (valueNode) {
				if (isBlock(value)) onError(valueNode.range, "BLOCK_IN_FLOW", blockMsg);
			} else if (valueProps.comment) if (keyNode.comment) keyNode.comment += "\n" + valueProps.comment;
			else keyNode.comment = valueProps.comment;
			const pair = new Pair(keyNode, valueNode);
			if (ctx.options.keepSourceTokens) pair.srcToken = collItem;
			if (isMap) {
				const map = coll;
				if (mapIncludes(ctx, map.items, keyNode)) onError(keyStart, "DUPLICATE_KEY", "Map keys must be unique");
				map.items.push(pair);
			} else {
				const map = new YAMLMap(ctx.schema);
				map.flow = true;
				map.items.push(pair);
				const endRange = (valueNode ?? keyNode).range;
				map.range = [
					keyNode.range[0],
					endRange[1],
					endRange[2]
				];
				coll.items.push(map);
			}
			offset = valueNode ? valueNode.range[2] : valueProps.end;
		}
	}
	const expectedEnd = isMap ? "}" : "]";
	const [ce, ...ee] = fc.end;
	let cePos = offset;
	if (ce?.source === expectedEnd) cePos = ce.offset + ce.source.length;
	else {
		const name = fcName[0].toUpperCase() + fcName.substring(1);
		const msg = atRoot ? `${name} must end with a ${expectedEnd}` : `${name} in block collection must be sufficiently indented and end with a ${expectedEnd}`;
		onError(offset, atRoot ? "MISSING_CHAR" : "BAD_INDENT", msg);
		if (ce && ce.source.length !== 1) ee.unshift(ce);
	}
	if (ee.length > 0) {
		const end = resolveEnd(ee, cePos, ctx.options.strict, onError);
		if (end.comment) if (coll.comment) coll.comment += "\n" + end.comment;
		else coll.comment = end.comment;
		coll.range = [
			fc.offset,
			cePos,
			end.offset
		];
	} else coll.range = [
		fc.offset,
		cePos,
		cePos
	];
	return coll;
}
//#endregion
//#region node_modules/yaml/browser/dist/compose/compose-collection.js
function resolveCollection(CN, ctx, token, onError, tagName, tag) {
	const coll = token.type === "block-map" ? resolveBlockMap(CN, ctx, token, onError, tag) : token.type === "block-seq" ? resolveBlockSeq(CN, ctx, token, onError, tag) : resolveFlowCollection(CN, ctx, token, onError, tag);
	const Coll = coll.constructor;
	if (tagName === "!" || tagName === Coll.tagName) {
		coll.tag = Coll.tagName;
		return coll;
	}
	if (tagName) coll.tag = tagName;
	return coll;
}
function composeCollection(CN, ctx, token, props, onError) {
	const tagToken = props.tag;
	const tagName = !tagToken ? null : ctx.directives.tagName(tagToken.source, (msg) => onError(tagToken, "TAG_RESOLVE_FAILED", msg));
	if (token.type === "block-seq") {
		const { anchor, newlineAfterProp: nl } = props;
		const lastProp = anchor && tagToken ? anchor.offset > tagToken.offset ? anchor : tagToken : anchor ?? tagToken;
		if (lastProp && (!nl || nl.offset < lastProp.offset)) onError(lastProp, "MISSING_CHAR", "Missing newline after block sequence props");
	}
	const expType = token.type === "block-map" ? "map" : token.type === "block-seq" ? "seq" : token.start.source === "{" ? "map" : "seq";
	if (!tagToken || !tagName || tagName === "!" || tagName === YAMLMap.tagName && expType === "map" || tagName === YAMLSeq.tagName && expType === "seq") return resolveCollection(CN, ctx, token, onError, tagName);
	let tag = ctx.schema.tags.find((t) => t.tag === tagName && t.collection === expType);
	if (!tag) {
		const kt = ctx.schema.knownTags[tagName];
		if (kt?.collection === expType) {
			ctx.schema.tags.push(Object.assign({}, kt, { default: false }));
			tag = kt;
		} else {
			if (kt) onError(tagToken, "BAD_COLLECTION_TYPE", `${kt.tag} used for ${expType} collection, but expects ${kt.collection ?? "scalar"}`, true);
			else onError(tagToken, "TAG_RESOLVE_FAILED", `Unresolved tag: ${tagName}`, true);
			return resolveCollection(CN, ctx, token, onError, tagName);
		}
	}
	const coll = resolveCollection(CN, ctx, token, onError, tagName, tag);
	const res = tag.resolve?.(coll, (msg) => onError(tagToken, "TAG_RESOLVE_FAILED", msg), ctx.options) ?? coll;
	const node = isNode(res) ? res : new Scalar(res);
	node.range = coll.range;
	node.tag = tagName;
	if (tag?.format) node.format = tag.format;
	return node;
}
//#endregion
//#region node_modules/yaml/browser/dist/compose/resolve-block-scalar.js
function resolveBlockScalar(ctx, scalar, onError) {
	const start = scalar.offset;
	const header = parseBlockScalarHeader(scalar, ctx.options.strict, onError);
	if (!header) return {
		value: "",
		type: null,
		comment: "",
		range: [
			start,
			start,
			start
		]
	};
	const type = header.mode === ">" ? Scalar.BLOCK_FOLDED : Scalar.BLOCK_LITERAL;
	const lines = scalar.source ? splitLines(scalar.source) : [];
	let chompStart = lines.length;
	for (let i = lines.length - 1; i >= 0; --i) {
		const content = lines[i][1];
		if (content === "" || content === "\r") chompStart = i;
		else break;
	}
	if (chompStart === 0) {
		const value = header.chomp === "+" && lines.length > 0 ? "\n".repeat(Math.max(1, lines.length - 1)) : "";
		let end = start + header.length;
		if (scalar.source) end += scalar.source.length;
		return {
			value,
			type,
			comment: header.comment,
			range: [
				start,
				end,
				end
			]
		};
	}
	let trimIndent = scalar.indent + header.indent;
	let offset = scalar.offset + header.length;
	let contentStart = 0;
	for (let i = 0; i < chompStart; ++i) {
		const [indent, content] = lines[i];
		if (content === "" || content === "\r") {
			if (header.indent === 0 && indent.length > trimIndent) trimIndent = indent.length;
		} else {
			if (indent.length < trimIndent) onError(offset + indent.length, "MISSING_CHAR", "Block scalars with more-indented leading empty lines must use an explicit indentation indicator");
			if (header.indent === 0) trimIndent = indent.length;
			contentStart = i;
			if (trimIndent === 0 && !ctx.atRoot) onError(offset, "BAD_INDENT", "Block scalar values in collections must be indented");
			break;
		}
		offset += indent.length + content.length + 1;
	}
	for (let i = lines.length - 1; i >= chompStart; --i) if (lines[i][0].length > trimIndent) chompStart = i + 1;
	let value = "";
	let sep = "";
	let prevMoreIndented = false;
	for (let i = 0; i < contentStart; ++i) value += lines[i][0].slice(trimIndent) + "\n";
	for (let i = contentStart; i < chompStart; ++i) {
		let [indent, content] = lines[i];
		offset += indent.length + content.length + 1;
		const crlf = content[content.length - 1] === "\r";
		if (crlf) content = content.slice(0, -1);
		/* istanbul ignore if already caught in lexer */
		if (content && indent.length < trimIndent) {
			const message = `Block scalar lines must not be less indented than their ${header.indent ? "explicit indentation indicator" : "first line"}`;
			onError(offset - content.length - (crlf ? 2 : 1), "BAD_INDENT", message);
			indent = "";
		}
		if (type === Scalar.BLOCK_LITERAL) {
			value += sep + indent.slice(trimIndent) + content;
			sep = "\n";
		} else if (indent.length > trimIndent || content[0] === "	") {
			if (sep === " ") sep = "\n";
			else if (!prevMoreIndented && sep === "\n") sep = "\n\n";
			value += sep + indent.slice(trimIndent) + content;
			sep = "\n";
			prevMoreIndented = true;
		} else if (content === "") if (sep === "\n") value += "\n";
		else sep = "\n";
		else {
			value += sep + content;
			sep = " ";
			prevMoreIndented = false;
		}
	}
	switch (header.chomp) {
		case "-": break;
		case "+":
			for (let i = chompStart; i < lines.length; ++i) value += "\n" + lines[i][0].slice(trimIndent);
			if (value[value.length - 1] !== "\n") value += "\n";
			break;
		default: value += "\n";
	}
	const end = start + header.length + scalar.source.length;
	return {
		value,
		type,
		comment: header.comment,
		range: [
			start,
			end,
			end
		]
	};
}
function parseBlockScalarHeader({ offset, props }, strict, onError) {
	/* istanbul ignore if should not happen */
	if (props[0].type !== "block-scalar-header") {
		onError(props[0], "IMPOSSIBLE", "Block scalar header not found");
		return null;
	}
	const { source } = props[0];
	const mode = source[0];
	let indent = 0;
	let chomp = "";
	let error = -1;
	for (let i = 1; i < source.length; ++i) {
		const ch = source[i];
		if (!chomp && (ch === "-" || ch === "+")) chomp = ch;
		else {
			const n = Number(ch);
			if (!indent && n) indent = n;
			else if (error === -1) error = offset + i;
		}
	}
	if (error !== -1) onError(error, "UNEXPECTED_TOKEN", `Block scalar header includes extra characters: ${source}`);
	let hasSpace = false;
	let comment = "";
	let length = source.length;
	for (let i = 1; i < props.length; ++i) {
		const token = props[i];
		switch (token.type) {
			case "space": hasSpace = true;
			case "newline":
				length += token.source.length;
				break;
			case "comment":
				if (strict && !hasSpace) onError(token, "MISSING_CHAR", "Comments must be separated from other tokens by white space characters");
				length += token.source.length;
				comment = token.source.substring(1);
				break;
			case "error":
				onError(token, "UNEXPECTED_TOKEN", token.message);
				length += token.source.length;
				break;
			/* istanbul ignore next should not happen */
			default: {
				onError(token, "UNEXPECTED_TOKEN", `Unexpected token in block scalar header: ${token.type}`);
				const ts = token.source;
				if (ts && typeof ts === "string") length += ts.length;
			}
		}
	}
	return {
		mode,
		indent,
		chomp,
		comment,
		length
	};
}
/** @returns Array of lines split up as `[indent, content]` */
function splitLines(source) {
	const split = source.split(/\n( *)/);
	const first = split[0];
	const m = first.match(/^( *)/);
	const lines = [m?.[1] ? [m[1], first.slice(m[1].length)] : ["", first]];
	for (let i = 1; i < split.length; i += 2) lines.push([split[i], split[i + 1]]);
	return lines;
}
//#endregion
//#region node_modules/yaml/browser/dist/compose/resolve-flow-scalar.js
function resolveFlowScalar(scalar, strict, onError) {
	const { offset, type, source, end } = scalar;
	let _type;
	let value;
	const _onError = (rel, code, msg) => onError(offset + rel, code, msg);
	switch (type) {
		case "scalar":
			_type = Scalar.PLAIN;
			value = plainValue(source, _onError);
			break;
		case "single-quoted-scalar":
			_type = Scalar.QUOTE_SINGLE;
			value = singleQuotedValue(source, _onError);
			break;
		case "double-quoted-scalar":
			_type = Scalar.QUOTE_DOUBLE;
			value = doubleQuotedValue(source, _onError);
			break;
		/* istanbul ignore next should not happen */
		default:
			onError(scalar, "UNEXPECTED_TOKEN", `Expected a flow scalar value, but found: ${type}`);
			return {
				value: "",
				type: null,
				comment: "",
				range: [
					offset,
					offset + source.length,
					offset + source.length
				]
			};
	}
	const valueEnd = offset + source.length;
	const re = resolveEnd(end, valueEnd, strict, onError);
	return {
		value,
		type: _type,
		comment: re.comment,
		range: [
			offset,
			valueEnd,
			re.offset
		]
	};
}
function plainValue(source, onError) {
	let badChar = "";
	switch (source[0]) {
		/* istanbul ignore next should not happen */
		case "	":
			badChar = "a tab character";
			break;
		case ",":
			badChar = "flow indicator character ,";
			break;
		case "%":
			badChar = "directive indicator character %";
			break;
		case "|":
		case ">":
			badChar = `block scalar indicator ${source[0]}`;
			break;
		case "@":
		case "`":
			badChar = `reserved character ${source[0]}`;
			break;
	}
	if (badChar) onError(0, "BAD_SCALAR_START", `Plain value cannot start with ${badChar}`);
	return foldLines(source);
}
function singleQuotedValue(source, onError) {
	if (source[source.length - 1] !== "'" || source.length === 1) onError(source.length, "MISSING_CHAR", "Missing closing 'quote");
	return foldLines(source.slice(1, -1)).replace(/''/g, "'");
}
function foldLines(source) {
	/**
	* The negative lookbehind here and in the `re` RegExp is to
	* prevent causing a polynomial search time in certain cases.
	*
	* The try-catch is for Safari, which doesn't support this yet:
	* https://caniuse.com/js-regexp-lookbehind
	*/
	let first, line;
	try {
		first = /* @__PURE__ */ new RegExp("(.*?)(?<![ 	])[ 	]*\r?\n", "sy");
		line = /* @__PURE__ */ new RegExp("[ 	]*(.*?)(?:(?<![ 	])[ 	]*)?\r?\n", "sy");
	} catch {
		first = /(.*?)[ \t]*\r?\n/sy;
		line = /[ \t]*(.*?)[ \t]*\r?\n/sy;
	}
	let match = first.exec(source);
	if (!match) return source;
	let res = match[1];
	let sep = " ";
	let pos = first.lastIndex;
	line.lastIndex = pos;
	while (match = line.exec(source)) {
		if (match[1] === "") if (sep === "\n") res += sep;
		else sep = "\n";
		else {
			res += sep + match[1];
			sep = " ";
		}
		pos = line.lastIndex;
	}
	const last = /[ \t]*(.*)/sy;
	last.lastIndex = pos;
	match = last.exec(source);
	return res + sep + (match?.[1] ?? "");
}
function doubleQuotedValue(source, onError) {
	let res = "";
	for (let i = 1; i < source.length - 1; ++i) {
		const ch = source[i];
		if (ch === "\r" && source[i + 1] === "\n") continue;
		if (ch === "\n") {
			const { fold, offset } = foldNewline(source, i);
			res += fold;
			i = offset;
		} else if (ch === "\\") {
			let next = source[++i];
			const cc = escapeCodes[next];
			if (cc) res += cc;
			else if (next === "\n") {
				next = source[i + 1];
				while (next === " " || next === "	") next = source[++i + 1];
			} else if (next === "\r" && source[i + 1] === "\n") {
				next = source[++i + 1];
				while (next === " " || next === "	") next = source[++i + 1];
			} else if (next === "x" || next === "u" || next === "U") {
				const length = next === "x" ? 2 : next === "u" ? 4 : 8;
				res += parseCharCode(source, i + 1, length, onError);
				i += length;
			} else {
				const raw = source.substr(i - 1, 2);
				onError(i - 1, "BAD_DQ_ESCAPE", `Invalid escape sequence ${raw}`);
				res += raw;
			}
		} else if (ch === " " || ch === "	") {
			const wsStart = i;
			let next = source[i + 1];
			while (next === " " || next === "	") next = source[++i + 1];
			if (next !== "\n" && !(next === "\r" && source[i + 2] === "\n")) res += i > wsStart ? source.slice(wsStart, i + 1) : ch;
		} else res += ch;
	}
	if (source[source.length - 1] !== "\"" || source.length === 1) onError(source.length, "MISSING_CHAR", "Missing closing \"quote");
	return res;
}
/**
* Fold a single newline into a space, multiple newlines to N - 1 newlines.
* Presumes `source[offset] === '\n'`
*/
function foldNewline(source, offset) {
	let fold = "";
	let ch = source[offset + 1];
	while (ch === " " || ch === "	" || ch === "\n" || ch === "\r") {
		if (ch === "\r" && source[offset + 2] !== "\n") break;
		if (ch === "\n") fold += "\n";
		offset += 1;
		ch = source[offset + 1];
	}
	if (!fold) fold = " ";
	return {
		fold,
		offset
	};
}
const escapeCodes = {
	"0": "\0",
	a: "\x07",
	b: "\b",
	e: "\x1B",
	f: "\f",
	n: "\n",
	r: "\r",
	t: "	",
	v: "\v",
	N: "",
	_: "\xA0",
	L: "\u2028",
	P: "\u2029",
	" ": " ",
	"\"": "\"",
	"/": "/",
	"\\": "\\",
	"	": "	"
};
function parseCharCode(source, offset, length, onError) {
	const cc = source.substr(offset, length);
	const code = cc.length === length && /^[0-9a-fA-F]+$/.test(cc) ? parseInt(cc, 16) : NaN;
	try {
		return String.fromCodePoint(code);
	} catch {
		const raw = source.substr(offset - 2, length + 2);
		onError(offset - 2, "BAD_DQ_ESCAPE", `Invalid escape sequence ${raw}`);
		return raw;
	}
}
//#endregion
//#region node_modules/yaml/browser/dist/compose/compose-scalar.js
function composeScalar(ctx, token, tagToken, onError) {
	const { value, type, comment, range } = token.type === "block-scalar" ? resolveBlockScalar(ctx, token, onError) : resolveFlowScalar(token, ctx.options.strict, onError);
	const tagName = tagToken ? ctx.directives.tagName(tagToken.source, (msg) => onError(tagToken, "TAG_RESOLVE_FAILED", msg)) : null;
	let tag;
	if (ctx.options.stringKeys && ctx.atKey) tag = ctx.schema[SCALAR$1];
	else if (tagName) tag = findScalarTagByName(ctx.schema, value, tagName, tagToken, onError);
	else if (token.type === "scalar") tag = findScalarTagByTest(ctx, value, token, onError);
	else tag = ctx.schema[SCALAR$1];
	let scalar;
	try {
		const res = tag.resolve(value, (msg) => onError(tagToken ?? token, "TAG_RESOLVE_FAILED", msg), ctx.options);
		scalar = isScalar$1(res) ? res : new Scalar(res);
	} catch (error) {
		const msg = error instanceof Error ? error.message : String(error);
		onError(tagToken ?? token, "TAG_RESOLVE_FAILED", msg);
		scalar = new Scalar(value);
	}
	scalar.range = range;
	scalar.source = value;
	if (type) scalar.type = type;
	if (tagName) scalar.tag = tagName;
	if (tag.format) scalar.format = tag.format;
	if (comment) scalar.comment = comment;
	return scalar;
}
function findScalarTagByName(schema, value, tagName, tagToken, onError) {
	if (tagName === "!") return schema[SCALAR$1];
	const matchWithTest = [];
	for (const tag of schema.tags) if (!tag.collection && tag.tag === tagName) if (tag.default && tag.test) matchWithTest.push(tag);
	else return tag;
	for (const tag of matchWithTest) if (tag.test?.test(value)) return tag;
	const kt = schema.knownTags[tagName];
	if (kt && !kt.collection) {
		schema.tags.push(Object.assign({}, kt, {
			default: false,
			test: void 0
		}));
		return kt;
	}
	onError(tagToken, "TAG_RESOLVE_FAILED", `Unresolved tag: ${tagName}`, tagName !== "tag:yaml.org,2002:str");
	return schema[SCALAR$1];
}
function findScalarTagByTest({ atKey, directives, schema }, value, token, onError) {
	const tag = schema.tags.find((tag) => (tag.default === true || atKey && tag.default === "key") && tag.test?.test(value)) || schema[SCALAR$1];
	if (schema.compat) {
		const compat = schema.compat.find((tag) => tag.default && tag.test?.test(value)) ?? schema[SCALAR$1];
		if (tag.tag !== compat.tag) onError(token, "TAG_RESOLVE_FAILED", `Value may be parsed as either ${directives.tagString(tag.tag)} or ${directives.tagString(compat.tag)}`, true);
	}
	return tag;
}
//#endregion
//#region node_modules/yaml/browser/dist/compose/util-empty-scalar-position.js
function emptyScalarPosition(offset, before, pos) {
	if (before) {
		pos ?? (pos = before.length);
		for (let i = pos - 1; i >= 0; --i) {
			let st = before[i];
			switch (st.type) {
				case "space":
				case "comment":
				case "newline":
					offset -= st.source.length;
					continue;
			}
			st = before[++i];
			while (st?.type === "space") {
				offset += st.source.length;
				st = before[++i];
			}
			break;
		}
	}
	return offset;
}
//#endregion
//#region node_modules/yaml/browser/dist/compose/compose-node.js
const CN = {
	composeNode,
	composeEmptyNode
};
function composeNode(ctx, token, props, onError) {
	const atKey = ctx.atKey;
	const { spaceBefore, comment, anchor, tag } = props;
	let node;
	let isSrcToken = true;
	switch (token.type) {
		case "alias":
			node = composeAlias(ctx, token, onError);
			if (anchor || tag) onError(token, "ALIAS_PROPS", "An alias node must not specify any properties");
			break;
		case "scalar":
		case "single-quoted-scalar":
		case "double-quoted-scalar":
		case "block-scalar":
			node = composeScalar(ctx, token, tag, onError);
			if (anchor) node.anchor = anchor.source.substring(1);
			break;
		case "block-map":
		case "block-seq":
		case "flow-collection":
			try {
				node = composeCollection(CN, ctx, token, props, onError);
				if (anchor) node.anchor = anchor.source.substring(1);
			} catch (error) {
				onError(token, "RESOURCE_EXHAUSTION", error instanceof Error ? error.message : String(error));
			}
			break;
		default:
			onError(token, "UNEXPECTED_TOKEN", token.type === "error" ? token.message : `Unsupported token (type: ${token.type})`);
			isSrcToken = false;
	}
	node ?? (node = composeEmptyNode(ctx, token.offset, void 0, null, props, onError));
	if (anchor && node.anchor === "") onError(anchor, "BAD_ALIAS", "Anchor cannot be an empty string");
	if (atKey && ctx.options.stringKeys && (!isScalar$1(node) || typeof node.value !== "string" || node.tag && node.tag !== "tag:yaml.org,2002:str")) onError(tag ?? token, "NON_STRING_KEY", "With stringKeys, all keys must be strings");
	if (spaceBefore) node.spaceBefore = true;
	if (comment) if (token.type === "scalar" && token.source === "") node.comment = comment;
	else node.commentBefore = comment;
	if (ctx.options.keepSourceTokens && isSrcToken) node.srcToken = token;
	return node;
}
function composeEmptyNode(ctx, offset, before, pos, { spaceBefore, comment, anchor, tag, end }, onError) {
	const node = composeScalar(ctx, {
		type: "scalar",
		offset: emptyScalarPosition(offset, before, pos),
		indent: -1,
		source: ""
	}, tag, onError);
	if (anchor) {
		node.anchor = anchor.source.substring(1);
		if (node.anchor === "") onError(anchor, "BAD_ALIAS", "Anchor cannot be an empty string");
	}
	if (spaceBefore) node.spaceBefore = true;
	if (comment) {
		node.comment = comment;
		node.range[2] = end;
	}
	return node;
}
function composeAlias({ options }, { offset, source, end }, onError) {
	const alias = new Alias(source.substring(1));
	if (alias.source === "") onError(offset, "BAD_ALIAS", "Alias cannot be an empty string");
	if (alias.source.endsWith(":")) onError(offset + source.length - 1, "BAD_ALIAS", "Alias ending in : is ambiguous", true);
	const valueEnd = offset + source.length;
	const re = resolveEnd(end, valueEnd, options.strict, onError);
	alias.range = [
		offset,
		valueEnd,
		re.offset
	];
	if (re.comment) alias.comment = re.comment;
	return alias;
}
//#endregion
//#region node_modules/yaml/browser/dist/compose/compose-doc.js
function composeDoc(options, directives, { offset, start, value, end }, onError) {
	const doc = new Document(void 0, Object.assign({ _directives: directives }, options));
	const ctx = {
		atKey: false,
		atRoot: true,
		directives: doc.directives,
		options: doc.options,
		schema: doc.schema
	};
	const props = resolveProps(start, {
		indicator: "doc-start",
		next: value ?? end?.[0],
		offset,
		onError,
		parentIndent: 0,
		startOnNewline: true
	});
	if (props.found) {
		doc.directives.docStart = true;
		if (value && (value.type === "block-map" || value.type === "block-seq") && !props.hasNewline) onError(props.end, "MISSING_CHAR", "Block collection cannot start on same line with directives-end marker");
	}
	doc.contents = value ? composeNode(ctx, value, props, onError) : composeEmptyNode(ctx, props.end, start, null, props, onError);
	const contentEnd = doc.contents.range[2];
	const re = resolveEnd(end, contentEnd, false, onError);
	if (re.comment) doc.comment = re.comment;
	doc.range = [
		offset,
		contentEnd,
		re.offset
	];
	return doc;
}
//#endregion
//#region node_modules/yaml/browser/dist/compose/composer.js
function getErrorPos(src) {
	if (typeof src === "number") return [src, src + 1];
	if (Array.isArray(src)) return src.length === 2 ? src : [src[0], src[1]];
	const { offset, source } = src;
	return [offset, offset + (typeof source === "string" ? source.length : 1)];
}
function parsePrelude(prelude) {
	let comment = "";
	let atComment = false;
	let afterEmptyLine = false;
	for (let i = 0; i < prelude.length; ++i) {
		const source = prelude[i];
		switch (source[0]) {
			case "#":
				comment += (comment === "" ? "" : afterEmptyLine ? "\n\n" : "\n") + (source.substring(1) || " ");
				atComment = true;
				afterEmptyLine = false;
				break;
			case "%":
				if (prelude[i + 1]?.[0] !== "#") i += 1;
				atComment = false;
				break;
			default:
				if (!atComment) afterEmptyLine = true;
				atComment = false;
		}
	}
	return {
		comment,
		afterEmptyLine
	};
}
/**
* Compose a stream of CST nodes into a stream of YAML Documents.
*
* ```ts
* import { Composer, Parser } from 'yaml'
*
* const src: string = ...
* const tokens = new Parser().parse(src)
* const docs = new Composer().compose(tokens)
* ```
*/
var Composer = class {
	constructor(options = {}) {
		this.doc = null;
		this.atDirectives = false;
		this.prelude = [];
		this.errors = [];
		this.warnings = [];
		this.onError = (source, code, message, warning) => {
			const pos = getErrorPos(source);
			if (warning) this.warnings.push(new YAMLWarning(pos, code, message));
			else this.errors.push(new YAMLParseError(pos, code, message));
		};
		this.directives = new Directives({ version: options.version || "1.2" });
		this.options = options;
	}
	decorate(doc, afterDoc) {
		const { comment, afterEmptyLine } = parsePrelude(this.prelude);
		if (comment) {
			const dc = doc.contents;
			if (afterDoc) doc.comment = doc.comment ? `${doc.comment}\n${comment}` : comment;
			else if (afterEmptyLine || doc.directives.docStart || !dc) doc.commentBefore = comment;
			else if (isCollection$1(dc) && !dc.flow && dc.items.length > 0) {
				let it = dc.items[0];
				if (isPair(it)) it = it.key;
				const cb = it.commentBefore;
				it.commentBefore = cb ? `${comment}\n${cb}` : comment;
			} else {
				const cb = dc.commentBefore;
				dc.commentBefore = cb ? `${comment}\n${cb}` : comment;
			}
		}
		if (afterDoc) {
			for (let i = 0; i < this.errors.length; ++i) doc.errors.push(this.errors[i]);
			for (let i = 0; i < this.warnings.length; ++i) doc.warnings.push(this.warnings[i]);
		} else {
			doc.errors = this.errors;
			doc.warnings = this.warnings;
		}
		this.prelude = [];
		this.errors = [];
		this.warnings = [];
	}
	/**
	* Current stream status information.
	*
	* Mostly useful at the end of input for an empty stream.
	*/
	streamInfo() {
		return {
			comment: parsePrelude(this.prelude).comment,
			directives: this.directives,
			errors: this.errors,
			warnings: this.warnings
		};
	}
	/**
	* Compose tokens into documents.
	*
	* @param forceDoc - If the stream contains no document, still emit a final document including any comments and directives that would be applied to a subsequent document.
	* @param endOffset - Should be set if `forceDoc` is also set, to set the document range end and to indicate errors correctly.
	*/
	*compose(tokens, forceDoc = false, endOffset = -1) {
		for (const token of tokens) yield* this.next(token);
		yield* this.end(forceDoc, endOffset);
	}
	/** Advance the composer by one CST token. */
	*next(token) {
		switch (token.type) {
			case "directive":
				this.directives.add(token.source, (offset, message, warning) => {
					const pos = getErrorPos(token);
					pos[0] += offset;
					this.onError(pos, "BAD_DIRECTIVE", message, warning);
				});
				this.prelude.push(token.source);
				this.atDirectives = true;
				break;
			case "document": {
				const doc = composeDoc(this.options, this.directives, token, this.onError);
				if (this.atDirectives && !doc.directives.docStart) this.onError(token, "MISSING_CHAR", "Missing directives-end/doc-start indicator line");
				this.decorate(doc, false);
				if (this.doc) yield this.doc;
				this.doc = doc;
				this.atDirectives = false;
				break;
			}
			case "byte-order-mark":
			case "space": break;
			case "comment":
			case "newline":
				this.prelude.push(token.source);
				break;
			case "error": {
				const msg = token.source ? `${token.message}: ${JSON.stringify(token.source)}` : token.message;
				const error = new YAMLParseError(getErrorPos(token), "UNEXPECTED_TOKEN", msg);
				if (this.atDirectives || !this.doc) this.errors.push(error);
				else this.doc.errors.push(error);
				break;
			}
			case "doc-end": {
				if (!this.doc) {
					this.errors.push(new YAMLParseError(getErrorPos(token), "UNEXPECTED_TOKEN", "Unexpected doc-end without preceding document"));
					break;
				}
				this.doc.directives.docEnd = true;
				const end = resolveEnd(token.end, token.offset + token.source.length, this.doc.options.strict, this.onError);
				this.decorate(this.doc, true);
				if (end.comment) {
					const dc = this.doc.comment;
					this.doc.comment = dc ? `${dc}\n${end.comment}` : end.comment;
				}
				this.doc.range[2] = end.offset;
				break;
			}
			default: this.errors.push(new YAMLParseError(getErrorPos(token), "UNEXPECTED_TOKEN", `Unsupported token ${token.type}`));
		}
	}
	/**
	* Call at end of input to yield any remaining document.
	*
	* @param forceDoc - If the stream contains no document, still emit a final document including any comments and directives that would be applied to a subsequent document.
	* @param endOffset - Should be set if `forceDoc` is also set, to set the document range end and to indicate errors correctly.
	*/
	*end(forceDoc = false, endOffset = -1) {
		if (this.doc) {
			this.decorate(this.doc, true);
			yield this.doc;
			this.doc = null;
		} else if (forceDoc) {
			const doc = new Document(void 0, Object.assign({ _directives: this.directives }, this.options));
			if (this.atDirectives) this.onError(endOffset, "MISSING_CHAR", "Missing directives-end indicator line");
			doc.range = [
				0,
				endOffset,
				endOffset
			];
			this.decorate(doc, false);
			yield doc;
		}
	}
};
//#endregion
//#region node_modules/yaml/browser/dist/parse/cst-scalar.js
function resolveAsScalar(token, strict = true, onError) {
	if (token) {
		const _onError = (pos, code, message) => {
			const offset = typeof pos === "number" ? pos : Array.isArray(pos) ? pos[0] : pos.offset;
			if (onError) onError(offset, code, message);
			else throw new YAMLParseError([offset, offset + 1], code, message);
		};
		switch (token.type) {
			case "scalar":
			case "single-quoted-scalar":
			case "double-quoted-scalar": return resolveFlowScalar(token, strict, _onError);
			case "block-scalar": return resolveBlockScalar({ options: { strict } }, token, _onError);
		}
	}
	return null;
}
/**
* Create a new scalar token with `value`
*
* Values that represent an actual string but may be parsed as a different type should use a `type` other than `'PLAIN'`,
* as this function does not support any schema operations and won't check for such conflicts.
*
* @param value The string representation of the value, which will have its content properly indented.
* @param context.end Comments and whitespace after the end of the value, or after the block scalar header. If undefined, a newline will be added.
* @param context.implicitKey Being within an implicit key may affect the resolved type of the token's value.
* @param context.indent The indent level of the token.
* @param context.inFlow Is this scalar within a flow collection? This may affect the resolved type of the token's value.
* @param context.offset The offset position of the token.
* @param context.type The preferred type of the scalar token. If undefined, the previous type of the `token` will be used, defaulting to `'PLAIN'`.
*/
function createScalarToken(value, context) {
	const { implicitKey = false, indent, inFlow = false, offset = -1, type = "PLAIN" } = context;
	const source = stringifyString({
		type,
		value
	}, {
		implicitKey,
		indent: indent > 0 ? " ".repeat(indent) : "",
		inFlow,
		options: {
			blockQuote: true,
			lineWidth: -1
		}
	});
	const end = context.end ?? [{
		type: "newline",
		offset: -1,
		indent,
		source: "\n"
	}];
	switch (source[0]) {
		case "|":
		case ">": {
			const he = source.indexOf("\n");
			const head = source.substring(0, he);
			const body = source.substring(he + 1) + "\n";
			const props = [{
				type: "block-scalar-header",
				offset,
				indent,
				source: head
			}];
			if (!addEndtoBlockProps(props, end)) props.push({
				type: "newline",
				offset: -1,
				indent,
				source: "\n"
			});
			return {
				type: "block-scalar",
				offset,
				indent,
				props,
				source: body
			};
		}
		case "\"": return {
			type: "double-quoted-scalar",
			offset,
			indent,
			source,
			end
		};
		case "'": return {
			type: "single-quoted-scalar",
			offset,
			indent,
			source,
			end
		};
		default: return {
			type: "scalar",
			offset,
			indent,
			source,
			end
		};
	}
}
/**
* Set the value of `token` to the given string `value`, overwriting any previous contents and type that it may have.
*
* Best efforts are made to retain any comments previously associated with the `token`,
* though all contents within a collection's `items` will be overwritten.
*
* Values that represent an actual string but may be parsed as a different type should use a `type` other than `'PLAIN'`,
* as this function does not support any schema operations and won't check for such conflicts.
*
* @param token Any token. If it does not include an `indent` value, the value will be stringified as if it were an implicit key.
* @param value The string representation of the value, which will have its content properly indented.
* @param context.afterKey In most cases, values after a key should have an additional level of indentation.
* @param context.implicitKey Being within an implicit key may affect the resolved type of the token's value.
* @param context.inFlow Being within a flow collection may affect the resolved type of the token's value.
* @param context.type The preferred type of the scalar token. If undefined, the previous type of the `token` will be used, defaulting to `'PLAIN'`.
*/
function setScalarValue(token, value, context = {}) {
	let { afterKey = false, implicitKey = false, inFlow = false, type } = context;
	let indent = "indent" in token ? token.indent : null;
	if (afterKey && typeof indent === "number") indent += 2;
	if (!type) switch (token.type) {
		case "single-quoted-scalar":
			type = "QUOTE_SINGLE";
			break;
		case "double-quoted-scalar":
			type = "QUOTE_DOUBLE";
			break;
		case "block-scalar": {
			const header = token.props[0];
			if (header.type !== "block-scalar-header") throw new Error("Invalid block scalar header");
			type = header.source[0] === ">" ? "BLOCK_FOLDED" : "BLOCK_LITERAL";
			break;
		}
		default: type = "PLAIN";
	}
	const source = stringifyString({
		type,
		value
	}, {
		implicitKey: implicitKey || indent === null,
		indent: indent !== null && indent > 0 ? " ".repeat(indent) : "",
		inFlow,
		options: {
			blockQuote: true,
			lineWidth: -1
		}
	});
	switch (source[0]) {
		case "|":
		case ">":
			setBlockScalarValue(token, source);
			break;
		case "\"":
			setFlowScalarValue(token, source, "double-quoted-scalar");
			break;
		case "'":
			setFlowScalarValue(token, source, "single-quoted-scalar");
			break;
		default: setFlowScalarValue(token, source, "scalar");
	}
}
function setBlockScalarValue(token, source) {
	const he = source.indexOf("\n");
	const head = source.substring(0, he);
	const body = source.substring(he + 1) + "\n";
	if (token.type === "block-scalar") {
		const header = token.props[0];
		if (header.type !== "block-scalar-header") throw new Error("Invalid block scalar header");
		header.source = head;
		token.source = body;
	} else {
		const { offset } = token;
		const indent = "indent" in token ? token.indent : -1;
		const props = [{
			type: "block-scalar-header",
			offset,
			indent,
			source: head
		}];
		if (!addEndtoBlockProps(props, "end" in token ? token.end : void 0)) props.push({
			type: "newline",
			offset: -1,
			indent,
			source: "\n"
		});
		for (const key of Object.keys(token)) if (key !== "type" && key !== "offset") delete token[key];
		Object.assign(token, {
			type: "block-scalar",
			indent,
			props,
			source: body
		});
	}
}
/** @returns `true` if last token is a newline */
function addEndtoBlockProps(props, end) {
	if (end) for (const st of end) switch (st.type) {
		case "space":
		case "comment":
			props.push(st);
			break;
		case "newline":
			props.push(st);
			return true;
	}
	return false;
}
function setFlowScalarValue(token, source, type) {
	switch (token.type) {
		case "scalar":
		case "double-quoted-scalar":
		case "single-quoted-scalar":
			token.type = type;
			token.source = source;
			break;
		case "block-scalar": {
			const end = token.props.slice(1);
			let oa = source.length;
			if (token.props[0].type === "block-scalar-header") oa -= token.props[0].source.length;
			for (const tok of end) tok.offset += oa;
			delete token.props;
			Object.assign(token, {
				type,
				source,
				end
			});
			break;
		}
		case "block-map":
		case "block-seq": {
			const nl = {
				type: "newline",
				offset: token.offset + source.length,
				indent: token.indent,
				source: "\n"
			};
			delete token.items;
			Object.assign(token, {
				type,
				source,
				end: [nl]
			});
			break;
		}
		default: {
			const indent = "indent" in token ? token.indent : -1;
			const end = "end" in token && Array.isArray(token.end) ? token.end.filter((st) => st.type === "space" || st.type === "comment" || st.type === "newline") : [];
			for (const key of Object.keys(token)) if (key !== "type" && key !== "offset") delete token[key];
			Object.assign(token, {
				type,
				indent,
				source,
				end
			});
		}
	}
}
//#endregion
//#region node_modules/yaml/browser/dist/parse/cst-stringify.js
/**
* Stringify a CST document, token, or collection item
*
* Fair warning: This applies no validation whatsoever, and
* simply concatenates the sources in their logical order.
*/
const stringify$1 = (cst) => "type" in cst ? stringifyToken(cst) : stringifyItem(cst);
function stringifyToken(token) {
	switch (token.type) {
		case "block-scalar": {
			let res = "";
			for (const tok of token.props) res += stringifyToken(tok);
			return res + token.source;
		}
		case "block-map":
		case "block-seq": {
			let res = "";
			for (const item of token.items) res += stringifyItem(item);
			return res;
		}
		case "flow-collection": {
			let res = token.start.source;
			for (const item of token.items) res += stringifyItem(item);
			for (const st of token.end) res += st.source;
			return res;
		}
		case "document": {
			let res = stringifyItem(token);
			if (token.end) for (const st of token.end) res += st.source;
			return res;
		}
		default: {
			let res = token.source;
			if ("end" in token && token.end) for (const st of token.end) res += st.source;
			return res;
		}
	}
}
function stringifyItem({ start, key, sep, value }) {
	let res = "";
	for (const st of start) res += st.source;
	if (key) res += stringifyToken(key);
	if (sep) for (const st of sep) res += st.source;
	if (value) res += stringifyToken(value);
	return res;
}
//#endregion
//#region node_modules/yaml/browser/dist/parse/cst-visit.js
const BREAK = Symbol("break visit");
const SKIP = Symbol("skip children");
const REMOVE = Symbol("remove item");
/**
* Apply a visitor to a CST document or item.
*
* Walks through the tree (depth-first) starting from the root, calling a
* `visitor` function with two arguments when entering each item:
*   - `item`: The current item, which included the following members:
*     - `start: SourceToken[]` – Source tokens before the key or value,
*       possibly including its anchor or tag.
*     - `key?: Token | null` – Set for pair values. May then be `null`, if
*       the key before the `:` separator is empty.
*     - `sep?: SourceToken[]` – Source tokens between the key and the value,
*       which should include the `:` map value indicator if `value` is set.
*     - `value?: Token` – The value of a sequence item, or of a map pair.
*   - `path`: The steps from the root to the current node, as an array of
*     `['key' | 'value', number]` tuples.
*
* The return value of the visitor may be used to control the traversal:
*   - `undefined` (default): Do nothing and continue
*   - `visit.SKIP`: Do not visit the children of this token, continue with
*      next sibling
*   - `visit.BREAK`: Terminate traversal completely
*   - `visit.REMOVE`: Remove the current item, then continue with the next one
*   - `number`: Set the index of the next step. This is useful especially if
*     the index of the current token has changed.
*   - `function`: Define the next visitor for this item. After the original
*     visitor is called on item entry, next visitors are called after handling
*     a non-empty `key` and when exiting the item.
*/
function visit(cst, visitor) {
	if ("type" in cst && cst.type === "document") cst = {
		start: cst.start,
		value: cst.value
	};
	_visit(Object.freeze([]), cst, visitor);
}
/** Terminate visit traversal completely */
visit.BREAK = BREAK;
/** Do not visit the children of the current item */
visit.SKIP = SKIP;
/** Remove the current item */
visit.REMOVE = REMOVE;
/** Find the item at `path` from `cst` as the root */
visit.itemAtPath = (cst, path) => {
	let item = cst;
	for (const [field, index] of path) {
		const tok = item?.[field];
		if (tok && "items" in tok) item = tok.items[index];
		else return void 0;
	}
	return item;
};
/**
* Get the immediate parent collection of the item at `path` from `cst` as the root.
*
* Throws an error if the collection is not found, which should never happen if the item itself exists.
*/
visit.parentCollection = (cst, path) => {
	const parent = visit.itemAtPath(cst, path.slice(0, -1));
	const field = path[path.length - 1][0];
	const coll = parent?.[field];
	if (coll && "items" in coll) return coll;
	throw new Error("Parent collection not found");
};
function _visit(path, item, visitor) {
	let ctrl = visitor(item, path);
	if (typeof ctrl === "symbol") return ctrl;
	for (const field of ["key", "value"]) {
		const token = item[field];
		if (token && "items" in token) {
			for (let i = 0; i < token.items.length; ++i) {
				const ci = _visit(Object.freeze(path.concat([[field, i]])), token.items[i], visitor);
				if (typeof ci === "number") i = ci - 1;
				else if (ci === BREAK) return BREAK;
				else if (ci === REMOVE) {
					token.items.splice(i, 1);
					i -= 1;
				}
			}
			if (typeof ctrl === "function" && field === "key") ctrl = ctrl(item, path);
		}
	}
	return typeof ctrl === "function" ? ctrl(item, path) : ctrl;
}
//#endregion
//#region node_modules/yaml/browser/dist/parse/cst.js
var cst_exports = /* @__PURE__ */ __exportAll({
	BOM: () => "﻿",
	DOCUMENT: () => "",
	FLOW_END: () => "",
	SCALAR: () => "",
	createScalarToken: () => createScalarToken,
	isCollection: () => isCollection,
	isScalar: () => isScalar,
	prettyToken: () => prettyToken,
	resolveAsScalar: () => resolveAsScalar,
	setScalarValue: () => setScalarValue,
	stringify: () => stringify$1,
	tokenType: () => tokenType,
	visit: () => visit
});
/** @returns `true` if `token` is a flow or block collection */
const isCollection = (token) => !!token && "items" in token;
/** @returns `true` if `token` is a flow or block scalar; not an alias */
const isScalar = (token) => !!token && (token.type === "scalar" || token.type === "single-quoted-scalar" || token.type === "double-quoted-scalar" || token.type === "block-scalar");
/* istanbul ignore next */
/** Get a printable representation of a lexer token */
function prettyToken(token) {
	switch (token) {
		case "﻿": return "<BOM>";
		case "": return "<DOC>";
		case "": return "<FLOW_END>";
		case "": return "<SCALAR>";
		default: return JSON.stringify(token);
	}
}
/** Identify the type of a lexer token. May return `null` for unknown tokens. */
function tokenType(source) {
	switch (source) {
		case "﻿": return "byte-order-mark";
		case "": return "doc-mode";
		case "": return "flow-error-end";
		case "": return "scalar";
		case "---": return "doc-start";
		case "...": return "doc-end";
		case "":
		case "\n":
		case "\r\n": return "newline";
		case "-": return "seq-item-ind";
		case "?": return "explicit-key-ind";
		case ":": return "map-value-ind";
		case "{": return "flow-map-start";
		case "}": return "flow-map-end";
		case "[": return "flow-seq-start";
		case "]": return "flow-seq-end";
		case ",": return "comma";
	}
	switch (source[0]) {
		case " ":
		case "	": return "space";
		case "#": return "comment";
		case "%": return "directive-line";
		case "*": return "alias";
		case "&": return "anchor";
		case "!": return "tag";
		case "'": return "single-quoted-scalar";
		case "\"": return "double-quoted-scalar";
		case "|":
		case ">": return "block-scalar-header";
	}
	return null;
}
//#endregion
//#region node_modules/yaml/browser/dist/parse/lexer.js
function isEmpty(ch) {
	switch (ch) {
		case void 0:
		case " ":
		case "\n":
		case "\r":
		case "	": return true;
		default: return false;
	}
}
const hexDigits = /* @__PURE__ */ new Set("0123456789ABCDEFabcdef");
const tagChars = /* @__PURE__ */ new Set("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz-#;/?:@&=+$_.!~*'()");
const flowIndicatorChars = /* @__PURE__ */ new Set(",[]{}");
const invalidAnchorChars = /* @__PURE__ */ new Set(" ,[]{}\n\r	");
const isNotAnchorChar = (ch) => !ch || invalidAnchorChars.has(ch);
/**
* Splits an input string into lexical tokens, i.e. smaller strings that are
* easily identifiable by `tokens.tokenType()`.
*
* Lexing starts always in a "stream" context. Incomplete input may be buffered
* until a complete token can be emitted.
*
* In addition to slices of the original input, the following control characters
* may also be emitted:
*
* - `\x02` (Start of Text): A document starts with the next token
* - `\x18` (Cancel): Unexpected end of flow-mode (indicates an error)
* - `\x1f` (Unit Separator): Next token is a scalar value
* - `\u{FEFF}` (Byte order mark): Emitted separately outside documents
*/
var Lexer = class {
	constructor() {
		/**
		* Flag indicating whether the end of the current buffer marks the end of
		* all input
		*/
		this.atEnd = false;
		/**
		* Explicit indent set in block scalar header, as an offset from the current
		* minimum indent, so e.g. set to 1 from a header `|2+`. Set to -1 if not
		* explicitly set.
		*/
		this.blockScalarIndent = -1;
		/**
		* Block scalars that include a + (keep) chomping indicator in their header
		* include trailing empty lines, which are otherwise excluded from the
		* scalar's contents.
		*/
		this.blockScalarKeep = false;
		/** Current input */
		this.buffer = "";
		/**
		* Flag noting whether the map value indicator : can immediately follow this
		* node within a flow context.
		*/
		this.flowKey = false;
		/** Count of surrounding flow collection levels. */
		this.flowLevel = 0;
		/**
		* Minimum level of indentation required for next lines to be parsed as a
		* part of the current scalar value.
		*/
		this.indentNext = 0;
		/** Indentation level of the current line. */
		this.indentValue = 0;
		/** Position of the next \n character. */
		this.lineEndPos = null;
		/** Stores the state of the lexer if reaching the end of incpomplete input */
		this.next = null;
		/** A pointer to `buffer`; the current position of the lexer. */
		this.pos = 0;
	}
	/**
	* Generate YAML tokens from the `source` string. If `incomplete`,
	* a part of the last line may be left as a buffer for the next call.
	*
	* @returns A generator of lexical tokens
	*/
	*lex(source, incomplete = false) {
		if (source) {
			if (typeof source !== "string") throw TypeError("source is not a string");
			this.buffer = this.buffer ? this.buffer + source : source;
			this.lineEndPos = null;
		}
		this.atEnd = !incomplete;
		let next = this.next ?? "stream";
		while (next && (incomplete || this.hasChars(1))) next = yield* this.parseNext(next);
	}
	atLineEnd() {
		let i = this.pos;
		let ch = this.buffer[i];
		while (ch === " " || ch === "	") ch = this.buffer[++i];
		if (!ch || ch === "#" || ch === "\n") return true;
		if (ch === "\r") return this.buffer[i + 1] === "\n";
		return false;
	}
	charAt(n) {
		return this.buffer[this.pos + n];
	}
	continueScalar(offset) {
		let ch = this.buffer[offset];
		if (this.indentNext > 0) {
			let indent = 0;
			while (ch === " ") ch = this.buffer[++indent + offset];
			if (ch === "\r") {
				const next = this.buffer[indent + offset + 1];
				if (next === "\n" || !next && !this.atEnd) return offset + indent + 1;
			}
			return ch === "\n" || indent >= this.indentNext || !ch && !this.atEnd ? offset + indent : -1;
		}
		if (ch === "-" || ch === ".") {
			const dt = this.buffer.substr(offset, 3);
			if ((dt === "---" || dt === "...") && isEmpty(this.buffer[offset + 3])) return -1;
		}
		return offset;
	}
	getLine() {
		let end = this.lineEndPos;
		if (typeof end !== "number" || end !== -1 && end < this.pos) {
			end = this.buffer.indexOf("\n", this.pos);
			this.lineEndPos = end;
		}
		if (end === -1) return this.atEnd ? this.buffer.substring(this.pos) : null;
		if (this.buffer[end - 1] === "\r") end -= 1;
		return this.buffer.substring(this.pos, end);
	}
	hasChars(n) {
		return this.pos + n <= this.buffer.length;
	}
	setNext(state) {
		this.buffer = this.buffer.substring(this.pos);
		this.pos = 0;
		this.lineEndPos = null;
		this.next = state;
		return null;
	}
	peek(n) {
		return this.buffer.substr(this.pos, n);
	}
	*parseNext(next) {
		switch (next) {
			case "stream": return yield* this.parseStream();
			case "line-start": return yield* this.parseLineStart();
			case "block-start": return yield* this.parseBlockStart();
			case "doc": return yield* this.parseDocument();
			case "flow": return yield* this.parseFlowCollection();
			case "quoted-scalar": return yield* this.parseQuotedScalar();
			case "block-scalar": return yield* this.parseBlockScalar();
			case "plain-scalar": return yield* this.parsePlainScalar();
		}
	}
	*parseStream() {
		let line = this.getLine();
		if (line === null) return this.setNext("stream");
		if (line[0] === "﻿") {
			yield* this.pushCount(1);
			line = line.substring(1);
		}
		if (line[0] === "%") {
			let dirEnd = line.length;
			let cs = line.indexOf("#");
			while (cs !== -1) {
				const ch = line[cs - 1];
				if (ch === " " || ch === "	") {
					dirEnd = cs - 1;
					break;
				} else cs = line.indexOf("#", cs + 1);
			}
			while (true) {
				const ch = line[dirEnd - 1];
				if (ch === " " || ch === "	") dirEnd -= 1;
				else break;
			}
			const n = (yield* this.pushCount(dirEnd)) + (yield* this.pushSpaces(true));
			yield* this.pushCount(line.length - n);
			this.pushNewline();
			return "stream";
		}
		if (this.atLineEnd()) {
			const sp = yield* this.pushSpaces(true);
			yield* this.pushCount(line.length - sp);
			yield* this.pushNewline();
			return "stream";
		}
		yield "";
		return yield* this.parseLineStart();
	}
	*parseLineStart() {
		const ch = this.charAt(0);
		if (!ch && !this.atEnd) return this.setNext("line-start");
		if (ch === "-" || ch === ".") {
			if (!this.atEnd && !this.hasChars(4)) return this.setNext("line-start");
			const s = this.peek(3);
			if ((s === "---" || s === "...") && isEmpty(this.charAt(3))) {
				yield* this.pushCount(3);
				this.indentValue = 0;
				this.indentNext = 0;
				return s === "---" ? "doc" : "stream";
			}
		}
		this.indentValue = yield* this.pushSpaces(false);
		if (this.indentNext > this.indentValue && !isEmpty(this.charAt(1))) this.indentNext = this.indentValue;
		return yield* this.parseBlockStart();
	}
	*parseBlockStart() {
		const [ch0, ch1] = this.peek(2);
		if (!ch1 && !this.atEnd) return this.setNext("block-start");
		if ((ch0 === "-" || ch0 === "?" || ch0 === ":") && isEmpty(ch1)) {
			const n = (yield* this.pushCount(1)) + (yield* this.pushSpaces(true));
			this.indentNext = this.indentValue + 1;
			this.indentValue += n;
			return "block-start";
		}
		return "doc";
	}
	*parseDocument() {
		yield* this.pushSpaces(true);
		const line = this.getLine();
		if (line === null) return this.setNext("doc");
		let n = yield* this.pushIndicators();
		switch (line[n]) {
			case "#": yield* this.pushCount(line.length - n);
			case void 0:
				yield* this.pushNewline();
				return yield* this.parseLineStart();
			case "{":
			case "[":
				yield* this.pushCount(1);
				this.flowKey = false;
				this.flowLevel = 1;
				return "flow";
			case "}":
			case "]":
				yield* this.pushCount(1);
				return "doc";
			case "*":
				yield* this.pushUntil(isNotAnchorChar);
				return "doc";
			case "\"":
			case "'": return yield* this.parseQuotedScalar();
			case "|":
			case ">":
				n += yield* this.parseBlockScalarHeader();
				n += yield* this.pushSpaces(true);
				yield* this.pushCount(line.length - n);
				yield* this.pushNewline();
				return yield* this.parseBlockScalar();
			default: return yield* this.parsePlainScalar();
		}
	}
	*parseFlowCollection() {
		let nl, sp;
		let indent = -1;
		do {
			nl = yield* this.pushNewline();
			if (nl > 0) {
				sp = yield* this.pushSpaces(false);
				this.indentValue = indent = sp;
			} else sp = 0;
			sp += yield* this.pushSpaces(true);
		} while (nl + sp > 0);
		const line = this.getLine();
		if (line === null) return this.setNext("flow");
		if (indent !== -1 && indent < this.indentNext && line[0] !== "#" || indent === 0 && (line.startsWith("---") || line.startsWith("...")) && isEmpty(line[3])) {
			if (!(indent === this.indentNext - 1 && this.flowLevel === 1 && (line[0] === "]" || line[0] === "}"))) {
				this.flowLevel = 0;
				yield "";
				return yield* this.parseLineStart();
			}
		}
		let n = 0;
		while (line[n] === ",") {
			n += yield* this.pushCount(1);
			n += yield* this.pushSpaces(true);
			this.flowKey = false;
		}
		n += yield* this.pushIndicators();
		switch (line[n]) {
			case void 0: return "flow";
			case "#":
				yield* this.pushCount(line.length - n);
				return "flow";
			case "{":
			case "[":
				yield* this.pushCount(1);
				this.flowKey = false;
				this.flowLevel += 1;
				return "flow";
			case "}":
			case "]":
				yield* this.pushCount(1);
				this.flowKey = true;
				this.flowLevel -= 1;
				return this.flowLevel ? "flow" : "doc";
			case "*":
				yield* this.pushUntil(isNotAnchorChar);
				return "flow";
			case "\"":
			case "'":
				this.flowKey = true;
				return yield* this.parseQuotedScalar();
			case ":": {
				const next = this.charAt(1);
				if (this.flowKey || isEmpty(next) || next === ",") {
					this.flowKey = false;
					yield* this.pushCount(1);
					yield* this.pushSpaces(true);
					return "flow";
				}
			}
			default:
				this.flowKey = false;
				return yield* this.parsePlainScalar();
		}
	}
	*parseQuotedScalar() {
		const quote = this.charAt(0);
		let end = this.buffer.indexOf(quote, this.pos + 1);
		if (quote === "'") while (end !== -1 && this.buffer[end + 1] === "'") end = this.buffer.indexOf("'", end + 2);
		else while (end !== -1) {
			let n = 0;
			while (this.buffer[end - 1 - n] === "\\") n += 1;
			if (n % 2 === 0) break;
			end = this.buffer.indexOf("\"", end + 1);
		}
		const qb = this.buffer.substring(0, end);
		let nl = qb.indexOf("\n", this.pos);
		if (nl !== -1) {
			while (nl !== -1) {
				const cs = this.continueScalar(nl + 1);
				if (cs === -1) break;
				nl = qb.indexOf("\n", cs);
			}
			if (nl !== -1) end = nl - (qb[nl - 1] === "\r" ? 2 : 1);
		}
		if (end === -1) {
			if (!this.atEnd) return this.setNext("quoted-scalar");
			end = this.buffer.length;
		}
		yield* this.pushToIndex(end + 1, false);
		return this.flowLevel ? "flow" : "doc";
	}
	*parseBlockScalarHeader() {
		this.blockScalarIndent = -1;
		this.blockScalarKeep = false;
		let i = this.pos;
		while (true) {
			const ch = this.buffer[++i];
			if (ch === "+") this.blockScalarKeep = true;
			else if (ch > "0" && ch <= "9") this.blockScalarIndent = Number(ch) - 1;
			else if (ch !== "-") break;
		}
		return yield* this.pushUntil((ch) => isEmpty(ch) || ch === "#");
	}
	*parseBlockScalar() {
		let nl = this.pos - 1;
		let indent = 0;
		let ch;
		loop: for (let i = this.pos; ch = this.buffer[i]; ++i) switch (ch) {
			case " ":
				indent += 1;
				break;
			case "\n":
				nl = i;
				indent = 0;
				break;
			case "\r": {
				const next = this.buffer[i + 1];
				if (!next && !this.atEnd) return this.setNext("block-scalar");
				if (next === "\n") break;
			}
			default: break loop;
		}
		if (!ch && !this.atEnd) return this.setNext("block-scalar");
		if (indent >= this.indentNext) {
			if (this.blockScalarIndent === -1) this.indentNext = indent;
			else this.indentNext = this.blockScalarIndent + (this.indentNext === 0 ? 1 : this.indentNext);
			do {
				const cs = this.continueScalar(nl + 1);
				if (cs === -1) break;
				nl = this.buffer.indexOf("\n", cs);
			} while (nl !== -1);
			if (nl === -1) {
				if (!this.atEnd) return this.setNext("block-scalar");
				nl = this.buffer.length;
			}
		}
		let i = nl + 1;
		ch = this.buffer[i];
		while (ch === " ") ch = this.buffer[++i];
		if (ch === "	") {
			while (ch === "	" || ch === " " || ch === "\r" || ch === "\n") ch = this.buffer[++i];
			nl = i - 1;
		} else if (!this.blockScalarKeep) do {
			let i = nl - 1;
			let ch = this.buffer[i];
			if (ch === "\r") ch = this.buffer[--i];
			const lastChar = i;
			while (ch === " ") ch = this.buffer[--i];
			if (ch === "\n" && i >= this.pos && i + 1 + indent > lastChar) nl = i;
			else break;
		} while (true);
		yield "";
		yield* this.pushToIndex(nl + 1, true);
		return yield* this.parseLineStart();
	}
	*parsePlainScalar() {
		const inFlow = this.flowLevel > 0;
		let end = this.pos - 1;
		let i = this.pos - 1;
		let ch;
		while (ch = this.buffer[++i]) if (ch === ":") {
			const next = this.buffer[i + 1];
			if (isEmpty(next) || inFlow && flowIndicatorChars.has(next)) break;
			end = i;
		} else if (isEmpty(ch)) {
			let next = this.buffer[i + 1];
			if (ch === "\r") if (next === "\n") {
				i += 1;
				ch = "\n";
				next = this.buffer[i + 1];
			} else end = i;
			if (next === "#" || inFlow && flowIndicatorChars.has(next)) break;
			if (ch === "\n") {
				const cs = this.continueScalar(i + 1);
				if (cs === -1) break;
				i = Math.max(i, cs - 2);
			}
		} else {
			if (inFlow && flowIndicatorChars.has(ch)) break;
			end = i;
		}
		if (!ch && !this.atEnd) return this.setNext("plain-scalar");
		yield "";
		yield* this.pushToIndex(end + 1, true);
		return inFlow ? "flow" : "doc";
	}
	*pushCount(n) {
		if (n > 0) {
			yield this.buffer.substr(this.pos, n);
			this.pos += n;
			return n;
		}
		return 0;
	}
	*pushToIndex(i, allowEmpty) {
		const s = this.buffer.slice(this.pos, i);
		if (s) {
			yield s;
			this.pos += s.length;
			return s.length;
		} else if (allowEmpty) yield "";
		return 0;
	}
	*pushIndicators() {
		let n = 0;
		loop: while (true) {
			switch (this.charAt(0)) {
				case "!":
					n += yield* this.pushTag();
					n += yield* this.pushSpaces(true);
					continue loop;
				case "&":
					n += yield* this.pushUntil(isNotAnchorChar);
					n += yield* this.pushSpaces(true);
					continue loop;
				case "-":
				case "?":
				case ":": {
					const inFlow = this.flowLevel > 0;
					const ch1 = this.charAt(1);
					if (isEmpty(ch1) || inFlow && flowIndicatorChars.has(ch1)) {
						if (!inFlow) this.indentNext = this.indentValue + 1;
						else if (this.flowKey) this.flowKey = false;
						n += yield* this.pushCount(1);
						n += yield* this.pushSpaces(true);
						continue loop;
					}
				}
			}
			break loop;
		}
		return n;
	}
	*pushTag() {
		if (this.charAt(1) === "<") {
			let i = this.pos + 2;
			let ch = this.buffer[i];
			while (!isEmpty(ch) && ch !== ">") ch = this.buffer[++i];
			return yield* this.pushToIndex(ch === ">" ? i + 1 : i, false);
		} else {
			let i = this.pos + 1;
			let ch = this.buffer[i];
			while (ch) if (tagChars.has(ch)) ch = this.buffer[++i];
			else if (ch === "%" && hexDigits.has(this.buffer[i + 1]) && hexDigits.has(this.buffer[i + 2])) ch = this.buffer[i += 3];
			else break;
			return yield* this.pushToIndex(i, false);
		}
	}
	*pushNewline() {
		const ch = this.buffer[this.pos];
		if (ch === "\n") return yield* this.pushCount(1);
		else if (ch === "\r" && this.charAt(1) === "\n") return yield* this.pushCount(2);
		else return 0;
	}
	*pushSpaces(allowTabs) {
		let i = this.pos - 1;
		let ch;
		do
			ch = this.buffer[++i];
		while (ch === " " || allowTabs && ch === "	");
		const n = i - this.pos;
		if (n > 0) {
			yield this.buffer.substr(this.pos, n);
			this.pos = i;
		}
		return n;
	}
	*pushUntil(test) {
		let i = this.pos;
		let ch = this.buffer[i];
		while (!test(ch)) ch = this.buffer[++i];
		return yield* this.pushToIndex(i, false);
	}
};
//#endregion
//#region node_modules/yaml/browser/dist/parse/line-counter.js
/**
* Tracks newlines during parsing in order to provide an efficient API for
* determining the one-indexed `{ line, col }` position for any offset
* within the input.
*/
var LineCounter = class {
	constructor() {
		this.lineStarts = [];
		/**
		* Should be called in ascending order. Otherwise, call
		* `lineCounter.lineStarts.sort()` before calling `linePos()`.
		*/
		this.addNewLine = (offset) => this.lineStarts.push(offset);
		/**
		* Performs a binary search and returns the 1-indexed { line, col }
		* position of `offset`. If `line === 0`, `addNewLine` has never been
		* called or `offset` is before the first known newline.
		*/
		this.linePos = (offset) => {
			let low = 0;
			let high = this.lineStarts.length;
			while (low < high) {
				const mid = low + high >> 1;
				if (this.lineStarts[mid] < offset) low = mid + 1;
				else high = mid;
			}
			if (this.lineStarts[low] === offset) return {
				line: low + 1,
				col: 1
			};
			if (low === 0) return {
				line: 0,
				col: offset
			};
			const start = this.lineStarts[low - 1];
			return {
				line: low,
				col: offset - start + 1
			};
		};
	}
};
//#endregion
//#region node_modules/yaml/browser/dist/parse/parser.js
function includesToken(list, type) {
	for (let i = 0; i < list.length; ++i) if (list[i].type === type) return true;
	return false;
}
function findNonEmptyIndex(list) {
	for (let i = 0; i < list.length; ++i) switch (list[i].type) {
		case "space":
		case "comment":
		case "newline": break;
		default: return i;
	}
	return -1;
}
function isFlowToken(token) {
	switch (token?.type) {
		case "alias":
		case "scalar":
		case "single-quoted-scalar":
		case "double-quoted-scalar":
		case "flow-collection": return true;
		default: return false;
	}
}
function getPrevProps(parent) {
	switch (parent.type) {
		case "document": return parent.start;
		case "block-map": {
			const it = parent.items[parent.items.length - 1];
			return it.sep ?? it.start;
		}
		case "block-seq": return parent.items[parent.items.length - 1].start;
		/* istanbul ignore next should not happen */
		default: return [];
	}
}
/** Note: May modify input array */
function getFirstKeyStartProps(prev) {
	if (prev.length === 0) return [];
	let i = prev.length;
	loop: while (--i >= 0) switch (prev[i].type) {
		case "doc-start":
		case "explicit-key-ind":
		case "map-value-ind":
		case "seq-item-ind":
		case "newline": break loop;
	}
	while (prev[++i]?.type === "space");
	return prev.splice(i, prev.length);
}
function arrayPushArray(target, source) {
	if (source.length < 1e5) Array.prototype.push.apply(target, source);
	else for (let i = 0; i < source.length; ++i) target.push(source[i]);
}
function fixFlowSeqItems(fc) {
	if (fc.start.type === "flow-seq-start") {
		for (const it of fc.items) if (it.sep && !it.value && !includesToken(it.start, "explicit-key-ind") && !includesToken(it.sep, "map-value-ind")) {
			if (it.key) it.value = it.key;
			delete it.key;
			if (isFlowToken(it.value)) if (it.value.end) arrayPushArray(it.value.end, it.sep);
			else it.value.end = it.sep;
			else arrayPushArray(it.start, it.sep);
			delete it.sep;
		}
	}
}
/**
* A YAML concrete syntax tree (CST) parser
*
* ```ts
* const src: string = ...
* for (const token of new Parser().parse(src)) {
*   // token: Token
* }
* ```
*
* To use the parser with a user-provided lexer:
*
* ```ts
* function* parse(source: string, lexer: Lexer) {
*   const parser = new Parser()
*   for (const lexeme of lexer.lex(source))
*     yield* parser.next(lexeme)
*   yield* parser.end()
* }
*
* const src: string = ...
* const lexer = new Lexer()
* for (const token of parse(src, lexer)) {
*   // token: Token
* }
* ```
*/
var Parser = class {
	/**
	* @param onNewLine - If defined, called separately with the start position of
	*   each new line (in `parse()`, including the start of input).
	*/
	constructor(onNewLine) {
		/** If true, space and sequence indicators count as indentation */
		this.atNewLine = true;
		/** If true, next token is a scalar value */
		this.atScalar = false;
		/** Current indentation level */
		this.indent = 0;
		/** Current offset since the start of parsing */
		this.offset = 0;
		/** On the same line with a block map key */
		this.onKeyLine = false;
		/** Top indicates the node that's currently being built */
		this.stack = [];
		/** The source of the current token, set in parse() */
		this.source = "";
		/** The type of the current token, set in parse() */
		this.type = "";
		this.lexer = new Lexer();
		this.onNewLine = onNewLine;
	}
	/**
	* Parse `source` as a YAML stream.
	* If `incomplete`, a part of the last line may be left as a buffer for the next call.
	*
	* Errors are not thrown, but yielded as `{ type: 'error', message }` tokens.
	*
	* @returns A generator of tokens representing each directive, document, and other structure.
	*/
	*parse(source, incomplete = false) {
		if (this.onNewLine && this.offset === 0) this.onNewLine(0);
		for (const lexeme of this.lexer.lex(source, incomplete)) yield* this.next(lexeme);
		if (!incomplete) yield* this.end();
	}
	/**
	* Advance the parser by the `source` of one lexical token.
	*/
	*next(source) {
		this.source = source;
		if (this.atScalar) {
			this.atScalar = false;
			yield* this.step();
			this.offset += source.length;
			return;
		}
		const type = tokenType(source);
		if (!type) {
			const message = `Not a YAML token: ${source}`;
			yield* this.pop({
				type: "error",
				offset: this.offset,
				message,
				source
			});
			this.offset += source.length;
		} else if (type === "scalar") {
			this.atNewLine = false;
			this.atScalar = true;
			this.type = "scalar";
		} else {
			this.type = type;
			yield* this.step();
			switch (type) {
				case "newline":
					this.atNewLine = true;
					this.indent = 0;
					if (this.onNewLine) this.onNewLine(this.offset + source.length);
					break;
				case "space":
					if (this.atNewLine && source[0] === " ") this.indent += source.length;
					break;
				case "explicit-key-ind":
				case "map-value-ind":
				case "seq-item-ind":
					if (this.atNewLine) this.indent += source.length;
					break;
				case "doc-mode":
				case "flow-error-end": return;
				default: this.atNewLine = false;
			}
			this.offset += source.length;
		}
	}
	/** Call at end of input to push out any remaining constructions */
	*end() {
		while (this.stack.length > 0) yield* this.pop();
	}
	get sourceToken() {
		return {
			type: this.type,
			offset: this.offset,
			indent: this.indent,
			source: this.source
		};
	}
	*step() {
		const top = this.peek(1);
		if (this.type === "doc-end" && top?.type !== "doc-end") {
			while (this.stack.length > 0) yield* this.pop();
			this.stack.push({
				type: "doc-end",
				offset: this.offset,
				source: this.source
			});
			return;
		}
		if (!top) return yield* this.stream();
		switch (top.type) {
			case "document": return yield* this.document(top);
			case "alias":
			case "scalar":
			case "single-quoted-scalar":
			case "double-quoted-scalar": return yield* this.scalar(top);
			case "block-scalar": return yield* this.blockScalar(top);
			case "block-map": return yield* this.blockMap(top);
			case "block-seq": return yield* this.blockSequence(top);
			case "flow-collection": return yield* this.flowCollection(top);
			case "doc-end": return yield* this.documentEnd(top);
		}
		/* istanbul ignore next should not happen */
		yield* this.pop();
	}
	peek(n) {
		return this.stack[this.stack.length - n];
	}
	*pop(error) {
		const token = error ?? this.stack.pop();
		/* istanbul ignore if should not happen */
		if (!token) yield {
			type: "error",
			offset: this.offset,
			source: "",
			message: "Tried to pop an empty stack"
		};
		else if (this.stack.length === 0) yield token;
		else {
			const top = this.peek(1);
			if (token.type === "block-scalar") token.indent = "indent" in top ? top.indent : 0;
			else if (token.type === "flow-collection" && top.type === "document") token.indent = 0;
			if (token.type === "flow-collection") fixFlowSeqItems(token);
			switch (top.type) {
				case "document":
					top.value = token;
					break;
				case "block-scalar":
					top.props.push(token);
					break;
				case "block-map": {
					const it = top.items[top.items.length - 1];
					if (it.value) {
						top.items.push({
							start: [],
							key: token,
							sep: []
						});
						this.onKeyLine = true;
						return;
					} else if (it.sep) it.value = token;
					else {
						Object.assign(it, {
							key: token,
							sep: []
						});
						this.onKeyLine = !it.explicitKey;
						return;
					}
					break;
				}
				case "block-seq": {
					const it = top.items[top.items.length - 1];
					if (it.value) top.items.push({
						start: [],
						value: token
					});
					else it.value = token;
					break;
				}
				case "flow-collection": {
					const it = top.items[top.items.length - 1];
					if (!it || it.value) top.items.push({
						start: [],
						key: token,
						sep: []
					});
					else if (it.sep) it.value = token;
					else Object.assign(it, {
						key: token,
						sep: []
					});
					return;
				}
				/* istanbul ignore next should not happen */
				default:
					yield* this.pop();
					yield* this.pop(token);
			}
			if ((top.type === "document" || top.type === "block-map" || top.type === "block-seq") && (token.type === "block-map" || token.type === "block-seq")) {
				const last = token.items[token.items.length - 1];
				if (last && !last.sep && !last.value && last.start.length > 0 && findNonEmptyIndex(last.start) === -1 && (token.indent === 0 || last.start.every((st) => st.type !== "comment" || st.indent < token.indent))) {
					if (top.type === "document") top.end = last.start;
					else top.items.push({ start: last.start });
					token.items.splice(-1, 1);
				}
			}
		}
	}
	*stream() {
		switch (this.type) {
			case "directive-line":
				yield {
					type: "directive",
					offset: this.offset,
					source: this.source
				};
				return;
			case "byte-order-mark":
			case "space":
			case "comment":
			case "newline":
				yield this.sourceToken;
				return;
			case "doc-mode":
			case "doc-start": {
				const doc = {
					type: "document",
					offset: this.offset,
					start: []
				};
				if (this.type === "doc-start") doc.start.push(this.sourceToken);
				this.stack.push(doc);
				return;
			}
		}
		yield {
			type: "error",
			offset: this.offset,
			message: `Unexpected ${this.type} token in YAML stream`,
			source: this.source
		};
	}
	*document(doc) {
		if (doc.value) return yield* this.lineEnd(doc);
		switch (this.type) {
			case "doc-start":
				if (findNonEmptyIndex(doc.start) !== -1) {
					yield* this.pop();
					yield* this.step();
				} else doc.start.push(this.sourceToken);
				return;
			case "anchor":
			case "tag":
			case "space":
			case "comment":
			case "newline":
				doc.start.push(this.sourceToken);
				return;
		}
		const bv = this.startBlockValue(doc);
		if (bv) this.stack.push(bv);
		else yield {
			type: "error",
			offset: this.offset,
			message: `Unexpected ${this.type} token in YAML document`,
			source: this.source
		};
	}
	*scalar(scalar) {
		if (this.type === "map-value-ind") {
			const start = getFirstKeyStartProps(getPrevProps(this.peek(2)));
			let sep;
			if (scalar.end) {
				sep = scalar.end;
				sep.push(this.sourceToken);
				delete scalar.end;
			} else sep = [this.sourceToken];
			const map = {
				type: "block-map",
				offset: scalar.offset,
				indent: scalar.indent,
				items: [{
					start,
					key: scalar,
					sep
				}]
			};
			this.onKeyLine = true;
			this.stack[this.stack.length - 1] = map;
		} else yield* this.lineEnd(scalar);
	}
	*blockScalar(scalar) {
		switch (this.type) {
			case "space":
			case "comment":
			case "newline":
				scalar.props.push(this.sourceToken);
				return;
			case "scalar":
				scalar.source = this.source;
				this.atNewLine = true;
				this.indent = 0;
				if (this.onNewLine) {
					let nl = this.source.indexOf("\n") + 1;
					while (nl !== 0) {
						this.onNewLine(this.offset + nl);
						nl = this.source.indexOf("\n", nl) + 1;
					}
				}
				yield* this.pop();
				break;
			/* istanbul ignore next should not happen */
			default:
				yield* this.pop();
				yield* this.step();
		}
	}
	*blockMap(map) {
		const it = map.items[map.items.length - 1];
		switch (this.type) {
			case "newline":
				this.onKeyLine = false;
				if (it.value) {
					const end = "end" in it.value ? it.value.end : void 0;
					if ((Array.isArray(end) ? end[end.length - 1] : void 0)?.type === "comment") end?.push(this.sourceToken);
					else map.items.push({ start: [this.sourceToken] });
				} else if (it.sep) it.sep.push(this.sourceToken);
				else it.start.push(this.sourceToken);
				return;
			case "space":
			case "comment":
				if (it.value) map.items.push({ start: [this.sourceToken] });
				else if (it.sep) it.sep.push(this.sourceToken);
				else {
					if (this.atIndentedComment(it.start, map.indent)) {
						const end = map.items[map.items.length - 2]?.value?.end;
						if (Array.isArray(end)) {
							arrayPushArray(end, it.start);
							end.push(this.sourceToken);
							map.items.pop();
							return;
						}
					}
					it.start.push(this.sourceToken);
				}
				return;
		}
		if (this.indent >= map.indent) {
			const atMapIndent = !this.onKeyLine && this.indent === map.indent;
			const atNextItem = atMapIndent && (it.sep || it.explicitKey) && this.type !== "seq-item-ind";
			let start = [];
			if (atNextItem && it.sep && !it.value) {
				const nl = [];
				for (let i = 0; i < it.sep.length; ++i) {
					const st = it.sep[i];
					switch (st.type) {
						case "newline":
							nl.push(i);
							break;
						case "space": break;
						case "comment":
							if (st.indent > map.indent) nl.length = 0;
							break;
						default: nl.length = 0;
					}
				}
				if (nl.length >= 2) start = it.sep.splice(nl[1]);
			}
			switch (this.type) {
				case "anchor":
				case "tag":
					if (atNextItem || it.value) {
						start.push(this.sourceToken);
						map.items.push({ start });
						this.onKeyLine = true;
					} else if (it.sep) it.sep.push(this.sourceToken);
					else it.start.push(this.sourceToken);
					return;
				case "explicit-key-ind":
					if (!it.sep && !it.explicitKey) {
						it.start.push(this.sourceToken);
						it.explicitKey = true;
					} else if (atNextItem || it.value) {
						start.push(this.sourceToken);
						map.items.push({
							start,
							explicitKey: true
						});
					} else this.stack.push({
						type: "block-map",
						offset: this.offset,
						indent: this.indent,
						items: [{
							start: [this.sourceToken],
							explicitKey: true
						}]
					});
					this.onKeyLine = true;
					return;
				case "map-value-ind":
					if (it.explicitKey) if (!it.sep) if (includesToken(it.start, "newline")) Object.assign(it, {
						key: null,
						sep: [this.sourceToken]
					});
					else {
						const start = getFirstKeyStartProps(it.start);
						this.stack.push({
							type: "block-map",
							offset: this.offset,
							indent: this.indent,
							items: [{
								start,
								key: null,
								sep: [this.sourceToken]
							}]
						});
					}
					else if (it.value) map.items.push({
						start: [],
						key: null,
						sep: [this.sourceToken]
					});
					else if (includesToken(it.sep, "map-value-ind")) this.stack.push({
						type: "block-map",
						offset: this.offset,
						indent: this.indent,
						items: [{
							start,
							key: null,
							sep: [this.sourceToken]
						}]
					});
					else if (isFlowToken(it.key) && !includesToken(it.sep, "newline")) {
						const start = getFirstKeyStartProps(it.start);
						const key = it.key;
						const sep = it.sep;
						sep.push(this.sourceToken);
						delete it.key;
						delete it.sep;
						this.stack.push({
							type: "block-map",
							offset: this.offset,
							indent: this.indent,
							items: [{
								start,
								key,
								sep
							}]
						});
					} else if (start.length > 0) it.sep = it.sep.concat(start, this.sourceToken);
					else it.sep.push(this.sourceToken);
					else if (!it.sep) Object.assign(it, {
						key: null,
						sep: [this.sourceToken]
					});
					else if (it.value || atNextItem) map.items.push({
						start,
						key: null,
						sep: [this.sourceToken]
					});
					else if (includesToken(it.sep, "map-value-ind")) this.stack.push({
						type: "block-map",
						offset: this.offset,
						indent: this.indent,
						items: [{
							start: [],
							key: null,
							sep: [this.sourceToken]
						}]
					});
					else it.sep.push(this.sourceToken);
					this.onKeyLine = true;
					return;
				case "alias":
				case "scalar":
				case "single-quoted-scalar":
				case "double-quoted-scalar": {
					const fs = this.flowScalar(this.type);
					if (atNextItem || it.value) {
						map.items.push({
							start,
							key: fs,
							sep: []
						});
						this.onKeyLine = true;
					} else if (it.sep) this.stack.push(fs);
					else {
						Object.assign(it, {
							key: fs,
							sep: []
						});
						this.onKeyLine = true;
					}
					return;
				}
				default: {
					const bv = this.startBlockValue(map);
					if (bv) {
						if (bv.type === "block-seq") {
							if (!it.explicitKey && it.sep && !includesToken(it.sep, "newline")) {
								yield* this.pop({
									type: "error",
									offset: this.offset,
									message: "Unexpected block-seq-ind on same line with key",
									source: this.source
								});
								return;
							}
						} else if (atMapIndent) map.items.push({ start });
						this.stack.push(bv);
						return;
					}
				}
			}
		}
		yield* this.pop();
		yield* this.step();
	}
	*blockSequence(seq) {
		const it = seq.items[seq.items.length - 1];
		switch (this.type) {
			case "newline":
				if (it.value) {
					const end = "end" in it.value ? it.value.end : void 0;
					if ((Array.isArray(end) ? end[end.length - 1] : void 0)?.type === "comment") end?.push(this.sourceToken);
					else seq.items.push({ start: [this.sourceToken] });
				} else it.start.push(this.sourceToken);
				return;
			case "space":
			case "comment":
				if (it.value) seq.items.push({ start: [this.sourceToken] });
				else {
					if (this.atIndentedComment(it.start, seq.indent)) {
						const end = seq.items[seq.items.length - 2]?.value?.end;
						if (Array.isArray(end)) {
							arrayPushArray(end, it.start);
							end.push(this.sourceToken);
							seq.items.pop();
							return;
						}
					}
					it.start.push(this.sourceToken);
				}
				return;
			case "anchor":
			case "tag":
				if (it.value || this.indent <= seq.indent) break;
				it.start.push(this.sourceToken);
				return;
			case "seq-item-ind":
				if (this.indent !== seq.indent) break;
				if (it.value || includesToken(it.start, "seq-item-ind")) seq.items.push({ start: [this.sourceToken] });
				else it.start.push(this.sourceToken);
				return;
		}
		if (this.indent > seq.indent) {
			const bv = this.startBlockValue(seq);
			if (bv) {
				this.stack.push(bv);
				return;
			}
		}
		yield* this.pop();
		yield* this.step();
	}
	*flowCollection(fc) {
		const it = fc.items[fc.items.length - 1];
		if (this.type === "flow-error-end") {
			let top;
			do {
				yield* this.pop();
				top = this.peek(1);
			} while (top?.type === "flow-collection");
		} else if (fc.end.length === 0) {
			switch (this.type) {
				case "comma":
				case "explicit-key-ind":
					if (!it || it.sep) fc.items.push({ start: [this.sourceToken] });
					else it.start.push(this.sourceToken);
					return;
				case "map-value-ind":
					if (!it || it.value) fc.items.push({
						start: [],
						key: null,
						sep: [this.sourceToken]
					});
					else if (it.sep) it.sep.push(this.sourceToken);
					else Object.assign(it, {
						key: null,
						sep: [this.sourceToken]
					});
					return;
				case "space":
				case "comment":
				case "newline":
				case "anchor":
				case "tag":
					if (!it || it.value) fc.items.push({ start: [this.sourceToken] });
					else if (it.sep) it.sep.push(this.sourceToken);
					else it.start.push(this.sourceToken);
					return;
				case "alias":
				case "scalar":
				case "single-quoted-scalar":
				case "double-quoted-scalar": {
					const fs = this.flowScalar(this.type);
					if (!it || it.value) fc.items.push({
						start: [],
						key: fs,
						sep: []
					});
					else if (it.sep) this.stack.push(fs);
					else Object.assign(it, {
						key: fs,
						sep: []
					});
					return;
				}
				case "flow-map-end":
				case "flow-seq-end":
					fc.end.push(this.sourceToken);
					return;
			}
			const bv = this.startBlockValue(fc);
			/* istanbul ignore else should not happen */
			if (bv) this.stack.push(bv);
			else {
				yield* this.pop();
				yield* this.step();
			}
		} else {
			const parent = this.peek(2);
			if (parent.type === "block-map" && (this.type === "map-value-ind" && parent.indent === fc.indent || this.type === "newline" && !parent.items[parent.items.length - 1].sep)) {
				yield* this.pop();
				yield* this.step();
			} else if (this.type === "map-value-ind" && parent.type !== "flow-collection") {
				const start = getFirstKeyStartProps(getPrevProps(parent));
				fixFlowSeqItems(fc);
				const sep = fc.end.splice(1, fc.end.length);
				sep.push(this.sourceToken);
				const map = {
					type: "block-map",
					offset: fc.offset,
					indent: fc.indent,
					items: [{
						start,
						key: fc,
						sep
					}]
				};
				this.onKeyLine = true;
				this.stack[this.stack.length - 1] = map;
			} else yield* this.lineEnd(fc);
		}
	}
	flowScalar(type) {
		if (this.onNewLine) {
			let nl = this.source.indexOf("\n") + 1;
			while (nl !== 0) {
				this.onNewLine(this.offset + nl);
				nl = this.source.indexOf("\n", nl) + 1;
			}
		}
		return {
			type,
			offset: this.offset,
			indent: this.indent,
			source: this.source
		};
	}
	startBlockValue(parent) {
		switch (this.type) {
			case "alias":
			case "scalar":
			case "single-quoted-scalar":
			case "double-quoted-scalar": return this.flowScalar(this.type);
			case "block-scalar-header": return {
				type: "block-scalar",
				offset: this.offset,
				indent: this.indent,
				props: [this.sourceToken],
				source: ""
			};
			case "flow-map-start":
			case "flow-seq-start": return {
				type: "flow-collection",
				offset: this.offset,
				indent: this.indent,
				start: this.sourceToken,
				items: [],
				end: []
			};
			case "seq-item-ind": return {
				type: "block-seq",
				offset: this.offset,
				indent: this.indent,
				items: [{ start: [this.sourceToken] }]
			};
			case "explicit-key-ind": {
				this.onKeyLine = true;
				const start = getFirstKeyStartProps(getPrevProps(parent));
				start.push(this.sourceToken);
				return {
					type: "block-map",
					offset: this.offset,
					indent: this.indent,
					items: [{
						start,
						explicitKey: true
					}]
				};
			}
			case "map-value-ind": {
				this.onKeyLine = true;
				const start = getFirstKeyStartProps(getPrevProps(parent));
				return {
					type: "block-map",
					offset: this.offset,
					indent: this.indent,
					items: [{
						start,
						key: null,
						sep: [this.sourceToken]
					}]
				};
			}
		}
		return null;
	}
	atIndentedComment(start, indent) {
		if (this.type !== "comment") return false;
		if (this.indent <= indent) return false;
		return start.every((st) => st.type === "newline" || st.type === "space");
	}
	*documentEnd(docEnd) {
		if (this.type !== "doc-mode") {
			if (docEnd.end) docEnd.end.push(this.sourceToken);
			else docEnd.end = [this.sourceToken];
			if (this.type === "newline") yield* this.pop();
		}
	}
	*lineEnd(token) {
		switch (this.type) {
			case "comma":
			case "doc-start":
			case "doc-end":
			case "flow-seq-end":
			case "flow-map-end":
			case "map-value-ind":
				yield* this.pop();
				yield* this.step();
				break;
			case "newline": this.onKeyLine = false;
			default:
				if (token.end) token.end.push(this.sourceToken);
				else token.end = [this.sourceToken];
				if (this.type === "newline") yield* this.pop();
		}
	}
};
//#endregion
//#region node_modules/yaml/browser/dist/public-api.js
function parseOptions(options) {
	const prettyErrors = options.prettyErrors !== false;
	return {
		lineCounter: options.lineCounter || prettyErrors && new LineCounter() || null,
		prettyErrors
	};
}
/**
* Parse the input as a stream of YAML documents.
*
* Documents should be separated from each other by `...` or `---` marker lines.
*
* @returns If an empty `docs` array is returned, it will be of type
*   EmptyStream and contain additional stream information. In
*   TypeScript, you should use `'empty' in docs` as a type guard for it.
*/
function parseAllDocuments(source, options = {}) {
	const { lineCounter, prettyErrors } = parseOptions(options);
	const parser = new Parser(lineCounter?.addNewLine);
	const composer = new Composer(options);
	const docs = Array.from(composer.compose(parser.parse(source)));
	if (prettyErrors && lineCounter) for (const doc of docs) {
		doc.errors.forEach(prettifyError(source, lineCounter));
		doc.warnings.forEach(prettifyError(source, lineCounter));
	}
	if (docs.length > 0) return docs;
	return Object.assign([], { empty: true }, composer.streamInfo());
}
/** Parse an input string into a single YAML.Document */
function parseDocument(source, options = {}) {
	const { lineCounter, prettyErrors } = parseOptions(options);
	const parser = new Parser(lineCounter?.addNewLine);
	const composer = new Composer(options);
	let doc = null;
	for (const _doc of composer.compose(parser.parse(source), true, source.length)) if (!doc) doc = _doc;
	else if (doc.options.logLevel !== "silent") {
		doc.errors.push(new YAMLParseError(_doc.range.slice(0, 2), "MULTIPLE_DOCS", "Source contains multiple documents; please use YAML.parseAllDocuments()"));
		break;
	}
	if (prettyErrors && lineCounter) {
		doc.errors.forEach(prettifyError(source, lineCounter));
		doc.warnings.forEach(prettifyError(source, lineCounter));
	}
	return doc;
}
function parse(src, reviver, options) {
	let _reviver = void 0;
	if (typeof reviver === "function") _reviver = reviver;
	else if (options === void 0 && reviver && typeof reviver === "object") options = reviver;
	const doc = parseDocument(src, options);
	if (!doc) return null;
	doc.warnings.forEach((warning) => warn(doc.options.logLevel, warning));
	if (doc.errors.length > 0) if (doc.options.logLevel !== "silent") throw doc.errors[0];
	else doc.errors = [];
	return doc.toJS(Object.assign({ reviver: _reviver }, options));
}
function stringify(value, replacer, options) {
	let _replacer = null;
	if (typeof replacer === "function" || Array.isArray(replacer)) _replacer = replacer;
	else if (options === void 0 && replacer) options = replacer;
	if (typeof options === "string") options = options.length;
	if (typeof options === "number") {
		const indent = Math.round(options);
		options = indent < 1 ? void 0 : indent > 8 ? { indent: 8 } : { indent };
	}
	if (value === void 0) {
		const { keepUndefined } = options ?? replacer ?? {};
		if (!keepUndefined) return void 0;
	}
	if (isDocument(value) && !_replacer) return value.toString(options);
	return new Document(value, _replacer, options).toString(options);
}
//#endregion
//#region node_modules/yaml/browser/index.js
var browser_default = /* @__PURE__ */ __exportAll({
	Alias: () => Alias,
	CST: () => cst_exports,
	Composer: () => Composer,
	Document: () => Document,
	Lexer: () => Lexer,
	LineCounter: () => LineCounter,
	Pair: () => Pair,
	Parser: () => Parser,
	Scalar: () => Scalar,
	Schema: () => Schema,
	YAMLError: () => YAMLError,
	YAMLMap: () => YAMLMap,
	YAMLParseError: () => YAMLParseError,
	YAMLSeq: () => YAMLSeq,
	YAMLWarning: () => YAMLWarning,
	isAlias: () => isAlias,
	isCollection: () => isCollection$1,
	isDocument: () => isDocument,
	isMap: () => isMap,
	isNode: () => isNode,
	isPair: () => isPair,
	isScalar: () => isScalar$1,
	isSeq: () => isSeq,
	parse: () => parse,
	parseAllDocuments: () => parseAllDocuments,
	parseDocument: () => parseDocument,
	stringify: () => stringify,
	visit: () => visit$1,
	visitAsync: () => visitAsync
});
//#endregion
//#region src/config.ts
const DEFAULT_PARAMETERS = {
	errorRate: .05,
	minLength: 2100,
	maxLength: 4300,
	primerTolerance: 1,
	primerWindow: 200,
	primerChop: 0,
	maxReadsPerSample: 1e5,
	familySizeThreshold: 1,
	ldaThreshold: .995,
	contaminationClusterThreshold: .015,
	contaminationProportionThreshold: .2,
	contaminationDistanceThreshold: .015,
	contaminationFilter: true,
	agreementThreshold: .6,
	artefactFraction: .25,
	outlierQuantile: .99,
	panelThreshold: 50,
	functionalMatchThreshold: .7,
	spoolPartitions: 64,
	deterministicSeed: 88301613631812n
};
const mapping = (value, context) => {
	if (!value || typeof value !== "object" || Array.isArray(value)) throw new Error(`${context} must be a YAML mapping.`);
	return value;
};
const text$1 = (value, context, optional = false) => {
	if (optional && value == null) return "";
	if (typeof value !== "string" || !optional && !value.trim()) throw new Error(`${context} must be text.`);
	return value;
};
const optionalNumber$1 = (value, context) => {
	if (value == null) return void 0;
	const parsed = Number(value);
	if (!Number.isFinite(parsed)) throw new Error(`${context} must be numeric.`);
	return parsed;
};
const pick = (object, ...keys) => keys.find((key) => object[key] !== void 0) ? object[keys.find((key) => object[key] !== void 0)] : void 0;
function parseParameters(source) {
	const result = { ...DEFAULT_PARAMETERS };
	if (!source) return result;
	const assignNumber = (property, ...keys) => {
		const value = pick(source, ...keys);
		if (value !== void 0) result[property] = Number(value);
	};
	assignNumber("errorRate", "errorRate", "error_rate");
	assignNumber("minLength", "minLength", "min_length");
	assignNumber("maxLength", "maxLength", "max_length");
	assignNumber("primerTolerance", "primerTolerance", "primer_tolerance", "primer_tol");
	assignNumber("primerWindow", "primerWindow", "primer_window");
	assignNumber("primerChop", "primerChop", "primer_chop");
	assignNumber("maxReadsPerSample", "maxReadsPerSample", "max_reads_per_sample", "max_reads");
	assignNumber("familySizeThreshold", "familySizeThreshold", "family_size_threshold", "fs_thresh");
	assignNumber("ldaThreshold", "ldaThreshold", "lda_threshold", "lda_thresh");
	assignNumber("contaminationClusterThreshold", "contaminationClusterThreshold", "contamination_cluster_threshold", "cluster_thresh");
	assignNumber("contaminationProportionThreshold", "contaminationProportionThreshold", "contamination_proportion_threshold", "proportion_thresh");
	assignNumber("contaminationDistanceThreshold", "contaminationDistanceThreshold", "contamination_distance_threshold", "dist_thresh");
	assignNumber("agreementThreshold", "agreementThreshold", "agreement_threshold", "agreement_thresh");
	assignNumber("artefactFraction", "artefactFraction", "artefact_fraction", "af_thresh");
	assignNumber("outlierQuantile", "outlierQuantile", "outlier_quantile", "q_thresh");
	assignNumber("panelThreshold", "panelThreshold", "panel_threshold", "panel_thresh");
	assignNumber("functionalMatchThreshold", "functionalMatchThreshold", "functional_match_threshold", "ff_match");
	assignNumber("spoolPartitions", "spoolPartitions", "spool_partitions");
	const contamination = pick(source, "contaminationFilter", "contamination_filter", "contam_toggle");
	if (contamination !== void 0) result.contaminationFilter = contamination === true || contamination === "on" || contamination === "true";
	const seed = pick(source, "deterministicSeed", "deterministic_seed");
	if (seed !== void 0) result.deterministicSeed = BigInt(String(seed));
	for (const [key, value] of Object.entries(result)) {
		if (key === "deterministicSeed" || key === "contaminationFilter") continue;
		if (!Number.isFinite(value)) throw new Error(`Parameter ${key} must be numeric.`);
	}
	if (result.deterministicSeed < 0n || result.deterministicSeed > 18446744073709551615n) throw new Error("deterministicSeed must be an unsigned 64-bit integer.");
	for (const key of [
		"ldaThreshold",
		"contaminationProportionThreshold",
		"agreementThreshold",
		"artefactFraction",
		"outlierQuantile",
		"functionalMatchThreshold"
	]) if (!(Number(result[key]) >= 0 && Number(result[key]) <= 1)) throw new Error(`${key} must be between 0 and 1.`);
	for (const key of [
		"minLength",
		"maxLength",
		"primerTolerance",
		"primerWindow",
		"primerChop",
		"maxReadsPerSample",
		"familySizeThreshold",
		"spoolPartitions"
	]) if (!Number.isSafeInteger(Number(result[key])) || Number(result[key]) < 0) throw new Error(`${key} must be a non-negative integer.`);
	if (!(result.errorRate >= 0 && result.errorRate < 1) || result.minLength >= result.maxLength) throw new Error("Invalid error-rate or read-length bounds.");
	if (result.familySizeThreshold < 1) throw new Error("familySizeThreshold must be at least 1.");
	if (result.contaminationClusterThreshold < 0 || result.contaminationDistanceThreshold < 0 || result.panelThreshold < 0) throw new Error("Distance and panel thresholds cannot be negative.");
	if (result.primerWindow < 16) throw new Error("primerWindow must be at least 16.");
	if (result.spoolPartitions < 1 || result.spoolPartitions > 256 || (result.spoolPartitions & result.spoolPartitions - 1) !== 0) throw new Error("spoolPartitions must be a power of two from 1 to 256.");
	return result;
}
function parseSamples(source) {
	return Object.entries(source).map(([name, raw]) => {
		const value = mapping(raw, `Sample ${name}`);
		return {
			name,
			cdnaPrimer: text$1(pick(value, "cDNA_primer", "cdna_primer", "cdnaPrimer"), `${name}.cDNA_primer`),
			secondStrandPrimer: text$1(pick(value, "sec_str_primer", "secondStrandPrimer"), `${name}.sec_str_primer`),
			panel: text$1(value.panel, `${name}.panel`),
			functionalReference: text$1(pick(value, "ff_ref", "functionalReference"), `${name}.ff_ref`, true) || void 0,
			panelSequences: [],
			familySizeOverride: optionalNumber$1(pick(value, "fs_override", "familySizeOverride"), `${name}.fs_override`),
			artefactFractionOverride: optionalNumber$1(pick(value, "af_override", "artefactFractionOverride"), `${name}.af_override`),
			outlierQuantileOverride: optionalNumber$1(pick(value, "q_override", "outlierQuantileOverride"), `${name}.q_override`),
			agreementOverride: optionalNumber$1(pick(value, "ma_override", "agreementOverride"), `${name}.ma_override`),
			functionalMatchOverride: optionalNumber$1(pick(value, "ff_match_override", "functionalMatchOverride"), `${name}.ff_match_override`)
		};
	});
}
function parseConfigYaml(source) {
	const parsed = mapping(browser_default.parse(source), "Configuration");
	let dataset, samples, parameters;
	if (typeof parsed.dataset === "string" && parsed.samples) {
		dataset = parsed.dataset;
		samples = mapping(parsed.samples, "samples");
		parameters = parsed.parameters ? mapping(parsed.parameters, "parameters") : void 0;
	} else {
		const entries = Object.entries(parsed);
		if (entries.length !== 1) throw new Error("Original PORPID YAML must contain exactly one dataset.");
		dataset = entries[0][0];
		samples = mapping(entries[0][1], dataset);
	}
	const contaminationPanel = typeof parsed.contaminationPanel === "string" ? parsed.contaminationPanel : typeof parsed.contamination_panel === "string" ? parsed.contamination_panel : "panels/contam_panel.fasta";
	const config = {
		dataset,
		samples: parseSamples(samples),
		contaminationPanel,
		contaminationPanelSequences: [],
		parameters: parseParameters(parameters)
	};
	if (!config.dataset.trim() || !config.samples.length) throw new Error("The configuration needs a dataset name and at least one sample.");
	const assays = /* @__PURE__ */ new Set();
	for (const sample of config.samples) {
		if (!sample.name.trim()) throw new Error("Every sample needs a non-empty name.");
		if (!sample.cdnaPrimer.match(/[a-z]+/)?.[0].toUpperCase()) throw new Error(`${sample.name} cDNA_primer has no lower-case sample ID.`);
		const assay = `${sample.secondStrandPrimer.toUpperCase()}\0${sample.cdnaPrimer.toUpperCase()}`;
		if (assays.has(assay)) throw new Error(`${sample.name} duplicates another sample's complete primer pair.`);
		assays.add(assay);
		if (!/[Nn]+/.test(sample.cdnaPrimer)) throw new Error(`${sample.name} cDNA_primer has no N-marked UMI.`);
		if (sample.familySizeOverride !== void 0 && (!Number.isSafeInteger(sample.familySizeOverride) || sample.familySizeOverride < 1)) throw new Error(`${sample.name}.fs_override must be an integer of at least 1.`);
		for (const [label, value] of [
			["af_override", sample.artefactFractionOverride],
			["q_override", sample.outlierQuantileOverride],
			["ma_override", sample.agreementOverride],
			["ff_match_override", sample.functionalMatchOverride]
		]) if (value !== void 0 && !(value >= 0 && value <= 1)) throw new Error(`${sample.name}.${label} must be between 0 and 1.`);
	}
	return config;
}
function parseFasta(source) {
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
const normalizePath = (value) => value.replaceAll("\\", "/").replace(/^\.\//, "");
const basename$1 = (value) => normalizePath(value).split("/").at(-1);
async function resolveReferenceFiles(config, files) {
	const cache = /* @__PURE__ */ new Map();
	const find = (requested) => {
		const normalized = normalizePath(requested);
		const candidates = [...files.entries()].filter(([name]) => normalizePath(name) === normalized || basename$1(name) === basename$1(normalized));
		if (candidates.length !== 1) throw new Error(candidates.length ? `Reference ${requested} is ambiguous.` : `Reference ${requested} was not supplied.`);
		return candidates[0];
	};
	const load = async (requested) => {
		if (cache.has(requested)) return cache.get(requested);
		const [, content] = find(requested);
		const parsed = parseFasta(await content());
		cache.set(requested, parsed);
		return parsed;
	};
	for (const sample of config.samples) {
		sample.panelSequences = await load(sample.panel);
		if (sample.functionalReference) sample.functionalReferenceSequence = (await load(sample.functionalReference))[0];
	}
	config.contaminationPanelSequences = config.parameters.contaminationFilter ? await load(config.contaminationPanel) : [];
	return config;
}
const nan = (value) => value ?? NaN;
function compileConfig(config) {
	const writer = new BinaryWriter();
	const p = config.parameters;
	writer.magic("WPC1");
	writer.u32(1);
	writer.string(config.dataset);
	writer.f64(p.errorRate);
	writer.u32(p.minLength);
	writer.u32(p.maxLength);
	writer.u32(p.primerTolerance);
	writer.u32(p.primerWindow);
	writer.u32(p.primerChop);
	writer.u32(p.maxReadsPerSample);
	writer.u32(p.familySizeThreshold);
	writer.f64(p.ldaThreshold);
	writer.f64(p.contaminationClusterThreshold);
	writer.f64(p.contaminationProportionThreshold);
	writer.f64(p.contaminationDistanceThreshold);
	writer.u8(p.contaminationFilter ? 1 : 0);
	writer.f64(p.agreementThreshold);
	writer.f64(p.artefactFraction);
	writer.f64(p.outlierQuantile);
	writer.f64(p.panelThreshold);
	writer.f64(p.functionalMatchThreshold);
	writer.u32(p.spoolPartitions);
	writer.u64(p.deterministicSeed);
	writer.u32(config.samples.length);
	for (const sample of config.samples) {
		writer.string(sample.name);
		writer.string(sample.cdnaPrimer);
		writer.string(sample.secondStrandPrimer);
		writer.string(sample.panel);
		writer.string(sample.functionalReference ?? "");
		writer.i32(sample.familySizeOverride ?? -1);
		writer.f64(nan(sample.artefactFractionOverride));
		writer.f64(nan(sample.outlierQuantileOverride));
		writer.f64(nan(sample.agreementOverride));
		writer.f64(nan(sample.functionalMatchOverride));
		writer.u32(sample.panelSequences.length);
		for (const sequence of sample.panelSequences) {
			writer.string(sequence.name);
			writer.string(sequence.sequence);
		}
		writer.string(sample.functionalReferenceSequence?.sequence ?? "");
	}
	return writer.finish();
}
function resultConfig(config) {
	return {
		dataset: config.dataset,
		samples: config.samples.map(({ panelSequences: _panel, functionalReferenceSequence: _reference, ...sample }) => sample),
		contaminationPanel: config.contaminationPanel,
		parameters: {
			...config.parameters,
			deterministicSeed: config.parameters.deterministicSeed.toString()
		}
	};
}
//#endregion
//#region src/alivibe-msa-codec.ts
const MAGIC$1 = [
	65,
	77,
	83,
	65
];
const encoder = new TextEncoder();
const decoder = new TextDecoder("utf-8", { fatal: true });
function encodeAlivibeMsaSequences(sequences) {
	if (sequences.length > 1e4) throw new Error("Alivibe MSA accepts at most 10,000 rows.");
	const rows = sequences.map((sequence, row) => {
		if (!/^[\x21-\x7e]*$/.test(sequence)) throw new Error(`Alivibe MSA row ${row + 1} contains a non-ASCII character.`);
		return encoder.encode(sequence);
	});
	const byteLength = 8 + rows.reduce((total, row) => total + 4 + row.byteLength, 0);
	if (byteLength > 256 * 1024 * 1024) throw new Error("Alivibe MSA input is too large.");
	const buffer = new ArrayBuffer(byteLength), bytes = new Uint8Array(buffer), view = new DataView(buffer);
	bytes.set(MAGIC$1);
	view.setUint32(4, rows.length, true);
	let offset = 8;
	for (const row of rows) {
		view.setUint32(offset, row.byteLength, true);
		offset += 4;
		bytes.set(row, offset);
		offset += row.byteLength;
	}
	return buffer;
}
function decodeAlivibeMsaSequences(buffer) {
	const bytes = new Uint8Array(buffer);
	if (bytes.byteLength < 8 || MAGIC$1.some((value, index) => bytes[index] !== value)) throw new Error("Alivibe MSA returned an invalid result.");
	const view = new DataView(buffer), count = view.getUint32(4, true), rows = [];
	let offset = 8;
	if (count > 1e4) throw new Error("Alivibe MSA returned too many rows.");
	for (let row = 0; row < count; row++) {
		if (offset + 4 > bytes.byteLength) throw new Error("Alivibe MSA returned a truncated result.");
		const length = view.getUint32(offset, true);
		offset += 4;
		if (offset + length > bytes.byteLength) throw new Error("Alivibe MSA returned a truncated row.");
		rows.push(decoder.decode(bytes.subarray(offset, offset + length)));
		offset += length;
	}
	if (offset !== bytes.byteLength) throw new Error("Alivibe MSA returned trailing bytes.");
	return rows;
}
function assertAlivibeMsaResult(input, aligned) {
	if (aligned.length !== input.length) throw new Error("Alivibe MSA returned the wrong number of rows.");
	const width = aligned[0]?.length ?? 0;
	aligned.forEach((sequence) => {
		if (sequence.length !== width) throw new Error("Alivibe MSA returned a non-rectangular alignment.");
	});
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
function pairwise(reference, query) {
	const columns = query.length + 1, trace = new Uint8Array((reference.length + 1) * columns);
	let previous = new Int32Array(columns), current = new Int32Array(columns);
	for (let column = 1; column < columns; column++) {
		previous[column] = previous[column - 1] - 99;
		trace[column] = 2;
	}
	for (let row = 1; row <= reference.length; row++) {
		current[0] = previous[0] - 99;
		trace[row * columns] = 3;
		for (let column = 1; column < columns; column++) {
			const diagonal = previous[column - 1] + (reference[row - 1] === query[column - 1] ? 100 : -100);
			const left = current[column - 1] - (row === reference.length ? 99 : 100), up = previous[column] - (column === query.length ? 99 : 100);
			if (diagonal >= left && diagonal >= up) {
				current[column] = diagonal;
				trace[row * columns + column] = 1;
			} else if (left >= up) {
				current[column] = left;
				trace[row * columns + column] = 2;
			} else {
				current[column] = up;
				trace[row * columns + column] = 3;
			}
		}
		[previous, current] = [current, previous];
	}
	let left = "", right = "", row = reference.length, column = query.length;
	while (row || column) {
		const op = trace[row * columns + column];
		if (op === 1) {
			left += reference[--row];
			right += query[--column];
		} else if (op === 2) {
			left += "-";
			right += query[--column];
		} else {
			left += reference[--row];
			right += "-";
		}
	}
	return {
		reference: [...left].reverse().join(""),
		query: [...right].reverse().join("")
	};
}
function functionalFilter(reference, sequence, threshold) {
	sequence = degap(sequence);
	const reasons = [];
	if (/[^ACGT]/.test(sequence)) return {
		passed: false,
		reasons: ["ambiguousSymbols-reject"]
	};
	const coding = longestOrf(sequence);
	if (!coding || coding.length % 3) return {
		passed: false,
		reasons: ["frameshift-reject"]
	};
	const aligned = pairwise(degap(reference), coding), first = aligned.reference.search(/[^-]/);
	let last = aligned.reference.length - 1;
	while (last >= first && aligned.reference[last] === "-") last--;
	const referenceRegion = aligned.reference.slice(first, last + 1), queryRegion = aligned.query.slice(first, last + 1);
	if ([referenceRegion, queryRegion].some((row) => {
		for (const match of row.matchAll(/-+/g)) {
			const start = match.index, end = start + match[0].length;
			if (start > 0 && end < row.length && match[0].length % 3 !== 0) return true;
		}
		return false;
	})) reasons.push("frameshift-reject");
	const trimmed = degap(queryRegion), aa = translate(trimmed);
	if (queryRegion.slice(0, 3) !== "ATG") reasons.push("lateStart-reject");
	if (queryRegion.slice(-3) === "---") reasons.push("earlyStop-reject");
	let matches = 0;
	for (let index = 0; index < referenceRegion.length; index++) if (referenceRegion[index] === queryRegion[index]) matches++;
	const rawRatio = matches / Math.max(1, referenceRegion.length);
	const digits = rawRatio === 0 ? 3 : 3 - Math.floor(Math.log10(Math.abs(rawRatio))) - 1;
	const ratio = Number(rawRatio.toFixed(Math.max(0, digits)));
	if (ratio < threshold) reasons.push(`badMatch-reject (match=${ratio})`);
	return {
		passed: !reasons.length,
		reasons,
		nt: trimmed,
		aa
	};
}
async function postprocess(consensuses, contamination, config, signal, runMsa = runAlivibeMsa) {
	const discarded = new Set(contamination.filter((call) => call.discarded).map((call) => call.sequenceId));
	const records = [], summaries = [], alignments = {};
	for (const sample of config.samples) {
		if (signal?.aborted) throw new DOMException("Analysis cancelled.", "AbortError");
		const source = consensuses.filter((record) => record.sample === sample.name), sizes = source.filter((record) => !discarded.has(record.id)).map((record) => record.familySize);
		const artefactCutoff = Math.ceil(quantile(sizes, sample.outlierQuantileOverride ?? config.parameters.outlierQuantile) * (sample.artefactFractionOverride ?? config.parameters.artefactFraction));
		const agreementThreshold = sample.agreementOverride ?? config.parameters.agreementThreshold;
		const preliminary = source.map((record, index) => ({
			record,
			index
		})).filter(({ record }) => record.familySize >= artefactCutoff && record.minimumAgreement >= agreementThreshold && !discarded.has(record.id));
		const scores = Array(source.length).fill(0), panelPass = Array(source.length).fill(true), extracted = /* @__PURE__ */ new Map();
		if (preliminary.length && sample.panelSequences.length) {
			const panelResult = extractAndScorePanel(preliminary.length > 1 ? await runScalableMsa(preliminary.map(({ record }) => degap(record.sequence)), runMsa, signal, 3, "nucleotide") : preliminary.map(({ record }) => degap(record.sequence)), sample.panelSequences.map((record) => record.sequence));
			preliminary.forEach(({ index }, candidate) => {
				extracted.set(index, degap(panelResult.sequences[candidate]));
				scores[index] = panelResult.scores[candidate];
				panelPass[index] = scores[index] < config.parameters.panelThreshold;
			});
		} else preliminary.forEach(({ record, index }) => extracted.set(index, degap(record.sequence)));
		const accepted = preliminary.filter(({ index }) => panelPass[index]);
		const acceptedAlignment = accepted.length > 1 ? await runScalableMsa(accepted.map(({ index }) => extracted.get(index)), runMsa, signal, 3, "nucleotide") : accepted.map(({ index }) => extracted.get(index));
		const alignmentByIndex = new Map(accepted.map(({ index }, position) => [index, acceptedAlignment[position]]));
		const consensus = alignmentConsensus(acceptedAlignment), nucleotideRows = [], proteinRows = [];
		let functionalPassed = 0;
		source.forEach((record, index) => {
			const artefactPass = record.familySize >= artefactCutoff, agreementPass = record.minimumAgreement >= agreementThreshold;
			const contaminationPass = !discarded.has(record.id), acceptedRow = alignmentByIndex.get(index), rejectionReasons = [];
			if (!artefactPass) rejectionReasons.push(`ccs_count < artefact cutoff (${artefactCutoff})`);
			if (!agreementPass) rejectionReasons.push(`minimum_agreement < ${agreementThreshold}`);
			if (!contaminationPass) rejectionReasons.push("contamination filter");
			if (!panelPass[index]) rejectionReasons.push(`distance_from_panel >= ${config.parameters.panelThreshold}`);
			let trimmedNt, trimmedAa, functionalPass;
			if (sample.functionalReferenceSequence) if (acceptedRow) {
				const outcome = functionalFilter(sample.functionalReferenceSequence.sequence, extracted.get(index), sample.functionalMatchOverride ?? config.parameters.functionalMatchThreshold);
				functionalPass = outcome.passed;
				trimmedNt = outcome.nt;
				trimmedAa = outcome.aa;
				rejectionReasons.push(...outcome.reasons);
				if (outcome.passed) {
					functionalPassed++;
					proteinRows.push({
						name: record.id,
						sequence: outcome.aa
					});
				}
			} else functionalPass = false;
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
		});
		if (nucleotideRows.length) alignments[`${sample.name}/nucleotide`] = fasta$1(nucleotideRows);
		if (proteinRows.length) {
			const proteinAlignment = proteinRows.length > 1 ? await runScalableMsa(proteinRows.map((row) => row.sequence), runMsa, signal, 3, "amino-acid") : proteinRows.map((row) => row.sequence);
			alignments[`${sample.name}/protein`] = fasta$1(proteinRows.map((row, index) => ({
				...row,
				sequence: proteinAlignment[index]
			})));
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
	}
	return {
		records,
		summaries,
		alignments
	};
}
//#endregion
//#region node_modules/@msgpack/msgpack/dist.esm/utils/utf8.mjs
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
//#region node_modules/@msgpack/msgpack/dist.esm/ExtData.mjs
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
//#region node_modules/@msgpack/msgpack/dist.esm/DecodeError.mjs
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
//#region node_modules/@msgpack/msgpack/dist.esm/ExtensionCodec.mjs
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
//#region node_modules/@msgpack/msgpack/dist.esm/utils/typedArrays.mjs
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
//#region node_modules/@msgpack/msgpack/dist.esm/encode.mjs
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
//#region node_modules/@msgpack/msgpack/dist.esm/utils/prettyByte.mjs
function prettyByte(byte) {
	return `${byte < 0 ? "-" : ""}0x${Math.abs(byte).toString(16).padStart(2, "0")}`;
}
//#endregion
//#region node_modules/@msgpack/msgpack/dist.esm/CachedKeyDecoder.mjs
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
//#region node_modules/@msgpack/msgpack/dist.esm/Decoder.mjs
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
//#region node_modules/@msgpack/msgpack/dist.esm/decode.mjs
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
//#region node_modules/fflate/esm/browser.js
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
			"consensusSequences",
			"contaminationPassed",
			"postprocPassed",
			"artefactCutoff"
		]) count(row[key], `summaries[${index}].${key}`);
		if (row.functionalPassed != null) count(row.functionalPassed, `summaries[${index}].functionalPassed`);
	});
	if (summarySamples.size !== samples.length) throw new Error("Result summaries are missing a configured sample.");
	const familyKeys = /* @__PURE__ */ new Set();
	array(bundle.umiFamilies, "umiFamilies").forEach((entry, index) => {
		const row = object(entry, `umiFamilies[${index}]`), sample = text(row.sample, `umiFamilies[${index}].sample`), sampleIndex = count(row.sampleIndex, `umiFamilies[${index}].sampleIndex`);
		knownSample(sample, `umiFamilies[${index}]`);
		if (sampleIndices.get(sample) !== sampleIndex) throw new Error("A UMI family has an inconsistent sample index.");
		const familyKey = `${sampleIndex}\0${text(row.umi, `umiFamilies[${index}].umi`)}`;
		if (familyKeys.has(familyKey)) throw new Error("UMI family identifiers must be unique within a sample.");
		familyKeys.add(familyKey);
		count(row.familySize, `umiFamilies[${index}].familySize`);
		text(row.mostLikelyParent, `umiFamilies[${index}].mostLikelyParent`);
		numeric(row.posteriorProbability, `umiFamilies[${index}].posteriorProbability`);
		numeric(row.logOffspringProbability, `umiFamilies[${index}].logOffspringProbability`, false);
		if (!DISPOSITIONS.has(text(row.disposition, `umiFamilies[${index}].disposition`))) throw new Error("A UMI family has an unknown disposition.");
		optionalNumber(row.minimumAgreement, `umiFamilies[${index}].minimumAgreement`);
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
	const recordIds = /* @__PURE__ */ new Set();
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
		optionalText(row.alignedNt, "postproc aligned sequence");
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
	});
	if (recordIds.size !== consensusIds.size || [...consensusIds].some((id) => !recordIds.has(id))) throw new Error("Consensus and post-processing records are inconsistent.");
	for (const [label, entries] of [["alignments", object(bundle.alignments, "alignments")], ["trees", object(bundle.trees, "trees")]]) for (const [name, contents] of Object.entries(entries)) {
		text(name, `${label} name`);
		text(contents, `${label}.${name}`);
		const sample = name.split("/", 1)[0];
		knownSample(sample, `${label}.${name}`);
	}
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
			], bundle.contamination.filter((row) => !sample || row.sample === sample).map((row) => [
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
		case "nucleotide-alignment": {
			const selected = alignmentSample(bundle, sample);
			return {
				extension: "nucleotide-alignment.fasta",
				mime: "text/x-fasta",
				text: bundle.alignments[`${selected}/nucleotide`] ?? ""
			};
		}
		case "protein-alignment": {
			const selected = alignmentSample(bundle, sample);
			return {
				extension: "protein-alignment.fasta",
				mime: "text/x-fasta",
				text: bundle.alignments[`${selected}/protein`] ?? ""
			};
		}
		case "newick": {
			const selected = alignmentSample(bundle, sample);
			return {
				extension: "tree.newick",
				mime: "text/plain",
				text: bundle.trees[`${selected}/nucleotide`] ?? ""
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
//#region cli-src/direct-fasttree.mjs
function completeNewick(output) {
	const candidates = output.split(/\r?\n/).filter((line) => line.includes("(") && line.includes(";")).reverse();
	for (const candidate of candidates) {
		const start = candidate.indexOf("("), end = candidate.lastIndexOf(";");
		if (start >= 0 && end > start) return candidate.slice(start, end + 1);
	}
	throw new Error("FastTree did not return a complete Newick tree.");
}
function createFastTreeRunner(javascriptPath, wasmPath) {
	const factory = createRequire(import.meta.url)(javascriptPath);
	let output = [], errors = [], invocation = 0;
	const modulePromise = readFile(wasmPath).then((wasmBinary) => factory({
		wasmBinary,
		print: (line) => output.push(String(line)),
		printErr: (line) => errors.push(String(line))
	}));
	return async (alignedFasta) => {
		const records = parseFasta(alignedFasta);
		const safeName = (name, index) => name.replace(/[^A-Za-z0-9_.|*+\-]/g, "_") || `tip_${index + 1}`;
		if (records.length < 3) return records.length === 2 ? `(${safeName(records[0].name, 0)}:0.0,${safeName(records[1].name, 1)}:0.0);` : records.length === 1 ? `(${safeName(records[0].name, 0)}:0.0);` : ";";
		const width = records[0].sequence.length;
		if (!width || records.some((record) => record.sequence.length !== width)) throw new Error("FastTree input must be a rectangular alignment.");
		const numeric = records.map((record, index) => `>${index}\n${record.sequence}`).join("\n") + "\n";
		const runtime = await modulePromise, path = `/webporpid-fasttree-${invocation++}.fa`;
		output = [];
		errors = [];
		runtime.FS.writeFile(path, numeric);
		try {
			runtime.callMain([
				"-nosupport",
				"-nt",
				"-gtr",
				path
			]);
		} finally {
			try {
				runtime.FS.unlink(path);
			} catch {}
		}
		let tree = completeNewick(output.join("\n"));
		records.forEach((record, index) => {
			const safe = safeName(record.name, index);
			tree = tree.replace(new RegExp(`([,(])${index}(?=[:),])`, "g"), `$1${safe}`);
		});
		return tree;
	};
}
//#endregion
//#region cli-src/direct-msa.mjs
function createMsaRunner(wasmPath) {
	let runtimePromise;
	const initialize = async () => {
		const wasi = new WASI([], [], []), module = await WebAssembly.compile(await readFile(wasmPath));
		const instance = await WebAssembly.instantiate(module, { wasi_snapshot_preview1: wasi.wasiImport });
		wasi.initialize(instance);
		return instance.exports;
	};
	return async (sequences, signal, iterations = 3, scoringMode = "nucleotide") => {
		if (signal?.aborted) throw new Error("Analysis cancelled.");
		if (sequences.length < 2) return [...sequences];
		const inputSequences = sequences.map(String), input = new Uint8Array(encodeAlivibeMsaSequences(inputSequences));
		const runtime = await (runtimePromise ??= initialize()), pointer = runtime.alivibe_msa_alloc(input.byteLength);
		if (!pointer && input.byteLength) throw new Error("Alivibe MSA ran out of WebAssembly memory.");
		try {
			new Uint8Array(runtime.memory.buffer, pointer, input.byteLength).set(input);
			if ((scoringMode === "nucleotide" ? runtime.alivibe_msa_run_nucleotide : scoringMode === "amino-acid" ? runtime.alivibe_msa_run_amino_acid : runtime.alivibe_msa_run)(pointer, input.byteLength, iterations) < 0) {
				const message = new TextDecoder().decode(new Uint8Array(runtime.memory.buffer, runtime.alivibe_msa_error_ptr(), runtime.alivibe_msa_error_len()));
				throw new Error(message || "Alivibe MSA failed.");
			}
			const aligned = decodeAlivibeMsaSequences(new Uint8Array(runtime.memory.buffer, runtime.alivibe_msa_result_ptr(), runtime.alivibe_msa_result_len()).slice().buffer);
			assertAlivibeMsaResult(inputSequences, aligned);
			return aligned;
		} finally {
			runtime.alivibe_msa_free(pointer);
		}
	};
}
//#endregion
//#region cli-src/porpid-cli.mjs
const VERSION = "0.1.1";
const UPSTREAM_COMMIT = "201af7942029cfb7974880e41674be9f0ddfaf3b";
const CLI_DIRECTORY = dirname(new URL(import.meta.url).pathname);
function defaultCliAssets() {
	const directory = join(CLI_DIRECTORY, "assets");
	return {
		wasmPath: join(directory, "webporpid.wasm"),
		msaPath: join(directory, "alivibe-msa.wasm"),
		fastTreeJavascriptPath: join(directory, "fasttree.cjs"),
		fastTreeWasmPath: join(directory, "fasttree.wasm")
	};
}
function usage() {
	return `porpid-cli ${VERSION}\n\nRun the complete nanopore/PacBio pipeline:\n  porpid-cli run reads.fastq.gz --config config.yaml --output results.webporpid [--workers N]\n\nInspect or export a saved analysis:\n  porpid-cli inspect results.webporpid\n  porpid-cli export results.webporpid --component consensus-fasta [--sample NAME] --output consensus.fasta\n\nWorkers default to all logical CPUs (${availableParallelism()}). Temporary read partitions are streamed to disk and removed after consensus.\nComponents: consensus-fasta, passed-consensus-fasta, rejected-consensus-fasta, trimmed-nt-fasta, trimmed-aa-fasta,\n            family-csv, low-agreement-csv, contamination-csv, postproc-csv, apobec-csv,\n            nucleotide-alignment, protein-alignment, newick, log`;
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
	const config = parseConfigYaml(await readFile(path, "utf8")), base = dirname(resolve(path));
	const requested = new Set([config.contaminationPanel]);
	for (const sample of config.samples) {
		requested.add(sample.panel);
		if (sample.functionalReference) requested.add(sample.functionalReference);
	}
	return resolveReferenceFiles(config, new Map([...requested].map((name) => [name, () => readFile(isAbsolute(name) ? name : resolve(base, name), "utf8")])));
}
async function runPipeline({ inputPath, configPath, outputPath, workers, assets }) {
	const input = resolve(inputPath), configuration = resolve(configPath);
	await stat(input);
	const config = await loadConfiguration(configuration), compiledConfig = compileConfig(config);
	const configHash = createHash("sha256").update(compiledConfig).digest("hex"), inputHash = createHash("sha256");
	status(`starting ${config.dataset} with ${workers} workers`);
	const pool = await WorkerPool.create(workers, assets.wasmPath, compiledConfig), store = await DiskPartitions.create(config.parameters.spoolPartitions);
	const log = [
		`${now()} webPORPID ${VERSION} started`,
		`${now()} execution: ${workers} WASM workers; disk-backed partition spool`,
		`${now()} parameters: error_rate=${config.parameters.errorRate}, lengths=(${config.parameters.minLength},${config.parameters.maxLength}), lda=${config.parameters.ldaThreshold}`
	];
	try {
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
		const mergedCounts = mergeFamilyCounts(countParts), selectedReads = decodeFamilyCounts(mergedCounts).reduce((sum, entry) => sum + entry.count, 0);
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
		const contamination = classifyContamination(consensuses, config);
		log.push(`${now()} contamination: ${contamination.filter((call) => call.discarded).length} discarded; ${contamination.filter((call) => call.suspectOnly).length} suspect calls`);
		const downstream = await postprocess(consensuses, contamination, config, void 0, createMsaRunner(assets.msaPath));
		downstream.summaries.forEach((summary, index) => {
			summary.demultiplexedReads = quality.perSample[index] ?? 0;
			summary.observedUmis = umiFamilies.filter((family) => family.sampleIndex === index && family.disposition !== "BPB-rejects").length;
			summary.likelyRealUmis = umiFamilies.filter((family) => family.sampleIndex === index && family.disposition === "likely_real").length;
		});
		const trees = {}, fastTree = createFastTreeRunner(assets.fastTreeJavascriptPath, assets.fastTreeWasmPath);
		for (const [name, alignment] of Object.entries(downstream.alignments).filter(([name]) => name.endsWith("/nucleotide"))) {
			status(`FastTree: ${name.split("/")[0]}`);
			trees[name] = await fastTree(alignment);
		}
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
			records: downstream.records,
			alignments: downstream.alignments,
			trees,
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
		components: {
			consensuses: result.consensuses.length,
			families: result.umiFamilies.length,
			contaminationCalls: result.contamination.length,
			records: result.records.length,
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
		assets: overrideAssets ?? defaultCliAssets()
	});
}
//#endregion
//#region cli-src/porpid-cli-node.mjs
runCli().catch((cause) => {
	process.stderr.write(`porpid-cli: ${cause instanceof Error ? cause.message : String(cause)}\n`);
	process.exitCode = 1;
});
//#endregion
