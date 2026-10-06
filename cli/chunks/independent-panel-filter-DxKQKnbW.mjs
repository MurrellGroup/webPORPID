//#region src/independent-panel-filter.ts
const ASCII_SYMBOLS = 128;
const GAP_OPEN = -2;
const GAP_EXTEND = -.2;
const NEGATIVE = -1e30;
const LOG_QUARTER = Math.log(.25);
function normalizeQuery(sequence) {
	const result = sequence.replaceAll("-", "").replace(/\s/g, "").toUpperCase();
	if (!result.length) throw new Error("The independent panel filter received an empty sequence.");
	return result;
}
/** Build the immutable panel profile once for every query in a worker chunk. */
function buildPanelProfile(panelRows) {
	if (!panelRows.length) throw new Error("The independent panel filter requires a reference panel.");
	const normalized = panelRows.map((row) => row.replace(/\s/g, "").toUpperCase());
	const width = normalized[0].length;
	if (!width || normalized.some((row) => row.length !== width)) throw new Error("The reference panel must be a non-empty rectangular alignment.");
	const probabilities = new Float32Array(width * ASCII_SYMBOLS), weight = 1 / normalized.length;
	for (const row of normalized) for (let column = 0; column < width; column++) {
		const code = row.charCodeAt(column);
		if (code < ASCII_SYMBOLS) probabilities[column * ASCII_SYMBOLS + code] += weight;
	}
	return {
		width,
		probabilities
	};
}
/**
* Global affine profile/query alignment using O(panelLength * band) work.
*
* The band follows the length-scaled diagonal, so a net length difference is
* absorbed without widening it. If the traceback approaches a band edge the
* caller retries with a wider band; this preserves large local indels without
* paying full quadratic cost for ordinary amplicons.
*/
function bandedProfileAlignment(profile, rawQuery, band) {
	const query = normalizeQuery(rawQuery), panelLength = profile.width, queryLength = query.length;
	const stride = Math.min(queryLength + 1, 2 * band + 3), cells = (panelLength + 1) * stride;
	if (!Number.isSafeInteger(cells) || cells > 768 * 1024 * 1024) throw new Error("The independent panel alignment exceeds the browser memory safety limit.");
	const trace = new Uint8Array(cells), starts = new Int32Array(panelLength + 1), ends = new Int32Array(panelLength + 1);
	let previousM = new Float32Array(stride), previousX = new Float32Array(stride), previousY = new Float32Array(stride);
	let currentM = new Float32Array(stride), currentX = new Float32Array(stride), currentY = new Float32Array(stride);
	previousM.fill(NEGATIVE);
	previousX.fill(NEGATIVE);
	previousY.fill(NEGATIVE);
	starts[0] = 0;
	ends[0] = Math.min(queryLength, band);
	previousM[0] = 0;
	for (let column = 1; column <= ends[0]; column++) {
		previousY[column] = GAP_OPEN + GAP_EXTEND * column;
		trace[column] = column > 1 ? 8 : 0;
	}
	let previousStart = 0, previousEnd = ends[0];
	for (let row = 1; row <= panelLength; row++) {
		const center = Math.round(row * queryLength / panelLength);
		const start = Math.max(0, center - band), end = Math.min(queryLength, center + band);
		starts[row] = start;
		ends[row] = end;
		currentM.fill(NEGATIVE);
		currentX.fill(NEGATIVE);
		currentY.fill(NEGATIVE);
		if (start === 0) {
			currentX[0] = GAP_OPEN + GAP_EXTEND * row;
			trace[row * stride] = row > 1 ? 4 : 0;
		}
		for (let column = Math.max(1, start); column <= end; column++) {
			const index = column - start, packedIndex = row * stride + index;
			const diagonal = column - 1;
			let bestM = NEGATIVE, fromM = 0;
			if (diagonal >= previousStart && diagonal <= previousEnd) {
				const previousIndex = diagonal - previousStart;
				bestM = previousM[previousIndex];
				if (previousX[previousIndex] > bestM) {
					bestM = previousX[previousIndex];
					fromM = 1;
				}
				if (previousY[previousIndex] > bestM) {
					bestM = previousY[previousIndex];
					fromM = 2;
				}
				const code = query.charCodeAt(column - 1), probability = code < ASCII_SYMBOLS ? profile.probabilities[(row - 1) * ASCII_SYMBOLS + code] : 0;
				currentM[index] = bestM + 2 * probability - 1;
			}
			let xContinues = false;
			if (column >= previousStart && column <= previousEnd) {
				const previousIndex = column - previousStart;
				const opened = previousM[previousIndex] + GAP_OPEN, extended = previousX[previousIndex] + GAP_EXTEND;
				xContinues = extended > opened;
				currentX[index] = xContinues ? extended : opened;
			}
			let yContinues = false;
			if (column - 1 >= start) {
				const leftIndex = index - 1, opened = currentM[leftIndex] + GAP_OPEN, extended = currentY[leftIndex] + GAP_EXTEND;
				yContinues = extended > opened;
				currentY[index] = yContinues ? extended : opened;
			}
			trace[packedIndex] = fromM | (xContinues ? 4 : 0) | (yContinues ? 8 : 0);
		}
		[previousM, currentM] = [currentM, previousM];
		[previousX, currentX] = [currentX, previousX];
		[previousY, currentY] = [currentY, previousY];
		previousStart = start;
		previousEnd = end;
	}
	if (queryLength < previousStart || queryLength > previousEnd) throw new Error("The independent panel alignment did not reach its endpoint.");
	const endpoint = queryLength - previousStart;
	let state = 0, maximum = previousM[endpoint];
	if (previousX[endpoint] > maximum) {
		maximum = previousX[endpoint];
		state = 1;
	}
	if (previousY[endpoint] > maximum) {
		maximum = previousY[endpoint];
		state = 2;
	}
	if (maximum <= NEGATIVE / 2) return {
		query,
		queryToPanel: new Int32Array(queryLength).fill(-1),
		touchedBand: true
	};
	const queryToPanel = new Int32Array(queryLength);
	queryToPanel.fill(-1);
	let row = panelLength, column = queryLength, touchedBand = false;
	while (row > 0 || column > 0) {
		if (row === 0) {
			queryToPanel[--column] = -1;
			state = 2;
			continue;
		}
		if (column === 0) {
			row--;
			state = 1;
			continue;
		}
		const start = starts[row], end = ends[row];
		if (column < start || column > end) return {
			query,
			queryToPanel,
			touchedBand: true
		};
		if (start > 0 && column <= start + 2 || end < queryLength && column >= end - 2) touchedBand = true;
		const packed = trace[row * stride + column - start];
		if (state === 0) {
			queryToPanel[column - 1] = row - 1;
			row--;
			column--;
			state = packed & 3;
		} else if (state === 1) {
			row--;
			state = packed & 4 ? 1 : 0;
		} else {
			queryToPanel[column - 1] = -1;
			column--;
			state = packed & 8 ? 2 : 0;
		}
	}
	return {
		query,
		queryToPanel,
		touchedBand
	};
}
function alignAdaptively(profile, sequence) {
	const queryLength = sequence.replaceAll("-", "").replace(/\s/g, "").length;
	const maximumBand = Math.max(profile.width, queryLength);
	let band = Math.min(maximumBand, Math.max(64, Math.ceil(Math.sqrt(Math.max(profile.width, queryLength)))));
	while (true) {
		const result = bandedProfileAlignment(profile, sequence, band);
		if (!result.touchedBand || band >= maximumBand) return result;
		band = Math.min(maximumBand, band * 2);
	}
}
/**
* Add one sequence to a fixed alignment without ever realigning its existing
* rows. Query insertions become all-gap columns in the existing profile. This
* is the profile-add operation needed when a biological sample MSA must remain
* invariant while a reference is added for coordinates and scoring.
*/
function addSequenceToProfile(profileRows, sequence) {
	const normalizedRows = profileRows.map((row) => row.replace(/\s/g, "").toUpperCase());
	const profile = buildPanelProfile(normalizedRows), aligned = alignAdaptively(profile, sequence);
	const insertions = Array.from({ length: profile.width + 1 }, () => "");
	const mapped = Array(profile.width).fill("-");
	let previousPanelColumn = -1;
	for (let queryColumn = 0; queryColumn < aligned.query.length; queryColumn++) {
		const panelColumn = aligned.queryToPanel[queryColumn], residue = aligned.query[queryColumn];
		if (panelColumn < 0) insertions[previousPanelColumn + 1] += residue;
		else {
			if (panelColumn <= previousPanelColumn) throw new Error("Profile-add alignment produced non-monotonic columns.");
			mapped[panelColumn] = residue;
			previousPanelColumn = panelColumn;
		}
	}
	const expandProfileRow = (row) => {
		let output = "-".repeat(insertions[0].length);
		for (let column = 0; column < profile.width; column++) output += row[column] + "-".repeat(insertions[column + 1].length);
		return output;
	};
	let alignedSequence = insertions[0];
	for (let column = 0; column < profile.width; column++) alignedSequence += mapped[column] + insertions[column + 1];
	const expandedRows = normalizedRows.map(expandProfileRow);
	if (expandedRows.some((row) => row.length !== alignedSequence.length)) throw new Error("Profile-add alignment did not produce a rectangular alignment.");
	if (alignedSequence.replaceAll("-", "") !== aligned.query) throw new Error("Profile-add alignment did not preserve the added sequence.");
	return {
		profileRows: expandedRows,
		sequence: alignedSequence
	};
}
function maximumSubarrayScore(profile, query, mapping) {
	let current = 0, best = 1;
	for (let index = 0; index < query.length; index++) {
		const panelColumn = mapping[index];
		if (panelColumn < 0) continue;
		const code = query.charCodeAt(index), probability = code < ASCII_SYMBOLS ? profile.probabilities[panelColumn * ASCII_SYMBOLS + code] : 0;
		current += -Math.log(probability + .01) + LOG_QUARTER;
		if (current > best) best = current;
		if (current <= 0) current = 0;
	}
	return best;
}
function filterQueriesAgainstPanel(sequences, panelRows, onProgress) {
	const profile = buildPanelProfile(panelRows), extracted = Array(sequences.length), scores = new Array(sequences.length);
	for (let index = 0; index < sequences.length; index++) {
		const aligned = alignAdaptively(profile, sequences[index]);
		let first = -1, last = -1;
		for (let position = 0; position < aligned.queryToPanel.length; position++) if (aligned.queryToPanel[position] >= 0) {
			if (first < 0) first = position;
			last = position;
		}
		if (first < 0) throw new Error("A query sequence did not overlap the reference panel.");
		extracted[index] = aligned.query.slice(first, last + 1);
		scores[index] = maximumSubarrayScore(profile, aligned.query, aligned.queryToPanel);
		if ((index & 7) === 7 || index + 1 === sequences.length) onProgress?.({
			completed: index + 1,
			total: sequences.length
		});
	}
	return {
		sequences: extracted,
		scores
	};
}
//#endregion
export { filterQueriesAgainstPanel as n, addSequenceToProfile as t };
