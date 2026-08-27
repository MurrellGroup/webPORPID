type Column = Map<string, number>;

function profile(rows: readonly string[]) {
  if (!rows.length) return [] as Column[];
  const width = rows[0].length;
  if (!width || rows.some((row) => row.length !== width)) throw new Error("Panel and sample profiles require rectangular alignments.");
  return Array.from({ length: width }, (_, column) => {
    const counts = new Map<string, number>();
    for (const row of rows) { const value = row[column].toUpperCase(); counts.set(value, (counts.get(value) ?? 0) + 1); }
    return new Map([...counts].map(([value, count]) => [value, count / rows.length]));
  });
}

function cost(left: Column, right: Column) {
  let equal = 0;
  for (const [base, probability] of left) equal += probability * (right.get(base) ?? 0);
  return 2 * equal - 1;
}

function alignProfiles(panel: Column[], sample: Column[]) {
  const rows = panel.length + 1, columns = sample.length + 1, cells = rows * columns;
  if (!Number.isSafeInteger(cells) || cells > 512 * 1024 * 1024)
    throw new Error("The panel/profile alignment is too large for browser memory.");
  const trace = new Uint8Array(cells), gapOpen = -2, gapExtend = -0.2, negative = Number.NEGATIVE_INFINITY;
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
      if (fromX > bestM) { bestM = fromX; fromM = 1; }
      if (fromY > bestM) { bestM = fromY; fromM = 2; }
      currentM[column] = bestM;
      const openX = previousM[column] + gapOpen, extendX = previousX[column] + gapExtend;
      const xContinues = extendX > openX; currentX[column] = xContinues ? extendX : openX;
      const openY = currentM[column - 1] + gapOpen, extendY = currentY[column - 1] + gapExtend;
      const yContinues = extendY > openY; currentY[column] = yContinues ? extendY : openY;
      trace[row * columns + column] = fromM | (xContinues ? 4 : 0) | (yContinues ? 8 : 0);
    }
    [previousM, currentM] = [currentM, previousM];
    [previousX, currentX] = [currentX, previousX];
    [previousY, currentY] = [currentY, previousY];
  }
  let state = 0, maximum = previousM[columns - 1];
  if (previousX[columns - 1] > maximum) { maximum = previousX[columns - 1]; state = 1; }
  if (previousY[columns - 1] > maximum) state = 2;
  let row = panel.length, column = sample.length;
  const reverse: number[] = [];
  while (row > 0 && column > 0) {
    const packed = trace[row * columns + column];
    if (state === 0) { reverse.push(row - 1); row--; column--; state = packed & 3; }
    else if (state === 1) { row--; state = packed & 4 ? 1 : 0; }
    else { reverse.push(-1); column--; state = packed & 8 ? 2 : 0; }
  }
  while (column-- > 0) reverse.push(-1);
  reverse.reverse();
  if (reverse.length !== sample.length) throw new Error("The panel/profile alignment produced an invalid sample map.");
  return reverse;
}

function maximumSubarray(values: number[]) {
  let current = 0, best = 1;
  for (const value of values) { current += value; best = Math.max(best, current); if (current <= 0) current = 0; }
  return best;
}

export function extractAndScorePanel(sampleRows: readonly string[], panelRows: readonly string[]) {
  if (!sampleRows.length) return { sequences: [] as string[], scores: [] as number[] };
  const sampleProfile = profile(sampleRows), panelProfile = profile(panelRows), sampleToPanel = alignProfiles(panelProfile, sampleProfile);
  const start = sampleToPanel.findIndex((index) => index >= 0);
  let end = sampleToPanel.length - 1; while (end >= 0 && sampleToPanel[end] < 0) end--;
  if (start < 0 || end < 0) throw new Error("The sample alignment does not overlap the reference panel.");
  const sequences = sampleRows.map((row) => row.slice(start, end + 1));
  const scores = sequences.map((sequence) => {
    const transformed: number[] = [];
    for (let offset = 0; offset < sequence.length; offset++) {
      const base = sequence[offset].toUpperCase(), panelIndex = sampleToPanel[start + offset];
      if (base === "-" || panelIndex < 0) continue;
      const probability = panelProfile[panelIndex].get(base) ?? 0;
      transformed.push(-Math.log(probability + 0.01) + Math.log(0.25));
    }
    return maximumSubarray(transformed);
  });
  return { sequences, scores };
}
