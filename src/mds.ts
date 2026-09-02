function pairwiseDistance(left: string, right: string): number {
  let differences = 0, compared = 0;
  for (let index = 0; index < Math.min(left.length, right.length); index += 1) if (left[index] !== "-" && right[index] !== "-") {
    compared += 1; if (left[index] !== right[index]) differences += 1;
  }
  return compared ? differences / compared : 0;
}

/** Classical metric MDS using the two largest eigenpairs of the centered distance matrix. */
export function classicalMds(sequences: string[]): Array<[number, number]> {
  const count = sequences.length; if (count < 2) return sequences.map(() => [0, 0]);
  const squared = Array.from({ length: count }, () => new Float64Array(count)), means = new Float64Array(count); let total = 0;
  for (let left = 0; left < count; left += 1) for (let right = left + 1; right < count; right += 1) {
    const value = pairwiseDistance(sequences[left], sequences[right]) ** 2; squared[left][right] = squared[right][left] = value; means[left] += value; means[right] += value; total += value * 2;
  }
  for (let index = 0; index < count; index += 1) means[index] /= count; total /= count * count;
  const multiply = (vector: Float64Array) => {
    const output = new Float64Array(count);
    for (let row = 0; row < count; row += 1) for (let column = 0; column < count; column += 1)
      output[row] += -.5 * (squared[row][column] - means[row] - means[column] + total) * vector[column];
    return output;
  };
  let spectralShift = 0;
  for (let row = 0; row < count; row += 1) {
    let absoluteRow = 0;
    for (let column = 0; column < count; column += 1) absoluteRow += Math.abs(-.5 * (squared[row][column] - means[row] - means[column] + total));
    spectralShift = Math.max(spectralShift, absoluteRow);
  }
  const vectors: Float64Array[] = [], eigenvalues: number[] = [];
  for (let axis = 0; axis < 2; axis += 1) {
    let vector = Float64Array.from({ length: count }, (_, index) => Math.sin((index + 1) * (axis + 1) * 1.61803398875));
    for (let iteration = 0; iteration < 160; iteration += 1) {
      let mean = vector.reduce((sum, value) => sum + value, 0) / count; for (let index = 0; index < count; index += 1) vector[index] -= mean;
      let next = multiply(vector); for (let index = 0; index < count; index += 1) next[index] += spectralShift * vector[index];
      for (const prior of vectors) { let projection = 0; for (let index = 0; index < count; index += 1) projection += next[index] * prior[index]; for (let index = 0; index < count; index += 1) next[index] -= projection * prior[index]; }
      mean = next.reduce((sum, value) => sum + value, 0) / count; for (let index = 0; index < count; index += 1) next[index] -= mean;
      let normSquared = 0; for (const value of next) normSquared += value * value;
      const norm = Math.sqrt(normSquared); if (!norm) break; for (let index = 0; index < count; index += 1) next[index] /= norm;
      vector = next;
    }
    const product = multiply(vector); let eigenvalue = 0; for (let index = 0; index < count; index += 1) eigenvalue += vector[index] * product[index];
    vectors.push(vector); eigenvalues.push(Math.max(0, eigenvalue));
  }
  return sequences.map((_, index) => [vectors[0][index] * Math.sqrt(eigenvalues[0]), vectors[1][index] * Math.sqrt(eigenvalues[1])]);
}
