import { readFile } from "node:fs/promises";

function parse(source) {
  const rows = new Map();
  for (const line of source.trim().split(/\r?\n/)) {
    const [name, minimum, sequence] = line.split("\t");
    if (!name || sequence === undefined) throw new Error(`Invalid parity row: ${line}`);
    rows.set(name, { minimum: Number(minimum), sequence });
  }
  return rows;
}

function editDistance(left, right) {
  let previous = Uint32Array.from({ length: right.length + 1 }, (_, index) => index), current = new Uint32Array(right.length + 1);
  for (let row = 1; row <= left.length; row++) {
    current[0] = row;
    for (let column = 1; column <= right.length; column++) current[column] = Math.min(
      previous[column] + 1, current[column - 1] + 1, previous[column - 1] + (left[row - 1] === right[column - 1] ? 0 : 1),
    );
    [previous, current] = [current, previous];
  }
  return previous[right.length];
}

const [nativePath, juliaPath] = process.argv.slice(2);
if (!nativePath || !juliaPath) throw new Error("Usage: compare-consensus-parity.mjs native.tsv julia.tsv");
const native = parse(await readFile(nativePath, "utf8")), julia = parse(await readFile(juliaPath, "utf8"));
let failed = false;
console.log("case\tedit_distance\tnormalized\tminimum_agreement_delta\tstatus");
for (const [name, expected] of julia) {
  const actual = native.get(name); if (!actual) throw new Error(`Native output is missing ${name}.`);
  const edits = editDistance(expected.sequence, actual.sequence), normalized = edits / Math.max(expected.sequence.length, actual.sequence.length, 1);
  const agreement = Math.abs(expected.minimum - actual.minimum), passed = normalized <= 0.005 && agreement <= 0.03;
  console.log(`${name}\t${edits}\t${normalized.toFixed(6)}\t${agreement.toFixed(3)}\t${passed ? "PASS" : "FAIL"}`);
  failed ||= !passed;
}
if (failed) process.exitCode = 1;
