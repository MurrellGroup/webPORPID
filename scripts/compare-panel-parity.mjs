import { readFile } from "node:fs/promises";

function parse(source) {
  return new Map(source.trim().split(/\r?\n/).filter(Boolean).map((line) => {
    const [caseIndex, row, score, sequence] = line.split("\t"); return [`${caseIndex}/${row}`, { score: Number(score), sequence }];
  }));
}

const [nativePath, juliaPath] = process.argv.slice(2);
if (!nativePath || !juliaPath) throw new Error("Usage: compare-panel-parity.mjs native.tsv julia.tsv");
const native = parse(await readFile(nativePath, "utf8")), julia = parse(await readFile(juliaPath, "utf8"));
let failed = false;
console.log("case\tsequence_equal\tscore_delta\tstatus");
for (const [key, expected] of julia) {
  const actual = native.get(key); if (!actual) throw new Error(`Native output is missing ${key}.`);
  const delta = Math.abs(expected.score - actual.score), equal = expected.sequence === actual.sequence, passed = equal && delta < 1e-10;
  console.log(`${key}\t${equal}\t${delta.toExponential(3)}\t${passed ? "PASS" : "FAIL"}`); failed ||= !passed;
}
if (failed) process.exitCode = 1;
