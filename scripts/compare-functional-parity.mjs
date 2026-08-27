import { readFile } from "node:fs/promises";

function records(source) {
  return new Map(source.trim().split(/\r?\n/).filter(Boolean).map((line) => {
    const [name, passed, reason, nt = "", aa = ""] = line.split("\t"); return [name, { passed: passed === "true", reason, nt, aa }];
  }));
}

const [nativePath, juliaPath] = process.argv.slice(2);
if (!nativePath || !juliaPath) throw new Error("Usage: compare-functional-parity.mjs native.tsv julia.tsv");
const native = records(await readFile(nativePath, "utf8")), julia = records(await readFile(juliaPath, "utf8"));
let failed = false;
console.log("case\tclassification\treason\ttrimmed_nt\ttrimmed_aa\tstatus");
for (const [name, expected] of julia) {
  const actual = native.get(name); if (!actual) throw new Error(`Native output is missing ${name}.`);
  const classification = actual.passed === expected.passed, reason = actual.reason === expected.reason;
  const nt = !expected.passed || actual.nt === expected.nt, aa = !expected.passed || actual.aa === expected.aa;
  const status = classification && nt && aa ? "PASS" : "FAIL"; failed ||= status === "FAIL";
  console.log([name, classification, reason, nt, aa, status].join("\t"));
}
if (failed) process.exitCode = 1;
