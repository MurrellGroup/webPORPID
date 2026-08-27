import { readFile } from "node:fs/promises";

function parse(source) {
  return new Map(source.trim().split(/\r?\n/).map((line) => {
    const [umi, parent, count, probability, disposition] = line.split("\t");
    return [umi, { parent, count: Number(count), probability: Number(probability), disposition: Number(disposition) }];
  }));
}

const [nativePath, juliaPath] = process.argv.slice(2), native = parse(await readFile(nativePath, "utf8")), julia = parse(await readFile(juliaPath, "utf8"));
let failed = false;
console.log("UMI\tparent_equal\tposterior_delta\tdisposition_equal\tstatus");
for (const [umi, expected] of julia) {
  const actual = native.get(umi); if (!actual) throw new Error(`Native output is missing ${umi}.`);
  const delta = Math.abs(expected.probability - actual.probability);
  const passed = expected.parent === actual.parent && expected.count === actual.count && expected.disposition === actual.disposition && delta < 1e-11;
  console.log(`${umi}\t${expected.parent === actual.parent}\t${delta.toExponential(3)}\t${expected.disposition === actual.disposition}\t${passed ? "PASS" : "FAIL"}`);
  failed ||= !passed;
}
if (failed) process.exitCode = 1;
