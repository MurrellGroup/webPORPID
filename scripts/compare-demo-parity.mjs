import { readFile } from "node:fs/promises";
import { decodeResult } from "../src/result-file.ts";

const [juliaPath, resultPath] = process.argv.slice(2);
if (!juliaPath || !resultPath) throw new Error("Usage: compare-demo-parity.mjs julia.tsv result.webporpid");
const lines = (await readFile(juliaPath, "utf8")).trim().split(/\r?\n/), bundle = decodeResult(new Uint8Array(await readFile(resultPath)));
const counts = lines.find((line) => line.startsWith("counts\t"))?.split("\t").slice(1).map(Number);
const family = lines.find((line) => line.startsWith("family\t"))?.split("\t");
const consensus = lines.find((line) => line.startsWith("consensus\t"))?.split("\t");
if (!counts || !family || !consensus) throw new Error("Julia demo output is incomplete.");
const qualityKeys = ["totalReads", "qualityReads", "badReads", "shortReads", "longReads", "demultiplexedReads"];
const checks = qualityKeys.map((key, index) => [`${key.replace(/[A-Z]/g, (letter) => ` ${letter.toLowerCase()}`)}`, bundle.quality[key] === counts[index]]);
const actualFamily = bundle.umiFamilies.find((row) => row.umi === family[1]), actualConsensus = bundle.consensuses.find((row) => row.umi === consensus[1]);
checks.push(["UMI family size", actualFamily?.familySize === Number(family[2])]);
checks.push(["UMI posterior", Math.abs((actualFamily?.posteriorProbability ?? Number.NaN) - Number(family[3])) < 1e-11]);
checks.push(["minimum agreement", Math.abs((actualConsensus?.minimumAgreement ?? Number.NaN) - Number(consensus[2])) < 1e-12]);
checks.push(["consensus sequence", actualConsensus?.sequence === consensus[3]]);
const record = bundle.records.find((row) => row.id === actualConsensus?.id);
const expectedProtein = "MPWAIGPYVYDGQLTTDNRQFVSEK*";
checks.push(["functional filter", record?.functionalPass === true]);
checks.push(["trimmed nucleotide", record?.trimmedNt === consensus[3]]);
checks.push(["trimmed amino acid", record?.trimmedAa === expectedProtein]);
checks.push(["protein alignment", bundle.alignments["sample_1/protein"]?.includes(expectedProtein) === true]);
checks.push(["nucleotide tree", bundle.trees["sample_1/nucleotide"]?.includes(actualConsensus?.id.replace(/[^A-Za-z0-9_.|*+\-]/g, "_") ?? "") === true]);
let failed = false; console.log("check\tstatus");
for (const [name, passed] of checks) { console.log(`${name}\t${passed ? "PASS" : "FAIL"}`); failed ||= !passed; }
if (failed) process.exitCode = 1;
