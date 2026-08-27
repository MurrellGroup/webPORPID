import { extractAndScorePanel } from "../src/panel-profile.ts";

const cases = [
  { panel: ["ACGTACGT", "ACGTTCGT"], sample: ["ACGTACGT", "ACGTTCGT"] },
  { panel: ["AACCGGTT", "AACC-GTT"], sample: ["TTAACCGGTTAA", "TTAACC-GTTAA"] },
  { panel: ["ACGTACGT", "ACGTACGT"], sample: ["GGACGTTACGTCC", "GGACG-TACGTCC"] },
];

cases.forEach(({ panel, sample }, caseIndex) => {
  const result = extractAndScorePanel(sample, panel);
  result.sequences.forEach((sequence, row) => process.stdout.write(`${caseIndex + 1}\t${row + 1}\t${result.scores[row]}\t${sequence}\n`));
});
