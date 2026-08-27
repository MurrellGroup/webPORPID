import { functionalFilter } from "../src/postprocess";

const reference = "ATGCCTTGGGCCATCGGACCATATGTTTACGATGGGCAGCTGACTACCGACAACCGTCAATTCGTCTCAGAGAAGTAA";
const replaceAt = (sequence: string, position: number, base: string) => sequence.slice(0, position) + base + sequence.slice(position + 1);
const cases: Array<[string, string]> = [
  ["exact", reference], ["synonymous", replaceAt(reference, 11, "T")],
  ["codon_deletion", reference.slice(0, 30) + reference.slice(33)],
  ["codon_insertion", reference.slice(0, 33) + "GCT" + reference.slice(33)],
  ["frameshift_deletion", reference.slice(0, 31) + reference.slice(32)],
  ["ambiguous", replaceAt(reference, 19, "N")], ["late_start", "TTG" + reference.slice(3)],
  ["early_stop", reference.slice(0, -3)], ["bad_match", "ATG" + "GCT".repeat(26) + "TAA"],
  ["flanking_sequence", "CC" + reference + "GG"],
];

for (const [name, sequence] of cases) {
  const result = functionalFilter(reference, sequence, 0.7);
  const reason = result.passed ? "pass" : (result.reasons[0]?.split(" ")[0] ?? "reject");
  process.stdout.write([name, String(result.passed), reason, result.nt ?? "", result.aa ?? ""].join("\t") + "\n");
}
