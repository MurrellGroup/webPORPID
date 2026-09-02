import { inspectAlignment } from "./alignment-utils.ts";
import type { CollapseGroup } from "./types";

export interface CollapsedAlignment {
  fasta: string;
  groups: CollapseGroup[];
}

function exactFasta(rows: Array<{ name: string; sequence: string }>): string {
  return rows.map((row) => `>${row.name}\n${row.sequence}\n`).join("");
}

/**
 * Collapse identical aligned haplotypes while retaining one count per UMI
 * family. Read-family sizes are metadata only and never contribute to the
 * collapse multiplicity. Per-family agreement is deliberately not projected
 * onto the resulting haplotype.
 */
export function collapseAlignment(fasta: string, sample: string): CollapsedAlignment {
  const alignment = inspectAlignment(fasta, 1);
  const bySequence = new Map<string, { haplotype: string; sequence: string; members: string[] }>();
  for (const row of alignment.records) {
    // A haplotype is defined by nucleotide content. Equivalent residues can
    // acquire different but equally scoring gap placements in a large MSA;
    // those are not distinct biological variants.
    const haplotype = row.sequence.replaceAll("-", "");
    let group = bySequence.get(haplotype);
    if (!group) {
      group = { haplotype, sequence: row.sequence, members: [] };
      bySequence.set(haplotype, group);
    }
    group.members.push(row.name);
  }
  // Variant numbers are scientific identifiers, so make both abundance order
  // and tie-breaking deterministic. Multiplicity is UMI-family count only.
  const variants = [...bySequence.values()].sort((left, right) =>
    right.members.length - left.members.length || left.haplotype.localeCompare(right.haplotype));
  const groups: CollapseGroup[] = variants.map((group, index) => ({
    sample, representativeId: `${sample}_v${index + 1}_${group.members.length}`,
    memberIds: group.members, familyCount: group.members.length,
  }));
  return {
    fasta: exactFasta(variants.map((group, index) => ({ name: groups[index].representativeId, sequence: group.sequence }))),
    groups,
  };
}
