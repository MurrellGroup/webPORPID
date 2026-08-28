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
  const bySequence = new Map<string, { name: string; sequence: string; members: string[] }>();
  for (const row of alignment.records) {
    // A haplotype is defined by nucleotide content. Equivalent residues can
    // acquire different but equally scoring gap placements in a large MSA;
    // those are not distinct biological variants.
    const haplotype = row.sequence.replaceAll("-", "");
    let group = bySequence.get(haplotype);
    if (!group) {
      group = { name: row.name, sequence: row.sequence, members: [] };
      bySequence.set(haplotype, group);
    }
    group.members.push(row.name);
  }
  const groups: CollapseGroup[] = [...bySequence.values()].map((group) => ({
    sample,
    representativeId: group.name,
    memberIds: group.members,
    familyCount: group.members.length,
  }));
  return {
    fasta: exactFasta([...bySequence.values()].map((group) => ({ name: group.name, sequence: group.sequence }))),
    groups,
  };
}
