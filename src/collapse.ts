import { inspectAlignment } from "./alignment-utils.ts";
import type { CollapseGroup, PostprocRecord } from "./types";

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
 * collapse multiplicity.
 */
export function collapseAlignment(fasta: string, sample: string, records: readonly PostprocRecord[]): CollapsedAlignment {
  const alignment = inspectAlignment(fasta, 1);
  const metadata = new Map(records.map((record) => [record.id, record]));
  const bySequence = new Map<string, { name: string; sequence: string; members: string[]; agreements: number[] }>();
  for (const row of alignment.records) {
    // A haplotype is defined by nucleotide content. Equivalent residues can
    // acquire different but equally scoring gap placements in a large MSA;
    // those are not distinct biological variants.
    const haplotype = row.sequence.replaceAll("-", "");
    let group = bySequence.get(haplotype);
    if (!group) {
      group = { name: row.name, sequence: row.sequence, members: [], agreements: [] };
      bySequence.set(haplotype, group);
    }
    group.members.push(row.name);
    const agreement = metadata.get(row.name)?.minimumAgreement;
    if (agreement != null && Number.isFinite(agreement)) group.agreements.push(agreement);
  }
  const groups: CollapseGroup[] = [...bySequence.values()].map((group) => ({
    sample,
    representativeId: group.name,
    memberIds: group.members,
    familyCount: group.members.length,
    minimumAgreement: group.agreements.length ? Math.min(...group.agreements) : 0,
  }));
  return {
    fasta: exactFasta([...bySequence.values()].map((group) => ({ name: group.name, sequence: group.sequence }))),
    groups,
  };
}
