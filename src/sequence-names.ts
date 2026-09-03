/** The stable label used for the first row of every functional alignment. */
export const FUNCTIONAL_REFERENCE_NAME = "functional_reference";

interface UncollapsedNameFields {
  id: string;
  familySize: number;
  minimumAgreement: number;
}

interface FunctionalNameFields {
  representativeId: string;
  referenceMatch?: number;
}

/**
 * Human-readable UMI-family FASTA label. The internal consensus identifier
 * remains `sample_UMI`, so joins and contamination decisions never depend on
 * presentation metadata embedded in a FASTA header.
 */
export function uncollapsedSequenceName(record: UncollapsedNameFields): string {
  return `${record.id} fs=${record.familySize} minag=${record.minimumAgreement.toFixed(2)}`;
}

/** Recover the stable internal identifier from a webPORPID UMI-family label. */
export function uncollapsedSequenceId(name: string): string {
  return name.replace(/\s+fs=\d+\s+minag=[^\s]+$/, "");
}

/**
 * Functional variants retain sample and abundance-rank identity but replace
 * the collapsed family-count suffix with their reference-match annotation.
 */
export function functionalSequenceName(group: FunctionalNameFields): string {
  if (group.referenceMatch == null) return group.representativeId;
  const match = /^(.*_v\d+)_\d+$/.exec(group.representativeId);
  const stem = match?.[1] ?? group.representativeId;
  return `${stem} rm=${group.referenceMatch.toFixed(2)}`;
}
