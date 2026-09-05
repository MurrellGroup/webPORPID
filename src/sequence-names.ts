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

/** Functional variants retain their complete collapsed identity and count. */
export function functionalSequenceName(group: FunctionalNameFields): string {
  if (group.referenceMatch == null) return group.representativeId;
  return `${group.representativeId} rm=${group.referenceMatch.toFixed(2)}`;
}

/** Functional label written by webPORPID 0.3.10–0.3.14; read-only compatibility. */
export function legacyFunctionalSequenceName(group: FunctionalNameFields): string {
  if (group.referenceMatch == null) return group.representativeId;
  const match = /^(.*_v\d+)_\d+$/.exec(group.representativeId);
  return `${match?.[1] ?? group.representativeId} rm=${group.referenceMatch.toFixed(2)}`;
}
