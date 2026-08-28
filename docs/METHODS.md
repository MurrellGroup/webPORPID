# Methods and behavioral boundaries

The preprocessing, UMI grouping, and family-consensus stages are a clean C++20 port of the supplied PORPID `nanopore` source. The implementation changes data structures and execution strategy where that does not change the intended calculation. Downstream analysis is also implemented locally, with exact panel-profile parity and the documented differences below.

## FASTQ filtering and demultiplexing

- FASTQ and FASTQ.GZ are parsed as streams; record structure and sequence/quality lengths are checked before WASM receives a batch.
- Mean per-read error is computed from Phred qualities. A read passes only when error is below `errorRate` and length is strictly between `minLength` and `maxLength`, matching the supplied code's boundary tests.
- Second-strand and conserved cDNA adapter primers are searched in both orientations using the nanopore branch's bounded-window, error-tolerant matching and non-overlapping last-hit behavior.
- Lower-case bases in `cDNA_primer` define the sample ID. `N` bases define the UMI. The BPB/tag dynamic program preserves the original emission/gap scores, traceback priorities, and more-than-four-error reject boundary.
- Passing records retain the oriented sequence and quality needed by heteroduplex and consensus logic.

If a sample exceeds `maxReadsPerSample`, the supplied implementation retains each read with probability `maxReadsPerSample / observedReads` ("about" the requested maximum). webPORPID preserves that Bernoulli rule but replaces Julia's process-global RNG with a deterministic 64-bit hash threshold. The retained subset can therefore differ from a Julia run, while its inclusion probability is the same and its identity is reproducible across worker counts.

## UMI offspring model

For each sample, webPORPID constructs likelihoods only between observed UMI strings reachable by the original two-error insertion/deletion/substitution process. Identity probabilities, error weights, Dirichlet concentration, iterative LDA update, convergence criterion, posterior parent assignment, and decision ordering are preserved.

The likelihood matrix is stored as sparse rows instead of the original dense/allocation-heavy intermediates. This is a representation optimization: the parity suite finds identical parent and disposition assignments in all seven tested one- and two-edit cases, with a maximum posterior delta of `1.11e-16`.

The decision sequence is LDA posterior, eight-base UMI length, then configured family-size threshold/override. BPB rejects are retained in the family report. The first 50 quality positions are used for the same heteroduplex t-test and minimum-barcode-quality condition before an otherwise accepted family enters consensus.

## Indel-tolerant family consensus

The consensus path does not assume low indel rates.

1. A six-mer count centroid selects a real family read as the initial reference.
2. Exact 30-base seed anchors divide long alignments where they are unique and collinear. Seeded intervals use the supplied global scoring and edge-gap behavior. If no consistent seed exists, an adaptive band covers the full terminal length difference plus a square-root noise allowance; it falls back to the complete matrix if the endpoint cannot be reached.
3. Three reference-refinement passes identify poorly supported adjacent positions, replace uncertain blocks with modal aligned substrings, and extend both ends by majority-supported iterative alignment.
4. A seeded final alignment, matching `RobustAmpliconDenoising.get_matches`, calculates per-position agreement. If refinement has reached a fixed point, its exact alignments are reused. `minimumAgreement` is the rounded minimum over interior positions, and every tied low-agreement position records its 3′-relative coordinate, rounded agreement, complemented modal base, and modal homopolymer run length.
5. The cDNA/UMI-side primer remainder is aligned and removed, then the consensus is reverse-complemented into output orientation.

Sparse six-mer profiles retain exactly the observed bins, canonical 30-mers use collision-free rolling integers, and byte-identical families take a direct path. The shortcut returns the same sequence and agreement of one. Fixed-point detection skips refinement passes that cannot change the result and reuses their alignments for agreement. These changes remove repeated dense allocations and redundant dynamic programs without introducing an edit-distance consensus shortcut.

Five parity families covering substitutions, mixed indels, high indel load, terminal overhangs, and identical reads have zero consensus edit distance and zero minimum-agreement delta against Julia. On the corrected opaque simulated run, 3,040 of 3,043 consensus sequences are byte-exact; 3,040 minimum-agreement values are within `0.01` and 2,993 are exact.

## Contamination filtering

Consensus sequences are converted to raw six-mer count vectors. A DP-means pass builds per-sample run-derived clusters; clusters above the configured proportion and an all-sequence center form the primary database. A second zero-proportion database records wider suspect calls. Both are followed by the external contamination panel so equal-distance label selection has the same ordering as Julia.

Sparse integer vectors are used for individual sequences and dense vectors only for cluster means. Squared distances and mean centers are algebraically the same as dense 4,096-element vectors.

One intentional reproducibility difference remains: Julia resolves ambiguous IUPAC bases through its global RNG. webPORPID resolves them with a stable sequence/seed-derived generator. Unambiguous sequences are unaffected; ambiguous sequences can fall on a different side of a clustering or distance threshold.

## Post-processing

- The artefact cutoff is `ceil(quantile(non-contaminant family sizes, q) * artefactFraction)` with per-sample overrides.
- Agreement and contamination decisions are retained separately for audit and export.
- Candidate sequences are aligned, converted to a sample profile, and affine-profile-aligned to the supplied reference panel. Profile column cost, gap-open/extend behavior, traceback priority, extracted coordinates, indel-ignoring probability, and maximum-subarray panel score are ported from Julia. All six rows in the three-case panel parity suite match exactly in sequence and score.
- The functional filter finds candidate coding regions, translates them, performs one joint amino-acid MSA with the reference, and projects alignment gaps back as complete codon triplets. It then applies start/stop, ambiguity, and reference-match gates and stores trimmed nucleotide/amino-acid outputs. This avoids arbitrary one-base gaps from independent nucleotide alignments. All ten focused parity cases classify identically; on the corrected long-read run the aggregate passing count is 1,185 versus Julia's 1,180 among 1,455 evaluated (`0.34%` count difference).
- APOBEC summaries use the local four-state mutation grid and posterior integration; they are reporting fields rather than a filtering gate.
- Accepted nucleotide sequences are aligned with the packaged MSA and nucleotide trees use packaged FastTree. Protein exploration does not depend on the optional functional filter: it directly translates every active aligned nucleotide row, with `--- → -` and mixed-gap/ambiguous/incomplete codons → `X`, using a persisted frame offset of 0, 1, or 2.

## Large-alignment batching difference

Below the 8,000-row/128 MiB threshold, downstream MSA uses the monolithic packaged aligner. Larger inputs use shared-anchor batches. Every output retains row order and ungapped residues, but equally scoring indel placements in repetitive regions can differ from a single monolithic MSA. This affects downstream visualization/tree columns, not preprocessing, UMI classification, or the stored raw family consensus. A true streaming profile/POA alternative is recorded in `improvements.md` because changing it requires a dedicated downstream parity study.

## Stable intentional differences

- Deterministic hash-based Bernoulli downsampling replaces Julia's global random draws.
- Ambiguous-base contamination resolution is deterministic.
- The seedless consensus fallback uses an adaptive band; this can choose a different equal-score path from a full matrix in an unusually repetitive region, although all focused high-indel parity families are exact.
- Worker scheduling cannot change family or read order at consensus time.
- Large MSA uses shared-anchor batching above its documented threshold.
- Four rejected functional cases can use a different failure-stage label although classification and accepted outputs match.
- Reporting is one versioned `.webporpid` bundle rather than generated HTML and intermediate directories.

Potential speedups expected to change more substantive behavior are deliberately excluded and tracked in `improvements.md`.
