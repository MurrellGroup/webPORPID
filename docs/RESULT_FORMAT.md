# `.webporpid` result format

A result file is one portable, compressed bundle containing the state needed to inspect or export an analysis after UMI family consensus. It deliberately excludes raw reads, qualities, preprocessing batches, temporary partitions, and aligner working memory.

Schema identifier: `webporpid-results/1`.

## Binary framing

The file consists of:

1. eight magic/version bytes: `57 50 52 00 01 0d 0a 1a`;
2. a gzip member;
3. one MessagePack object inside the gzip member.

The browser and CLI share the same encoder/decoder. The current decoder caps compressed input at 1 GiB and the declared/decompressed payload at 2 GiB, rejects truncated gzip/MessagePack, and deeply validates the decoded object before exposing it.

## Stored components

- **Provenance:** webPORPID version, creation time, engine, worker count, input name and SHA-256, compiled-configuration SHA-256, deterministic seed, and supplied PORPID branch/commit identifier.
- **Configuration and resumable references:** dataset, sample primers, optional `donor_ID`, reference filenames, overrides, every effective pipeline parameter, the contamination panel, and the small per-sample panel/functional-reference records required to resume downstream work. Raw reads and the preprocessing spool are never retained.
- **Quality summary:** total, expected-error, length, primer, sample-ID, BPB, malformed, demultiplexed, downsampled, and per-sample counters, including exact selected/subsampled reads for each sample.
- **UMI family table:** observed UMI, size, most likely parent, posterior/log-offspring probability, disposition, and consensus minimum agreement where available.
- **Consensus records:** identifiers, sample, UMI, family size, nucleotide sequence, minimum agreement, and low-agreement site details.
- **Contamination calls:** one decision per consensus family, with nearest non-self label/distance and primary-discard versus wider-suspect status. A primary decision takes precedence over the wider suspect pass.
- **Post-processing records:** family consensus/aligned sequences, every family-level filter decision, panel score, and optional APOBEC posterior summary. New results do not attach functional decisions to these uncollapsed records.
- **Per-sample products:** collapsed and uncollapsed nucleotide FASTA alignments, reference-clipped codon-aware functional-pass collapsed-variant nucleotide/protein alignments, directly translated protein views, optional Newick trees, and the aligned reference rows used for coordinate display.
- **Collapse membership:** abundance-ranked `sample_vN_K` representative, member identifiers, a multiplicity that counts UMI families rather than reads, and collapsed-only functional decision/trimmed outputs. Minimum agreement remains exclusively attached to the member UMI families.
- **Input mapping/run options/stage state:** exact uploaded-file to YAML-slot assignments, every requested defer option, automatic versus external scratch selection, an explicit completed/deferred/skipped record for contamination, post-processing, collapse, and tree stages, and whether stored post-processing applied or bypassed contamination. Scratch directory handles and paths are never serialized. A consensus-only project is valid; missing post-processing records are not interpreted as passes. Contamination may be incomplete while later stages are complete, but those outputs must be marked bypassed and cannot contain contamination rejects.
- **Optional alignment edits:** exact corrected nucleotide FASTA, translation frame, baseline/edited/tree fingerprints, row/column dimensions, exact added/deleted/changed row identifiers, deterministic minimum-edit per-row and aggregate substitution/insertion/deletion counts, explicit shared-row reorder state, gap-only change flags, warning summary, edit source/time, stale-tree state, and optional recalculated Newick tree. An append-only action history records edits, frame changes, resets, and tree recalculations. The original alignment remains untouched.
- **Interactive threshold decisions:** each checkpoint identifier and phase, acceptance time, accepted global and sample-specific values, and a human-readable change list. The final values are also reflected in the stored effective configuration. Human review pauses are excluded from computational stage and total timing.
- **Run timing/log:** per-stage wall times plus persistent stage counts, input mappings, interactive tree runs, execution/storage choices, and fallbacks.

Validation requires unique sample, family, consensus, and post-processing identities; valid sample indices/references; complete sample summaries; known family dispositions; numeric counters; consensus-linked contamination calls; exact consensus metadata/sequence agreement between consensus and post-processing records; complete family-count collapse membership; reference/alignment width agreement; edit fingerprints; and known-sample alignment/tree keys.

## Component exports

| Component | Contents |
|---|---|
| `consensus-fasta` | Every raw family consensus |
| `passed-consensus-fasta` | Consensus sequences passing artefact, agreement, contamination, and panel gates |
| `rejected-consensus-fasta` | Consensus sequences failing at least one of those gates |
| `trimmed-nt-fasta` | Functionally accepted collapsed variants as trimmed nucleotide sequences |
| `trimmed-aa-fasta` | Functionally accepted collapsed variants as translated amino-acid sequences |
| `family-csv` | UMI parent, family size, posterior/log probability, disposition, and minimum agreement |
| `low-agreement-csv` | Per-site agreement, 3′ coordinate, modal base, and homopolymer run length |
| `contamination-csv` | One primary-or-suspect contamination decision per consensus family |
| `postproc-csv` | All uncollapsed family-level filter flags, panel score, and rejection reasons |
| `apobec-csv` | Stored APOBEC posterior summaries |
| `collapse-csv` | Collapsed representatives, UMI-family multiplicities, functional decisions/reasons, and member identifiers; per-family agreement remains in the post-processing export |
| `nucleotide-alignment` | Selected sample collapsed nucleotide FASTA alignment |
| `protein-alignment` | Direct translation of the selected sample's active collapsed nucleotide alignment and saved frame |
| `newick` | Selected sample active collapsed nucleotide tree, including a recalculated edited tree when present |
| `uncollapsed-nucleotide-alignment` | Every retained UMI-family consensus in its active alignment |
| `uncollapsed-protein-alignment` | Direct translation of the active uncollapsed alignment and saved frame |
| `uncollapsed-newick` | Optional on-demand family-level tree |
| `functional-nucleotide-alignment` | Every functional-filtered collapsed variant in the codon-aware alignment, clipped to the first/last functional-reference nucleotide |
| `functional-protein-alignment` | Direct translation of the active functional nucleotide alignment and saved frame |
| `functional-newick` | Optional on-demand functional-sequence tree |
| `log` | Persistent plain-text run log |

Browser exports always use the sample currently selected in the explorer. If Alivibe or an imported FASTA changed the alignment, alignment/protein/tree exports use that persisted edit. CLI FASTA/CSV exports accept optional `--sample`; alignment and tree exports require it when a result contains multiple samples.

The browser's **Export all (.tar.gz)** control stores the complete editable `.webporpid` project and global run log at the archive root. Every sample-specific component in the table above is placed below a directory named from that sample ID. `cross-sample-overview/` contains sample summary (including donor IDs), input filtering, all effective parameters, optional-stage status, interactive threshold decisions, timings, file mappings, and provenance tables. Empty alignment or Newick members explicitly represent unavailable or deferred components.

```bash
node cli/porpid-cli.mjs export results.webporpid \
  --component trimmed-aa-fasta \
  --sample sample_1 \
  --output sample_1.trimmed-aa.fasta
```

## Versioning

Readers reject unknown schema identifiers instead of guessing. A future incompatible structure should increment both the schema string and the binary magic version, with an explicit migration path if old results need to remain editable.
