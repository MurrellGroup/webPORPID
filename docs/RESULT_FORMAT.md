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
- **Configuration:** dataset, sample primers, reference filenames, overrides, and every pipeline parameter. Reference sequence bodies are not duplicated here.
- **Quality summary:** total, expected-error, length, primer, sample-ID, BPB, malformed, demultiplexed, downsampled, and per-sample counters.
- **UMI family table:** observed UMI, size, most likely parent, posterior/log-offspring probability, disposition, and consensus minimum agreement where available.
- **Consensus records:** identifiers, sample, UMI, family size, nucleotide sequence, minimum agreement, and low-agreement site details.
- **Contamination calls:** nearest non-self label/distance and primary-discard versus wider-suspect status.
- **Post-processing records:** consensus/aligned/trimmed sequences, every filter decision, panel score, functional rejection reasons, and optional APOBEC posterior summary.
- **Per-sample products:** nucleotide FASTA alignments, directly translated protein views, and nucleotide Newick trees.
- **Optional alignment edits:** exact corrected nucleotide FASTA, translation frame, baseline/edited fingerprints, edit source/time, and refreshed Newick tree.
- **Run timing/log:** per-stage wall times plus persistent stage counts, execution/storage choices, and fallbacks.

Validation requires unique sample, family, consensus, and post-processing identities; valid sample indices/references; complete sample summaries; known family dispositions; numeric counters; consensus-linked contamination calls; exact consensus metadata/sequence agreement between consensus and post-processing records; and known-sample alignment/tree keys.

## Component exports

| Component | Contents |
|---|---|
| `consensus-fasta` | Every raw family consensus |
| `passed-consensus-fasta` | Consensus sequences passing artefact, agreement, contamination, and panel gates |
| `rejected-consensus-fasta` | Consensus sequences failing at least one of those gates |
| `trimmed-nt-fasta` | Functionally accepted trimmed nucleotide sequences |
| `trimmed-aa-fasta` | Functionally accepted translated amino-acid sequences |
| `family-csv` | UMI parent, family size, posterior/log probability, disposition, and minimum agreement |
| `low-agreement-csv` | Per-site agreement, 3′ coordinate, modal base, and homopolymer run length |
| `contamination-csv` | Primary and wider-suspect contamination calls |
| `postproc-csv` | All filter flags, panel score, functional flag, and rejection reasons |
| `apobec-csv` | Stored APOBEC posterior summaries |
| `nucleotide-alignment` | Selected sample nucleotide FASTA alignment |
| `protein-alignment` | Direct translation of the selected sample's active nucleotide alignment and saved frame |
| `newick` | Selected sample active nucleotide tree, including a refreshed edited tree when present |
| `log` | Persistent plain-text run log |

Browser exports always use the sample currently selected in the explorer. If Alivibe or an imported FASTA changed the alignment, alignment/protein/tree exports use that persisted edit. CLI FASTA/CSV exports accept optional `--sample`; alignment and tree exports require it when a result contains multiple samples.

```bash
node cli/porpid-cli.mjs export results.webporpid \
  --component trimmed-aa-fasta \
  --sample sample_1 \
  --output sample_1.trimmed-aa.fasta
```

## Versioning

Readers reject unknown schema identifiers instead of guessing. A future incompatible structure should increment both the schema string and the binary magic version, with an explicit migration path if old results need to remain editable.
