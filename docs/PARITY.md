# Julia parity

Parity was measured with Julia 1.10.5 against the source tree supplied with this project handoff. The validation environment used its pinned manifest, including PORPID.jl, RobustAmpliconDenoising.jl, and SeededAlignment.jl. The original PORPIDpipeline repository was not retrieved online.

## Recorded results

| Suite | Cases checked | Result |
|---|---:|---|
| Family consensus | 5 families | All consensus edit distances `0`; all normalized distances `0`; all minimum-agreement deltas `0.000` |
| UMI offspring/LDA | 7 observed UMIs | All parent assignments and dispositions equal; maximum posterior delta `1.11e-16` |
| Panel profile/extraction | 3 alignments, 6 output rows | Every extracted sequence equal; every panel-score delta `0` |
| Functional filter | 10 sequences | Accept/reject classification `10/10`; accepted trimmed nucleotide/protein outputs exact |
| End-to-end simulated demo | 15 assertions | Read counters, family size/posterior, minimum agreement, consensus, functional trimming, protein alignment, and tree all pass |

The consensus cases cover identical reads, a substitution split, mixed insertion/deletion noise, a high-indel family, and terminal overhangs. The UMI cases cover identity and one-/two-edit insertion, deletion, and substitution neighborhoods.

## Corrected six-sample simulated run

The larger supplied dataset was treated as an opaque fixture. Julia and webPORPID produced the same 3,043 consensus keys. Consensus sequences are exact for 3,040 (99.90%); minimum agreement is exact for 2,993 and within `0.01` for 3,040. Contamination pass count and all six per-sample non-functional postprocessing pass counts are exact. The codon-aware functional batch filter passes 1,185 versus Julia's 1,180 among 1,455 evaluated sequences (`0.34%` aggregate count difference). Full stage timings and counts are in [BENCHMARKS.md](BENCHMARKS.md).

In functional parity, the `frameshift_deletion`, `late_start`, `early_stop`, and `bad_match` cases are rejected by both implementations but expose a different first/stage diagnostic label. Classification is the parity gate; reason-label equality is printed separately so this difference remains visible.

## Reproduce

Build the current native and WASM artifacts first:

```bash
export WASI_SDK=/path/to/wasi-sdk-25-or-newer
npm ci
npm run build:wasm
npm test
```

Then run every Julia comparison:

```bash
JULIA_BIN=/path/to/julia-1.10.5/bin/julia \
PORPID_JULIA_PROJECT=/path/to/supplied/PORPIDpipeline-source \
npm run parity
```

Individual suites are available as:

```bash
bash scripts/parity-consensus.sh
bash scripts/parity-umi.sh
bash scripts/parity-panel.sh
bash scripts/parity-functional.sh
bash scripts/parity-demo.sh
```

`JULIA_BIN` and `PORPID_JULIA_PROJECT` must be exported for those individual commands. Temporary parity files are made with `mktemp` and removed when each script exits.

## Non-Julia regression coverage

`npm test` additionally checks:

- native C++ preprocessing, tag extraction, UMI logic, consensus, and a non-unanimous agreement case;
- the optimized WASI reactor through a complete in-memory preprocess → count → model → consensus run;
- malformed spool-record rejection and header-only downsampling selection;
- deterministic large-MSA batching with 8,001 rows, including rectangularity, row order, and ungapped-residue preservation;
- result round-trip, export, framing rejection, truncation rejection, and consensus/postproc consistency validation;
- direct aligned translation, frame-aware protein export, Swig branch scaling, Alivibe bridge revision/snapshot/staleness guards, and persisted edited alignment/frame/tree round-trip;
- a gzipped CLI analysis, result inspection, FASTA export, and identical consensus with one versus two workers.

The parity fixtures are neutral synthetic strings. The bundled demo is simulated and contains no original biological demo material.
