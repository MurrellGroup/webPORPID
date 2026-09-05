# webPORPID

webPORPID is a standalone browser and command-line implementation of the PORPID `nanopore` workflow for long-read amplicon data. The read-to-consensus stages run in a shared C++20 WebAssembly core; TypeScript orchestrates streaming, downstream filtering, alignment, trees, result exploration, and exports.

The repository contains a neutral simulated demo only. Raw reads are processed locally and are never included in a saved result bundle.

## Included pipeline

- Incremental FASTQ and FASTQ.GZ decoding, expected-error and length filtering.
- Primer/orientation detection, sample demultiplexing, BPB/UMI extraction, and deterministic downsampling.
- Sparse two-edit UMI offspring likelihoods, LDA assignment, family-size gates, and heteroduplex detection.
- Indel-tolerant family consensus, minimum-agreement calculation, and low-agreement-site logging.
- Run-aware contamination clustering and filtering with three-pass DP-means, exact sparse six-mer kernels, bounded inverted posting indexes, donor-aware self groups, threshold-safe distance pruning, and time-budgeted progress updates for large runs.
- Artefact and agreement filters plus APOBEC summaries. Reference-panel filtering defaults to bundled MAFFT 7.520 FFT-NS-2 and also offers a multi-worker independent-query affine profile aligner for very large family sets.
- Post-filter haplotype collapse with abundance-ranked `sample_vN_K` names, one multiplicity per retained UMI family, collapsed-only functional filtering, collapsed-tree abundance bubbles, and optional on-demand family-level trees. Uncollapsed FASTA/display names are `sample_UMI fs=K minag=0.yy`; their internal join key remains `sample_UMI`.
- Exact abundance geometry: bubble area is strictly linear in retained UMI-family count with no radius cap or floor, and the display slider is expressed as area per family.
- Live, scrollable per-sample demultiplexing counts, a blue/loaded-red/consensus-green read-block grid, and phase-specific feedback with an independent working heartbeat for long operations.
- Dedicated pipeline workers and a best-effort Screen Wake Lock allow work to continue when another tab has focus. The default-checked “Run when not in focus” controls for both the initial run and a post-review rerun add a genuinely audible, very quiet synthesized engine hum while processing; this can reduce Chrome background throttling, but browser/OS suspension still cannot be overridden. Use the CLI for guaranteed unattended runs.
- Browser-history protection and unload warnings while selected inputs, an active run, or loaded results would otherwise be lost.
- Optional contamination, downstream-filtering, collapse, and tree stages can be deferred before a run or skipped while active. Contamination is independently bypassable: downstream work continues without excluding anything at that gate and is explicitly labelled unfiltered. Collapse still requires post-processing, and the default tree still requires collapse.
- An optional consensus-threshold checkpoint after agreement and contamination eligibility are known, with live plots, sliders, unrestricted direct numeric entry, per-sample overrides, and persisted audit records. Completed results can reopen that dialogue and replay only downstream filtering, collapse, and tree stages; if the global cutoffs did not change, only samples with changed effective cutoffs are recomputed and other samples' outputs/edits remain untouched.
- One-click `.tar.gz` export containing the editable project, every sample-prefixed component and static jitter/scatter SVG under its sample-ID directory, and a `cross-sample-overview/` directory of CSV status, parameter, provenance, timing, mapping, and summary tables. A visible export toggle adds balanced-width static phylogram + modal-highlighter SVGs.
- An across-sample sortable overview with demultiplexed and selected-read counts, family- and read/CCS-level percentages for UMI and consensus filters, collapsed-haplotype totals, and functional-pass counts.
- Direct frame-selectable translation plus a Swig-derived linked tree/alignment viewer with prominent nucleotide/amino-acid switching, explicitly applied reference-coordinate regions, modal highlighting, mutation mapping, hideable names, and both tree-only and coordinated tree+alignment SVG exports. Trees open rooted on the zero-length edge to the UMI-family-weighted modal tip, with a topology- and distance-preserving midpoint-root control.
- A contamination workbench with one decision per family and on-demand alignment/tree inference for contamination-panel references, discarded donor contaminants, and retained donor sequences using three categorical tip colors.
- Optional YAML/UI `donor_ID`, same-donor contamination protection, and a donor-level combined collapsed/functional alignment and phylogeny workbench with sample-colored tips.
- A bundled Alivibe pop-out editor with permissive biological-edit warnings, validated return, explicit tree recalculation, separately persisted alignment/frame/tree edits, and a detailed append-only edit audit.
- Interactive UMI, artefact, agreement, MDS/APOBEC, and dinucleotide figures with labelled axes and SVG export.
- Reference-clipped codon-aware nucleotide/protein alignments for every functional-filtered collapsed variant. Sample proteins are aligned first; the reference is then added to that fixed profile without rearranging sample rows and is retained first under its supplied FASTA name. Passing variants are labelled `sample_vN_K rm=0.yy`, retaining the family-count suffix. Protein rows are verified before backtranslation; a sequence-preserving fallback repairs a structurally invalid MSA, while an unrecoverable functional-filter error is isolated to that sample and does not stop other samples. Sequences without a complete start-to-stop ORF are rejected.
- A single compressed `.webporpid` result file, component exports, and a complete gzip-compressed tar bundle.
- A subtle package-derived version label in the page header, so a deployed build can be identified immediately.

## Browser application

Prebuilt WASM assets are committed, so a normal Pages build only needs Node.js 24:

```bash
npm ci
npm run build
```

The static site is written to `dist/`, including a linked Methods index and three detailed topic pages. `.github/workflows/deploy-pages.yml` rebuilds, tests, and deploys it to GitHub Pages. The application accepts either the original single-dataset PORPID YAML shape or webPORPID's editable `dataset`/`samples`/`parameters` shape; either form may add an optional `donor_ID` to a sample. Uploads accumulate across selections and drag/drop operations. Once YAML is present, its panel, contamination, and functional-reference paths become labelled slots; renamed files can be assigned explicitly and the exact mapping is stored in the run log and result file. Current Chromium browsers default to a user-selected external scratch directory, which bypasses browser-origin quota; automatic browser storage remains available as an explicit alternative.

To rebuild both SIMD WASM cores from source, install WASI SDK 25 or newer and run:

```bash
export WASI_SDK=/path/to/wasi-sdk
npm run build:wasm
```

Do not use an HTML preview inside constrained notebook/workspace viewers; use the production build through ordinary static hosting or GitHub Pages.

## Command-line application

The checked-in Node bundle uses the same WASM core and local MSA/FastTree assets:

```bash
node scripts/build-cli.mjs
node cli/porpid-cli.mjs run reads.fastq.gz --config config.yaml --output results.webporpid
node cli/porpid-cli.mjs inspect results.webporpid
node cli/porpid-cli.mjs export results.webporpid --component consensus-fasta --sample sample_1 --output consensus.fasta
```

Reference paths are resolved relative to the configuration file. Workers default to all logical CPUs; use `--workers N` to cap them. Use `--defer-phylogeny` to store collapsed alignments without running FastTree until requested in the browser. Temporary partitions default to the operating-system temporary directory and can be redirected with `WEBPORPID_TMPDIR`.

The reference-panel strategy is stored as `parameters.panelFilterMode` in YAML (`mafft-batch`, the default, or `independent-query`). The CLI can override it for one run with `--panel-filter independent-query`. The independent mode avoids a joint candidate MSA, distributes candidates over the requested workers, and is recommended when a sample contains a very large number of retained UMI families.

The release workflow compiles self-contained `porpid-cli` executables for Linux x64/arm64, macOS x64/arm64, and Windows x64 whenever a `v*` tag is pushed or the workflow is started manually.

## Simulated demo

```bash
node cli/porpid-cli.mjs run demo/synthetic_reads.fastq --config demo/synthetic_config.yaml --output demo.synthetic.webporpid --workers 2
node cli/porpid-cli.mjs inspect demo.synthetic.webporpid
node cli/porpid-cli.mjs export demo.synthetic.webporpid --component trimmed-aa-fasta --sample sample_1 --output demo.trimmed-aa.fasta
```

The supplied demo has one accepted UMI family and exercises indel-tolerant consensus and functional trimming. Its expected trimmed protein is `MPWAIGPYVYDGQLTTDNRQFVSEK*`.

## Scale and memory model

Input is decoded in bounded batches. Demultiplexed reads are hashed into disk partitions; the count and consensus passes scan fixed-size record headers and materialize one partition per active worker. Preprocessing WASM heaps are discarded before consensus, each worker receives only its current partition’s family model, and each completed consensus block overwrites its consumed raw-read partition. Workers are then released and compact result blocks are assembled sequentially instead of accumulating every block in memory. External scratch is the browser default when the directory API is available: the user chooses a writable directory, partitions are streamed sequentially outside origin quota, and the temporary subdirectory is removed after consensus. Automatic OPFS/browser storage remains selectable; its monotone deterministic cutoff bypasses records that cannot survive `maxReadsPerSample`, and periodic compaction removes stale early candidates without changing the final selected set. Count/consensus concurrency is bounded from measured partition sizes. The CLI always uses an ordinary temporary directory; the automatic browser fallback has an explicit 512 MiB in-memory limit if OPFS is unavailable.

The default panel filter runs one MAFFT FFT-NS-2 candidate MSA per sample. The optional independent-query strategy uses adaptive-banded affine alignment to the fixed panel profile and does not allocate an all-candidate MSA. Other downstream MSAs run monolithically up to 8,000 rows and 128 MiB of input bases; larger jobs use deterministic 2,000-row shared-anchor batches.

## Validation

```bash
npm test
npm run build
```

Detailed conformance evidence, behavioral boundaries, and simulated-data performance records are retained in the developer documentation.

## Documentation

- [Architecture](docs/ARCHITECTURE.md)
- [Methods and behavioral boundaries](docs/METHODS.md)
- Browser Methods pages: `methods.html`, `methods-preprocessing.html`, `methods-consensus.html`, and `methods-downstream.html`
- [Read-to-consensus conformance evidence](docs/PARITY.md)
- [Corrected simulated-data benchmarks](docs/BENCHMARKS.md)
- [0.3.8 performance audit](docs/PERFORMANCE_AUDIT.md)
- [Report figure mapping](docs/REPORT_VISUALS.md)
- [Result format and exports](docs/RESULT_FORMAT.md)
- [Possible behavior-changing speedups](improvements.md)
