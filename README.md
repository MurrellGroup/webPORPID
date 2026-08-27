# webPORPID

webPORPID is a standalone browser and command-line implementation of the PORPID `nanopore` workflow for long-read amplicon data. The read-to-consensus stages run in a shared C++20 WebAssembly core; TypeScript orchestrates streaming, downstream filtering, alignment, trees, result exploration, and exports.

The repository contains a neutral simulated demo only. Raw reads are processed locally and are never included in a saved result bundle.

## Included pipeline

- Incremental FASTQ and FASTQ.GZ decoding, expected-error and length filtering.
- Primer/orientation detection, sample demultiplexing, BPB/UMI extraction, and deterministic downsampling.
- Sparse two-edit UMI offspring likelihoods, LDA assignment, family-size gates, and heteroduplex detection.
- Indel-tolerant family consensus, minimum-agreement calculation, and low-agreement-site logging.
- Run-aware contamination clustering and filtering.
- Artefact, agreement, panel-profile, and functional filters, plus APOBEC summaries.
- Nucleotide/protein MSA, FastTree inference, and linked interactive alignment/tree exploration.
- A single compressed `.webporpid` result file and fourteen component export types.

## Browser application

Prebuilt WASM assets are committed, so a normal Pages build only needs Node.js 24:

```bash
npm ci
npm run build
```

The static site is written to `dist/`. `.github/workflows/deploy-pages.yml` rebuilds, tests, and deploys it to GitHub Pages. The application accepts either the original single-dataset PORPID YAML shape or webPORPID's editable `dataset`/`samples`/`parameters` shape. In the browser, select the referenced panel, contamination, and optional functional FASTA files alongside the YAML.

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

Reference paths are resolved relative to the configuration file. Workers default to all logical CPUs; use `--workers N` to cap them. Temporary partitions default to the operating-system temporary directory and can be redirected with `WEBPORPID_TMPDIR`.

The release workflow compiles self-contained `porpid-cli` executables for Linux x64/arm64, macOS x64/arm64, and Windows x64 whenever a `v*` tag is pushed or the workflow is started manually.

## Simulated demo

```bash
node cli/porpid-cli.mjs run demo/synthetic_reads.fastq --config demo/synthetic_config.yaml --output demo.synthetic.webporpid --workers 2
node cli/porpid-cli.mjs inspect demo.synthetic.webporpid
node cli/porpid-cli.mjs export demo.synthetic.webporpid --component trimmed-aa-fasta --sample sample_1 --output demo.trimmed-aa.fasta
```

The supplied demo has one accepted UMI family and exercises indel-tolerant consensus and functional trimming. Its expected trimmed protein is `MPWAIGPYVYDGQLTTDNRQFVSEK*`.

## Scale and memory model

Input is decoded in bounded batches. Demultiplexed reads are hashed into disk/OPFS partitions; the count and consensus passes scan fixed-size record headers and materialize only deterministically selected records from one partition per active worker. The CLI always uses temporary files. The browser uses Origin Private File System storage when available and otherwise has an explicit 512 MiB fallback limit.

Downstream MSA runs monolithically up to 8,000 rows and 128 MiB of input bases. Larger jobs use deterministic 2,000-row shared-anchor batches, retaining only the final alignment rather than a single enormous aligner workspace.

## Validation

```bash
npm test
npm run build
```

With Julia 1.10.5 and the supplied validation source tree installed:

```bash
JULIA_BIN=/path/to/julia \
PORPID_JULIA_PROJECT=/path/to/supplied/PORPIDpipeline-source \
npm run parity
```

The recorded parity outcomes and the small set of intentional differences are in [docs/PARITY.md](docs/PARITY.md) and [docs/METHODS.md](docs/METHODS.md).

## Documentation

- [Architecture](docs/ARCHITECTURE.md)
- [Methods and behavioral boundaries](docs/METHODS.md)
- [Julia parity evidence](docs/PARITY.md)
- [Result format and exports](docs/RESULT_FORMAT.md)
- [Possible behavior-changing speedups](improvements.md)
