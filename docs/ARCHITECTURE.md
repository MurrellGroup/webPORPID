# Architecture

webPORPID has one implementation of the read-to-consensus algorithms: a C++20 core compiled as a WASI reactor. Browser workers and `porpid-cli` both instantiate that same `public/webporpid.wasm` artifact.

## Execution flow

1. A FASTQ byte stream is hashed and, for `.gz` input, decompressed incrementally.
2. Bounded 256-record/4 MiB batches are distributed across C++/WASM workers for filtering, primer detection, orientation, demultiplexing, and UMI extraction.
3. Compact binary read records are routed by sample/UMI hash into spool partitions.
4. The partitions are scanned for family counts. Sparse UMI offspring likelihoods and LDA assignments form one run-wide family model.
5. The model is copied to every worker and the same selected partition records are scanned for family consensus.
6. Raw-read partitions are closed and deleted. Only consensus and summary data continue into contamination, post-processing, MSA, FastTree, result serialization, and exports.

## Bounded raw-read storage

The spool is the key memory boundary. A record has a fixed 24-byte header containing its complete record length, sample index, and deterministic sampling hash, followed by UMI, read name, sequence, and quality bytes. Count and consensus passes therefore inspect the header before reading a record body. Records excluded by `maxReadsPerSample` are never materialized again.

- The CLI writes each partition to an isolated temporary directory and removes the directory after consensus, including on errors.
- A browser analysis uses Origin Private File System synchronous access handles from its dedicated pipeline worker.
- If OPFS is unavailable, a bounded in-memory backend is used and stops with a clear error at 512 MiB.

Only one partition per active worker is transferred into WASM at a time. The default 64 hash partitions and default 100,000-read per-sample ceiling keep these working sets substantially smaller than the input dataset. A single unusually dominant UMI family still necessarily requires its selected reads together for faithful consensus.

## Parallelism and determinism

The default worker count is the runtime's logical CPU count. The browser control and CLI `--workers N` option can reduce it.

Preprocessing batches can finish out of order because spool order is not semantically relevant. Family counts are merged by sample/UMI key, and reads are stably sorted by parsed input ordinal before consensus. Hash-threshold downsampling is seed-, read-, sample-, and UMI-dependent rather than scheduling-dependent. The CLI smoke suite verifies identical consensus output with one and two workers.

## Consensus performance

The consensus core uses compact 2-bit six-mer codes, contiguous buffers, sparse UMI likelihood rows, seeded long-read alignment, and WASM SIMD (`simd128`) for six-mer distance accumulation. Six-mer work arrays live on the heap to avoid the small WASI stack. Identical families bypass alignment without changing their consensus or minimum agreement. Non-identical and high-indel families retain the unbanded alignment/refinement behavior needed for parity.

## Downstream scale

Alivibe-compatible MSA and FastTree WASM assets are packaged locally; there is no runtime dependency on Swig.

MSA is monolithic up to 8,000 sequences and 128 MiB of input residues. Beyond either threshold, webPORPID aligns deterministic batches of 2,000 rows against a shared anchor, decomposes every batch into anchor insertion slots, merges slot widths, and reconstructs a rectangular alignment in original order. This bounds each aligner's dynamic-programming workspace while retaining the final alignment required for export and visualization.

FastTree consumes the final nucleotide alignment. The browser uses the packaged BioWASM/Aioli runtime; the CLI invokes the same packaged FastTree WASM directly. Browser tree failure is logged and produces an explicit zero-branch star fallback so a completed analysis remains loadable.

## Result lifecycle

The post-consensus state is validated, MessagePack encoded, gzip compressed, and prefixed with a versioned `.webporpid` magic header. Loading performs size, framing, schema, type, uniqueness, sample-reference, consensus/postproc-consistency, and alignment/tree-reference checks before the UI receives the bundle. Raw reads and spool records are never serialized.

The same result object backs the browser explorer, CLI `inspect`, and all component exports, so exports do not rerun an algorithm or depend on hidden temporary files.
