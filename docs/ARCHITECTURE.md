# Architecture

webPORPID has one implementation of the read-to-consensus algorithms: a C++20 core compiled as a WASI reactor. Browser workers and `porpid-cli` both instantiate that same `public/webporpid.wasm` artifact.

## Execution flow

1. A FASTQ byte stream is hashed and, for `.gz` input, decompressed incrementally.
2. Bounded 256-record/4 MiB batches are distributed across C++/WASM workers for filtering, primer detection, orientation, demultiplexing, and UMI extraction.
3. Compact binary read records are routed by sample/UMI hash into spool partitions. With a positive `maxReadsPerSample`, a monotone online cutoff prevents records that cannot survive deterministic downsampling from reaching disk.
4. The partitions are scanned for family counts. Sparse UMI offspring likelihoods and LDA assignments form one run-wide family model.
5. The model is copied to every worker and the same selected partition records are scanned for family consensus.
6. Raw-read partitions are closed and deleted. Only consensus and summary data continue into contamination, post-processing, MSA, FastTree, result serialization, and exports.

## Bounded raw-read storage

The spool is the key memory boundary. A record has a fixed 24-byte header containing its complete record length, sample index, and deterministic sampling hash, followed by UMI, read name, sequence, and quality bytes. Count and consensus passes therefore inspect the header before reading a record body.

The final downsampling cutoff for a sample is `floor((2^64 - 1) * maximum / final_count)`. During streaming, webPORPID applies the same expression to the count seen so far. That provisional cutoff can only decrease, so it is guaranteed to retain every record that can pass the final cutoff, independent of worker completion order. Non-candidates bypass the spool, and the browser compacts stale early candidates in place every 512 MiB of admitted data and once at the final cutoff. This preserves the selected read set while preventing a highly compressed multi-million-read FASTQ from expanding into an all-read OPFS spool.

- The CLI writes each partition to an isolated temporary directory and removes the directory after consensus, including on errors.
- A browser analysis uses Origin Private File System synchronous access handles from its dedicated pipeline worker. Short reads and writes are completed in loops, and large inputs make a best-effort persistent-storage request. A real quota/I/O failure reports the run spool size and the browser's remaining origin quota instead of presenting a misleading truncated-frame error.
- If OPFS is unavailable, a bounded in-memory backend is used and stops with a clear error at 512 MiB.

Only one partition per active worker is transferred into WASM at a time. The default 64 hash partitions and default 100,000-read per-sample ceiling keep these working sets substantially smaller than the input dataset. Setting the ceiling to zero intentionally disables adaptive admission and retains every demultiplexed read, so the browser must then have enough origin storage for the fully expanded spool. A single unusually dominant UMI family still necessarily requires its selected reads together for faithful consensus.

## Parallelism and determinism

The default worker count is the runtime's logical CPU count. The browser control and CLI `--workers N` option can reduce it.

Preprocessing batches can finish out of order because spool order is not semantically relevant. Family counts are merged by sample/UMI key, and reads are stably sorted by parsed input ordinal before consensus. Hash-threshold downsampling is seed-, read-, sample-, and UMI-dependent rather than scheduling-dependent. The CLI smoke suite verifies identical consensus output with one and two workers.

## Consensus performance

The consensus hot path uses sparse exact six-mer profiles, collision-free rolling 60-bit encodings for canonical 30-base seeds, contiguous buffers, and seeded long-read alignment. A seedless region uses an adaptive band whose width includes the complete length difference plus a square-root noise allowance; seeded intervals retain the original full dynamic program. Ambiguous/lower-case seed input keeps the byte-exact fallback path. The core is compiled at `-O3 -msimd128` and the optimized binary contains SIMD vector operations.

Centroid profiles are counted once per read rather than allocating and rescanning a dense 4,096-bin vector. Refinement stops at its deterministic fixed point, and the final fixed-point alignments are reused for agreement counting. Agreement symbols use fixed counters instead of per-position strings. Identical families still bypass alignment without changing their sequence or agreement. The substitution, mixed-indel, high-indel, terminal-overhang, and identical-family Julia parity cases remain exact.

## Downstream scale

Alivibe-compatible MSA and FastTree WASM assets are packaged locally; there is no runtime dependency on Swig. Independent per-sample MSA and FastTree jobs use bounded worker pools, capped by both the configured worker count and the number of samples, so a machine with many logical CPUs does not instantiate idle runtimes.

MSA is monolithic up to 8,000 sequences and 128 MiB of input residues. Beyond either threshold, webPORPID aligns deterministic batches of 2,000 rows against a shared anchor, decomposes every batch into anchor insertion slots, merges slot widths, and reconstructs a rectangular alignment in original order. This bounds each aligner's dynamic-programming workspace while retaining the final alignment required for export and visualization.

After non-functional filtering, identical ungapped haplotypes are collapsed. Each representative stores its member UMI-family identifiers; bubble area is therefore proportional to family count rather than read count. FastTree consumes this collapsed nucleotide alignment by default. The uncollapsed alignment is retained for optional on-demand family-level inference. The browser uses the packaged double-precision BioWASM/Aioli FastTree runtime; the CLI invokes the same packaged FastTree WASM directly. Browser tree failure is logged and produces an explicit zero-branch star fallback so a completed analysis remains loadable.

## Result lifecycle

The post-consensus state is validated, MessagePack encoded, gzip compressed, and prefixed with a versioned `.webporpid` magic header. Loading performs size, framing, schema, type, uniqueness, sample-reference, consensus/postproc-consistency, collapse-membership, edit-fingerprint, and alignment/tree-reference checks before the UI receives the bundle. Raw reads and spool records are never serialized. A returned Alivibe/manual edit stores its exact nucleotide FASTA, translation frame, baseline/edited/tree fingerprints, warnings, source/time, and optional recalculated Newick tree as a separate copy in the same project file.

The same result object backs the browser explorer, CLI `inspect`, and all component exports, so exports do not rerun an algorithm or depend on hidden temporary files.
