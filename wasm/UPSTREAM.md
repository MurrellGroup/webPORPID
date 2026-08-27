# WebAssembly source provenance

The PORPID algorithms in `src/` are a clean C++20 port of the nanopore branch
at MurrellGroup/PORPIDpipeline commit
`201af7942029cfb7974880e41674be9f0ddfaf3b`. The compatibility target and every
intentional implementation difference are recorded in `docs/METHODS.md` and
`docs/PARITY.md`.

`vendor/alivibe_msa.cpp` is copied from MurrellGroup/swig commit
`e6bced8d1c21f021d3b855783fcfdc4299249b80`. It is the C++20 Alivibe-compatible
MSA path used by Swig, including its scoring and three refinement passes.

The browser supplies parallelism by instantiating the same WASI reactor in
independent Web Workers. The release CLI embeds that same reactor; it is not a
separate reimplementation of the read-to-consensus algorithms.
