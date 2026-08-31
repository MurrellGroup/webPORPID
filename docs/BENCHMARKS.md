# Corrected simulated-data benchmark

The 0.3.8 whole-pipeline hot-path audit and reproducible synthetic benchmark are in [PERFORMANCE_AUDIT.md](PERFORMANCE_AUDIT.md).

This benchmark uses the corrected six-sample simulated dataset supplied for the 2026-08-27 validation. The biological strings were treated as opaque fixtures: the comparison reads pipeline outputs and aggregate counts, not sequence content. The run contained 59,695 FASTQ records, 59,432 demultiplexed reads, 3,454 observed UMI families, and 3,043 consensus calls.

## Version 0.3.7 contamination stress check

A separate in-memory fixture was generated locally from 12 random 600-base founders, bounded substitutions, six samples, 6,000 consensus families, and 80 panel references. It contains no uploaded data. On Node 24 in the same process environment, the 0.3.6 contamination implementation took 55.667 s and the final 0.3.7 implementation took 3.920 s (14.2× faster). The complete, sorted contamination-call JSON—including cluster/panel label, discard/suspect decision, and floating distance—was identical. A second 320-family differential containing IUPAC and unknown symbols also matched every call and distance exactly.

This is a focused algorithm stress check rather than an end-to-end prediction for every run. Speedup depends on sample count, sequence length, and the number/occupancy of DP-means clusters; the temporary inverted index has a fixed posting cap and deliberately falls back to the exact sparse scan instead of growing without bound.

## Measurement interpretation

- Julia 1.10.5 ran the supplied local nanopore pipeline through Snakemake with five cores. No online PORPIDpipeline source was used.
- “Julia algorithm” is the timer printed inside the supplied Julia stage. Per-sample stages list the six individual values because those processes overlap.
- “Julia workflow wall” is the observed interval for the corresponding Snakemake wave and includes Julia process/package startup and scheduling. This is useful operationally but is not a pure kernel benchmark.
- webPORPID is the checked-in Node/WASI CLI with five workers. Its stored stage times are aggregate wall times around the complete stage, including all six samples.
- The Julia run completed every analysis, all 30 PNG figures, all six tree SVGs, all six HTML reports, and the index. Its final legacy tar-packaging rule then failed on an unrelated invalid `Base.-` method definition; the reported analysis timings stop before that packaging failure.

## Stage timings

| Stage | Julia algorithm timer | Julia workflow wall | webPORPID, 5 workers | Comparable result |
|---|---:|---:|---:|---|
| Setup | included in each process | — | 0.131 s | web-only isolated setup |
| FASTQ filter + demux | 39.975 s | 60 s | 4.051 s | 9.87× vs Julia algorithm |
| UMI grouping/model | 26.243–30.270 s per sample | 77 s | 2.496 s | 10.5–12.1× vs one Julia sample timer |
| Family consensus | 53.631–67.488 s per sample; about 111 s five-core critical path | 134 s | 5.553 s | about 20.0× vs Julia algorithm critical path; 24.1× vs workflow wall |
| Contamination | 11.741 s | 67 s | 2.351 s | 4.99× vs Julia algorithm |
| Downstream filters + MSA | 83.008–119.234 s per sample, including downstream tree work | 203 s | 13.364 s | see combined row |
| FastTree | not isolated by Julia | included above | 38.904 s | six trees, bounded parallel |
| Downstream + trees | included above | 203 s | 52.268 s | 3.88× vs workflow wall |
| Static reports + index | not isolated internally | 18 s + 21 s | not applicable | web figures render interactively from the result bundle |
| Analysis through reports/index | — | about 590 s | 66.852 s | 8.83× end-to-end |

The legacy workflow's external process continued to 605.034 seconds only because the final packaging rule ran and failed after the completed reports. That extra interval is excluded from the end-to-end comparison.

## Optimization impact

The freshly cloned GitHub baseline (`9338a12cc937c500294c2e90e283f362519eafb6`) took 759.62 seconds and spent 539.54 seconds in consensus on the same fixture. The optimized implementation takes 66.85 seconds overall and 5.553 seconds in consensus:

| Metric | GitHub baseline | Optimized | Speedup |
|---|---:|---:|---:|
| Consensus | 539.539 s | 5.553 s | 97.15× |
| Complete analysis | 759.623 s | 66.852 s | 11.36× |

The consensus reduction comes from exact sparse centroid arithmetic, collision-free rolling canonical seeds, adaptive seedless DP, fixed-point termination, alignment reuse, and allocation-free agreement/modal counting. Focused Julia substitution/indel parity remains exact, and all 3,043 sequences and minimum-agreement values are unchanged between the first corrected-semantics optimized run and the final fastest build.

## CPU scaling and determinism

The same CLI run with `--workers 1` took 192.204 seconds. Its consensus FASTA and postprocessing CSV have the same SHA-256 values as the five-worker outputs.

| Stage | 1 worker | 5 workers | Scaling |
|---|---:|---:|---:|
| Preprocessing | 22.351 s | 4.051 s | 5.52× |
| UMI grouping | 5.320 s | 2.496 s | 2.13× |
| Consensus | 22.870 s | 5.553 s | 4.12× |
| Contamination | 2.466 s | 2.351 s | 1.05× |
| Postprocessing | 25.046 s | 13.364 s | 1.87× |
| Six trees | 114.017 s | 38.904 s | 2.93× |
| Complete analysis | 192.204 s | 66.852 s | 2.88× |

Default worker count is the runtime's logical CPU count; browser and CLI controls can lower it. Per-sample MSA/tree pools are capped at the number of actual sample jobs to avoid allocating idle WASM runtimes on high-core machines.

## Output parity on this run

| Check | Julia | webPORPID | Agreement |
|---|---:|---:|---|
| Consensus keys | 3,043 | 3,043 | all shared |
| Exact consensus sequence | — | — | 3,040 / 3,043 (99.90%) |
| Exact minimum agreement | — | — | 2,993 / 3,043 (98.36%) |
| Minimum agreement within 0.01 | — | — | 3,040 / 3,043 (99.90%) |
| Contamination passed | 3,043 | 3,043 | exact |
| Postprocessing passed, six samples | 496, 492, 467, 494, 528, 466 | 496, 492, 467, 494, 528, 466 | exact |
| Functional passed, three configured samples | 459, 374, 347 | 464, 373, 348 | aggregate +5 (`0.34%`) across 1,455 evaluated |
