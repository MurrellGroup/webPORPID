# webPORPID 0.3.8 performance audit

This audit covers the browser/CLI pipeline from streamed FASTQ input through reporting. Measurements below use only the deterministic synthetic generator in `scripts/benchmark-hotpaths.mjs`; no uploaded sequence data is read.

## Implemented in 0.3.8

| Hot path | Previous cost | 0.3.8 change | Scientific behavior |
|---|---|---|---|
| Contamination DP-means | Up to 30 complete assignment/mean passes per sample | Hard cap of three passes, with early convergence | Deliberate requested approximation. Center assignments after pass three are used as-is. Threshold inequalities and exact distance kernels are unchanged. |
| Donor-aware contamination | Same-sample was the only self category | Precomputed donor/sample self-group key is used by both indexed and fallback nearest-neighbor paths | Samples with the same nonblank `donor_ID` cannot be non-self candidates for one another. Samples without an ID retain the previous per-sample behavior. |
| Consensus partition grouping | Ordered `std::map` insertion for every selected read: O(reads log families) | Reserved hash table followed by the existing deterministic within-family read sort | Consensus records are sorted by the caller, so output content and final order are preserved. |
| Post-processing dispatch | Every worker-scoped sample scanned the complete consensus array | One linear pre-bucketing pass by sample index/name | Exact same records per sample. |
| Pipeline/sample summaries | Repeated whole-array `filter` calls per sample and per statistic | Single-pass typed counters and per-sample buckets | Exact counts and percentages. |
| Alignment modal sequence | A new `Map` and sort at every alignment column | Reused 128-symbol typed histogram plus first-observed tie order | Exact modal sequence, including the existing first-observed tie rule. |
| APOBEC grid | Repeated `Math.log` of every transition probability for every retained sequence | Log transition matrices are cached once with the existing cached grid | Exact stored APOBEC values in the differential fixture. |
| Consensus-threshold replay | A changed downstream cutoff required rerunning the FASTQ pipeline | Reconstructs the decision plot from compact stored families/consensuses and recomputes only previously completed downstream stages | UMI classifications, consensus calls, and contamination calls are preserved; display thinning is never used for classification. |
| Interactive plot preparation | Repeated per-sample scans would scale as samples × families | One linear bucketing pass builds the consensus checkpoint inputs | Exact frequency tables and decisions are unchanged. |

The existing 0.3.7 optimizations remain in place: fused IUPAC resolution/six-mer counting, sparse vectors, exact inverted posting indexes, norm lower bounds, elapsed-time cooperative yields, selected-read spool compaction, bounded partition concurrency for external scratch, seeded/banded alignment with full-matrix fallback, fixed-point consensus exit, and reuse of fixed-point alignments for agreement.

## Reproducible synthetic timings

Environment: Node v24.19.0, Linux x64, one process. Timings are wall-clock and should be treated as directional rather than hardware-independent throughput claims.

| Fixture | 0.3.7 | 0.3.8 | Speedup |
|---|---:|---:|---:|
| Contamination, 6,000 × 600-nt consensuses, 6 samples, 80 panel rows | 3.963 s | 3.174 s | 1.25× |
| Downstream filtering/annotation, 250 × 600-nt consensuses | 0.222 s | 0.158 s | 1.41× |
| Across-sample overview, 100,000 families across 100 samples | 0.116 s | 0.011 s | 10.5× |

The downstream non-APOBEC record hash was identical (`5e203a…ba9`) and a direct APOBEC differential was also exact. The contamination hash intentionally differs because the requested three-pass cap can stop before the former 30-pass fixed point; the number of calls in this fixture happened to remain 5,832.

Run the current benchmark with:

```bash
node scripts/benchmark-hotpaths.mjs
```

The script accepts module URL environment variables so a prior source tree can be used for a local differential without copying data.

## Remaining opportunities, not silently changed

1. **Parallel contamination classification.** The exact database index is currently built and queried in one JavaScript worker. A shared immutable index plus query shards could use more cores, but transferring or reconstructing a large posting index may erase the gain. This should be implemented with a size-based planner and benchmarked in browsers, not assumed from Node.
2. **UMI two-edit likelihood construction.** This remains a major CPU/memory candidate when the number of distinct observed tags is very high. A packed two-bit eight-mer specialization could remove many temporary strings, while a general indel-aware path remains necessary for abnormal-length UMIs. It needs a large differential corpus before replacing the current behavior.
3. **Panel/profile and MSA work.** Joint alignment is inherently expensive. The current scalable MSA already bounds memory with anchor batches. More parallel batches are possible, but merge order and equally scoring gap placement must remain deterministic.
4. **APOBEC opt-in/defer control.** Annotation can dominate runs with very many retained consensuses. A future optional/deferred annotation stage would be behavior-preserving and likely more valuable than approximating the posterior grid. It is recorded as an opportunity rather than changing workflow semantics in this release.
5. **Worker-local partition decoding.** Large partitions are decoded once for UMI counting and again for consensus by design to keep peak memory bounded. Retaining decoded reads would save CPU but violate the large-dataset memory objective. External scratch plus more partitions remains the preferred trade-off.
