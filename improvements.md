# Potential improvements

These ideas may change behavior enough that they were not folded into the initial port.

- Replace pairwise reference refinement with a partial-order graph or wavefront consensus. This could further improve long, repetitive, very-high-indel families, but it would change tie paths and some local indel placement.
- Learn the seedless alignment band from each family's observed indel distribution, or use sparse wavefront expansion. The implemented conservative length-delta-plus-square-root band is deterministic; a tighter family-specific boundary would be faster but needs a much larger parity study around repetitive and structurally variant reads.
- Use minimizer/syncmer chaining instead of exact unique 30-mer anchors. This should find anchors in noisier families and shorten more dynamic-programming intervals, but changes the anchor chain and therefore can change equal-score indel placement.
- Replace the UMI posterior scan with a BK-tree or trie-indexed candidate search for very large UMI spaces. It should be faster when there are many observed UMIs, but the candidate boundary must be checked carefully.
- Evaluate a validated approximate-nearest-neighbor index for contamination databases with an extreme number of clusters. The current implementation already uses exact reverse-triangle lower-bound pruning and preserves first-argmin tie order; approximate indexing could skip more comparisons, but could also change a call near `dist_thresh` and therefore needs a dedicated sensitivity study.
- Replace shared-anchor large-MSA batching with a streaming profile or partial-order alignment that incrementally updates a run-wide graph. This could improve repetitive or highly divergent alignments and remove reliance on the first sequence as anchor, but equal-score indel placement and therefore downstream tree columns would change.
- Add optional GPU-backed tree/MSA paths for very large downstream alignments. The current standalone WebAssembly path is simpler and portable.
- Evaluate a double-precision, pthread-enabled FastTree build or a divide-and-merge tree strategy for one very large uncollapsed sample. Current runs parallelize independent sample trees across cores, but one FastTree invocation is single-threaded; alternative parallel inference can change heuristics, topology, and branch lengths and therefore needs a dedicated downstream study.
