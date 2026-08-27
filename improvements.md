# Potential improvements

These ideas may change behavior enough that they were not folded into the initial port.

- Replace the current pairwise consensus refinement with a banded partial-order graph consensus. This could greatly improve long, repetitive, high-indel families, but it would change tie paths and some local indel placement.
- Add adaptive consensus band sizing based on observed family indel density. This should be faster on clean families and safer on noisy families, but would need a new parity study.
- Replace the UMI posterior scan with a BK-tree or trie-indexed candidate search for very large UMI spaces. It should be faster when there are many observed UMIs, but the candidate boundary must be checked carefully.
- Add exact lower-bound pruning or a validated approximate-nearest-neighbor index to the contamination database search. The current sparse vectors reduce allocation and arithmetic, but nearest-cluster comparison is still quadratic in the number of sequences/clusters. Approximate indexing could be dramatically faster for very large runs and could also change a call near `dist_thresh`; even exact pruning must preserve Julia's first-argmin tie order.
- Replace shared-anchor large-MSA batching with a streaming profile or partial-order alignment that incrementally updates a run-wide graph. This could improve repetitive or highly divergent alignments and remove reliance on the first sequence as anchor, but equal-score indel placement and therefore downstream tree columns would change.
- Add optional GPU-backed tree/MSA paths for very large downstream alignments. The current standalone WebAssembly path is simpler and portable.
