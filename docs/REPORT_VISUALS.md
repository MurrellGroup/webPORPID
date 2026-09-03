# Julia report-to-web figure mapping

The supplied Julia run generated five PNG figures and one tree SVG for each of six samples. webPORPID renders the same report concepts from the portable result bundle rather than embedding static images in generated HTML. The goal is semantic and visual correspondence with responsive, inspectable graphics—not pixel-identical Matplotlib output.

| Julia report asset | webPORPID figure | Matched encoding |
|---|---|---|
| `*_qc_bins.png` | “UMI family size × UMI length” | x = log2 family size; y = UMI length with deterministic jitter; PORPID decision colors; orange non-outlier quantile and red artefact cutoff; Julia-style 7/8/9 UMI-length focus |
| `*_artefacts.png` | “Artefact-cutoff decision” | horizontal family-size strip; deterministic jitter; 5%, 15%, …, 95% guide lines; red artefact and orange quantile lines; decision colors and threshold summary |
| `*_minags.png` | “Low-agreement positions” | one longest modal-run minimum site per eligible family; x = 3′-relative sequence position; y = agreement; point size = modal run length; C/G/T/A colors; purple threshold |
| `*.fasta.mds.png` | “P(APOBEC\|mutations) · classical MDS” | classical MDS of normalized pairwise non-gap Hamming distance; point area follows family size; jet color follows stored APOBEC posterior; all points in this validation run are shown |
| `*_di_nuc_freq.png` | “UMI dinucleotide frequencies” | 4×4 family-unweighted and read-weighted matrices; A/C/G/T columns, reversed T/G/C/A rows, annotated values, fixed 0–0.20 YlOrBr-like scale and color bars |
| `*.tre.svg` | Alignment workbench phylogram | exact Swig Newick parsing/layout and branch normalization; rectilinear branches, scale bar, adjustable width, topology/branch-length modes, ladderization, tree SVG/Newick export |

The Julia report tables map to webPORPID's Overview, Families, Sequences, Contamination, and Log tabs. Counters and per-record decisions remain exportable as CSV or FASTA.

## Tree and alignment behavior

The previous tree component compressed ordinary fractional branch lengths against a denominator of one. The Swig layout normalizes against the observed maximum root-to-tip distance, so trees now fill the chosen width while retaining relative branch lengths. Alignment cells are drawn in a virtualized canvas and share vertical scroll/selection with the SVG tree, allowing large alignments without one DOM element per residue.

The Protein button always translates the active nucleotide alignment and never requires a functional-filter pass. Complete gap codons render as `-`; mixed-gap, ambiguous, and incomplete codons render as `X`. The user can choose nucleotide-column frame 1, 2, or 3.

The export-all static-figure toggle produces a publication-style SVG in which the phylogram and modal highlighter each occupy exactly half of a 1600-unit canvas. Tree leaves and alignment rows share one y coordinate; the alignment uses a top coordinate ruler, horizontal row guides, white modal matches, A/G/T/C mutation colors, grey gaps, a bottom nucleotide legend, and a branch-length scale. Abundance circles retain exact linear area (`25 px² × represented UMI-family count`) with no floor or cap. Functional-reference rows use a separate outlined marker and therefore do not imply a sampled-family abundance.

## Alivibe round trip

“Open in Alivibe” launches the bundled, same-origin Swig Alivibe build. webPORPID verifies the bridge revision, loads and re-reads the exact nucleotide FASTA, installs the local MSA runner, and adds a “Return alignment to webPORPID” control. On return it rejects stale sample/alignment sessions, verifies identifiers/order/residues, accepts gap editing plus row/base deletion, reruns FastTree, and stores the corrected FASTA, frame, fingerprints, edit provenance, and refreshed tree inside the next `.webporpid` save. A corrected FASTA can also be imported directly, and “Restore pipeline alignment” removes the edit.
