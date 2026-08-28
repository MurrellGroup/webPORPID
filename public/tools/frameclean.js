/*!
 * frameclean.js — Reference-anchored frame assignment + frame-correction post-processor for nucleotide MSAs
 *
 * This library assigns a reading frame to alignment coordinates and (optionally) produces a “frame-corrected”
 * alignment where backbone columns land on consistent codon phases.
 *
 * Key design (new implementation):
 *   - If a trusted reference is provided, then:
 *       * Columns where the reference has a base are ANCHORS: their phase is fixed by the reference ORF.
 *       * Columns where the reference has a gap are “unanchored”: the model decides whether each column
 *         should be counted as backbone (advances the frame coordinate; typically reference deletions),
 *         or treated as insertion/noise (does not advance the coordinate).
 *       * Between consecutive anchor columns, the number of “counted backbone” columns must be ≡ 0 (mod 3)
 *         in hard mode, or may be ≡ 1/2 with a penalty in soft mode.
 *
 * API compatibility:
 *   - FrameClean.inferFrame(alignment, options) -> { frameVec, backboneMask, meta }
 *   - FrameClean.buildCorrectedAlignment(alignment, inferred, options) -> { alignment, frameVec, meta }
 *   - FrameClean.models.STANDARD_CODE, FrameClean.models.translateCodon(codon, code)
 *
 * Alignment format:
 *   alignment = [{ id: "s0", seq: "ACG--T..." }, ...] (all seqs same length).
 *
 * Returned vectors:
 *   - backboneMask[j] = true if column j is backbone (advances coordinate), else false (insertion/noise).
 *   - frameVec[j] = 0/1/2 for backbone columns; null for insertion/noise columns.
 *
 * Notes:
 *   - With a reference, this returns a REFERENCE-ANCHORED frame coordinate system.
 *   - Without a reference, it falls back to an occupancy-only HMM that does not know “true ancestry”.
 */

(function (root, factory) {
    if (typeof module === "object" && typeof module.exports === "object") {
      module.exports = factory();
    } else {
      root.FrameClean = factory();
    }
  })(typeof self !== "undefined" ? self : this, function () {
    "use strict";

    /* ---------------------------
     * Utilities
     * ------------------------- */

    function assert(cond, msg) {
      if (!cond) throw new Error(msg || "Assertion failed");
    }

    function log(x) {
      return Math.log(x);
    }

    function clamp(x, lo, hi) {
      return Math.max(lo, Math.min(hi, x));
    }

    function safeProb(p) {
      // Avoid 0/1 exactly for logs.
      return clamp(p, 1e-9, 1 - 1e-9);
    }

    function isGap(c) {
      return c === "-";
    }

    function normalizeAlignment(alignment) {
      assert(Array.isArray(alignment) && alignment.length > 0, "alignment must be a non-empty array");
      const L = alignment[0].seq.length;
      for (const s of alignment) {
        assert(typeof s.id === "string", "each sequence needs an id");
        assert(typeof s.seq === "string", "each sequence needs a seq string");
        assert(s.seq.length === L, "all sequences must have the same aligned length");
      }
      // Uppercase; keep '-' as gaps.
      return alignment.map((s) => ({
        id: s.id,
        seq: s.seq.toUpperCase(),
      }));
    }

    function getColumnChars(aln, j) {
      const col = new Array(aln.length);
      for (let i = 0; i < aln.length; i++) col[i] = aln[i].seq[j];
      return col;
    }

    function columnNonGapCount(aln, j) {
      let k = 0;
      for (let i = 0; i < aln.length; i++) if (!isGap(aln[i].seq[j])) k++;
      return k;
    }

    function argmax3(a0, a1, a2) {
      if (a0 >= a1 && a0 >= a2) return 0;
      if (a1 >= a0 && a1 >= a2) return 1;
      return 2;
    }

    /* ---------------------------
     * Genetic code helpers
     * ------------------------- */

    // Standard genetic code (DNA codons -> AA one-letter; '*' stop)
    const STANDARD_CODE = {
      TTT: "F", TTC: "F", TTA: "L", TTG: "L",
      TCT: "S", TCC: "S", TCA: "S", TCG: "S",
      TAT: "Y", TAC: "Y", TAA: "*", TAG: "*",
      TGT: "C", TGC: "C", TGA: "*", TGG: "W",

      CTT: "L", CTC: "L", CTA: "L", CTG: "L",
      CCT: "P", CCC: "P", CCA: "P", CCG: "P",
      CAT: "H", CAC: "H", CAA: "Q", CAG: "Q",
      CGT: "R", CGC: "R", CGA: "R", CGG: "R",

      ATT: "I", ATC: "I", ATA: "I", ATG: "M",
      ACT: "T", ACC: "T", ACA: "T", ACG: "T",
      AAT: "N", AAC: "N", AAA: "K", AAG: "K",
      AGT: "S", AGC: "S", AGA: "R", AGG: "R",

      GTT: "V", GTC: "V", GTA: "V", GTG: "V",
      GCT: "A", GCC: "A", GCA: "A", GCG: "A",
      GAT: "D", GAC: "D", GAA: "E", GAG: "E",
      GGT: "G", GGC: "G", GGA: "G", GGG: "G",
    };

    function translateCodon(codon, code) {
      const c = (codon || "").toUpperCase();
      if (c.length !== 3) return "X";
      if (c.indexOf("-") !== -1) return "X";
      if (c.indexOf("N") !== -1) return "X";
      const aa = (code || STANDARD_CODE)[c];
      return aa ? aa : "X";
    }

    /* ---------------------------
     * Scoring model (occupancy emissions + simple priors)
     * ------------------------- */

    function makeColumnEmissionFn(aln, opts) {
      const N = aln.length;
      const pB = safeProb((opts.columnModelOptions && opts.columnModelOptions.pBackbone) != null ? opts.columnModelOptions.pBackbone : 0.9);
      const pI = safeProb((opts.columnModelOptions && opts.columnModelOptions.pInsertion) != null ? opts.columnModelOptions.pInsertion : 0.05);

      const logpB = log(pB), logqB = log(1 - pB);
      const logpI = log(pI), logqI = log(1 - pI);

      // Optional user hook: add additional per-column score
      const extra = (opts && typeof opts.columnScoreFn === "function") ? opts.columnScoreFn : null;

      // Return emission(j, state) where state: 0 insertion, 1 backbone
      return function emission(j, state) {
        const k = columnNonGapCount(aln, j);
        const base = (state === 1)
          ? (k * logpB + (N - k) * logqB)
          : (k * logpI + (N - k) * logqI);
        if (!extra) return base;
        return base + (extra({ colIndex: j, state, nonGapCount: k, N }) || 0);
      };
    }

    function getPriorParams(opts) {
      const prior = (opts && opts.priorModelOptions) ? opts.priorModelOptions : {};
      return {
        // In this implementation (with reference):
        //   These penalize selecting BACKBONE columns inside reference-gap segments (i.e. explaining ref gaps as ref deletions).
        //   This is a *policy choice* for API compatibility; emissions still dominate.
        backboneStartPenalty: (prior.insertionStartPenalty != null) ? prior.insertionStartPenalty : -2.0,
        backboneExtendPenalty: (prior.insertionExtendPenalty != null) ? prior.insertionExtendPenalty : -0.2,

        // Residue penalties for soft mode when #backbone_between_anchors mod3 != 0
        // Default: strongly discourage residues 1/2.
        residuePenalty: Array.isArray(prior.insertionEndResidPenalty) && prior.insertionEndResidPenalty.length === 3
          ? prior.insertionEndResidPenalty.slice(0, 3)
          : [0.0, -6.0, -6.0],
      };
    }

    /* ---------------------------
     * Reference-anchored inference
     * ------------------------- */

    function computeReferencePhases(aln, refIndex, frameOffset) {
      const ref = aln[refIndex].seq;
      const L = ref.length;

      const isAnchor = new Array(L);
      const refUngappedPos = new Array(L); // -1 for gaps, else 0-based ungapped index
      let u = -1;
      for (let j = 0; j < L; j++) {
        const c = ref[j];
        if (!isGap(c)) {
          u++;
          isAnchor[j] = true;
          refUngappedPos[j] = u;
        } else {
          isAnchor[j] = false;
          refUngappedPos[j] = -1;
        }
      }

      const anchorPhase = new Array(L).fill(null);
      for (let j = 0; j < L; j++) {
        if (!isAnchor[j]) continue;
        const idx = refUngappedPos[j]; // 0-based
        anchorPhase[j] = (idx + (frameOffset | 0)) % 3;
        if (anchorPhase[j] < 0) anchorPhase[j] += 3;
      }

      return { isAnchor, anchorPhase };
    }

    /**
     * Solve one reference-gap segment between anchors using DP.
     *
     * Segment columns are indices segCols[] such that reference has gaps there.
     * We decide for each column whether it is:
     *   - insertion/noise (state 0): does NOT advance coordinate; frameVec null
     *   - backbone (state 1): advances coordinate; frameVec assigned; contributes to residue count
     *
     * Hard mode: total backbone count in segment must be ≡ 0 (mod 3).
     * Soft mode: any residue allowed but penalized.
     */
    function solveGapSegmentDP(segCols, emission, prior, mode) {
      const n = segCols.length;
      if (n === 0) {
        return { backbone: [], endResid: 0, score: 0.0 };
      }

      // dp[i][res][s] for i in 0..n, res in 0..2, s in {0,1}
      // Use flat arrays for speed.
      const NEG = -1e300;

      const dp = new Float64Array((n + 1) * 3 * 2);
      const bt = new Int32Array((n + 1) * 3 * 2); // backtrace packed: prevRes*2 + prevS, plus 6 for "unset"
      for (let t = 0; t < dp.length; t++) { dp[t] = NEG; bt[t] = -1; }

      function idx(i, res, s) { return ((i * 3 + res) * 2 + s); }

      // Start: before segment, residue 0, previous state insertion (0)
      dp[idx(0, 0, 0)] = 0.0;
      bt[idx(0, 0, 0)] = -2;

      for (let i = 0; i < n; i++) {
        const colJ = segCols[i];

        for (let res = 0; res < 3; res++) {
          for (let s = 0; s < 2; s++) {
            const cur = dp[idx(i, res, s)];
            if (cur <= NEG / 2) continue;

            // Next state s' = 0 (insertion)
            {
              const res2 = res;
              const s2 = 0;
              const sc = cur + emission(colJ, 0);
              const k = idx(i + 1, res2, s2);
              if (sc > dp[k]) {
                dp[k] = sc;
                bt[k] = res * 2 + s;
              }
            }

            // Next state s' = 1 (backbone)
            {
              const res2 = (res + 1) % 3;
              const s2 = 1;
              let trans = 0.0;
              if (s === 0) trans += prior.backboneStartPenalty;
              else trans += prior.backboneExtendPenalty;

              const sc = cur + trans + emission(colJ, 1);
              const k = idx(i + 1, res2, s2);
              if (sc > dp[k]) {
                dp[k] = sc;
                bt[k] = res * 2 + s;
              }
            }
          }
        }
      }

      // Choose best ending (any prev state) with residue constraint/penalty
      const resPen = prior.residuePenalty;
      let bestRes = 0, bestS = 0, bestScore = NEG;

      for (let res = 0; res < 3; res++) {
        const pen = (mode === "soft") ? resPen[res] : (res === 0 ? 0.0 : NEG);
        for (let s = 0; s < 2; s++) {
          const sc = dp[idx(n, res, s)] + pen;
          if (sc > bestScore) {
            bestScore = sc;
            bestRes = res;
            bestS = s;
          }
        }
      }

      // Backtrace decisions
      const backbone = new Array(n).fill(false);
      let i = n, res = bestRes, s = bestS;

      while (i > 0) {
        const packed = bt[idx(i, res, s)];
        if (packed < 0) break;
        // Decision at i-1 is state s (current state at dp cell corresponds to state at column i-1)
        backbone[i - 1] = (s === 1);

        const prevRes = Math.floor(packed / 2);
        const prevS = packed % 2;

        // If current s is backbone, residue came from prevRes -> res = (prevRes+1)%3, so prevRes is consistent.
        // If current s is insertion, residue unchanged.
        res = prevRes;
        s = prevS;
        i--;
      }

      return { backbone, endResid: bestRes, score: bestScore };
    }

    function inferFrameWithReference(aln, opts) {
      const ref = opts.reference || {};
      const refIndex = clamp(ref.index | 0, 0, aln.length - 1);
      const mode = (ref.mode === "soft") ? "soft" : "hard";
      const frameOffset = ref.frameOffset | 0;

      const { isAnchor, anchorPhase } = computeReferencePhases(aln, refIndex, frameOffset);
      const L = aln[0].seq.length;

      const emission = makeColumnEmissionFn(aln, opts);
      const prior = getPriorParams(opts);

      const backboneMask = new Array(L).fill(false);
      const frameVec = new Array(L).fill(null);

      // Identify anchors
      const anchors = [];
      for (let j = 0; j < L; j++) {
        if (isAnchor[j]) anchors.push(j);
      }

      // If reference has no anchors, fall back.
      if (anchors.length === 0) {
        return inferFrameNoReference(aln, opts);
      }

      // Anchors: always backbone, phase fixed.
      for (const j of anchors) {
        backboneMask[j] = true;
        frameVec[j] = anchorPhase[j];
      }

      // Solve each gap segment between consecutive anchors independently
      for (let t = 0; t + 1 < anchors.length; t++) {
        const left = anchors[t];
        const right = anchors[t + 1];

        // Segment columns are (left, right) exclusive where ref has gaps (should all be gaps, but guard).
        const segCols = [];
        for (let j = left + 1; j <= right - 1; j++) {
          if (!isAnchor[j]) segCols.push(j);
          else {
            // If there is an unexpected anchor inside, we break segments at anchors anyway
            // (but this shouldn't happen because anchors list is consecutive).
          }
        }

        const sol = solveGapSegmentDP(segCols, emission, prior, mode);
        // Assign backbone vs insertion for segCols
        for (let i = 0; i < segCols.length; i++) {
          const j = segCols[i];
          if (sol.backbone[i]) backboneMask[j] = true;
        }

        // Assign phases in this region:
        // Let runningPhase be the phase that the NEXT backbone column after 'left' should take.
        // Since 'right' is the next reference base, its phase is fixed to (phase(left)+1) mod 3.
        let runningPhase = (anchorPhase[left] + 1) % 3;

        for (let i = 0; i < segCols.length; i++) {
          const j = segCols[i];
          if (!backboneMask[j]) continue;
          frameVec[j] = runningPhase;
          runningPhase = (runningPhase + 1) % 3;
        }

        // Ensure right anchor stays consistent (in hard mode it must; in soft mode, it is fixed anyway).
        // We don't “propagate” runningPhase into anchors; anchors are fixed.
        // But it’s useful to detect gross inconsistencies for debugging:
        // (We do not throw; we just note it in meta.)
      }

      // Leading/trailing reference-gap columns (before first anchor, after last): treat as insertion by default.
      // If you want something more elaborate, use columnScoreFn + your own policy.
      const firstA = anchors[0];
      for (let j = 0; j < firstA; j++) {
        if (!isAnchor[j]) { backboneMask[j] = false; frameVec[j] = null; }
      }
      const lastA = anchors[anchors.length - 1];
      for (let j = lastA + 1; j < L; j++) {
        if (!isAnchor[j]) { backboneMask[j] = false; frameVec[j] = null; }
      }

      // Meta
      const meta = {
        mode: "reference",
        reference: { index: refIndex, mode, frameOffset },
        anchorsCount: anchors.length,
      };

      return { frameVec, backboneMask, meta };
    }

    /* ---------------------------
     * No-reference inference (simple fallback)
     *
     * This is intentionally simple: an HMM over all columns with an occupancy emission and an insertion prior,
     * plus a global offset chosen by codon coherence (weak).
     * ------------------------- */

    function chooseOffsetByStopRate(aln, backboneMask, opts) {
      // Optional codon scoring to pick a delta such that translated codons (from backbone columns) minimize stops.
      const codOpts = opts.codonModelOptions || {};
      const stopPenalty = (codOpts.stopCodonPenalty != null) ? codOpts.stopCodonPenalty : 8.0;

      // Build backbone column indices
      const cols = [];
      for (let j = 0; j < backboneMask.length; j++) if (backboneMask[j]) cols.push(j);
      if (cols.length < 3) return 0;

      // Evaluate each delta in {0,1,2} by counting stop codons across sequences on backbone triplets.
      // This is a crude heuristic, but stable.
      function scoreDelta(delta) {
        let score = 0.0;
        // We need triplets in backbone coordinate: phase = (k + delta) % 3
        // So codon boundaries occur where (k+delta)%3==0.
        // We'll scan backbone index k.
        const M = cols.length;
        for (let i = 0; i + 2 < M; i++) {
          const k = i; // backbone index
          if (((k + delta) % 3) !== 0) continue;
          const j0 = cols[i], j1 = cols[i + 1], j2 = cols[i + 2];
          for (let s = 0; s < aln.length; s++) {
            const c0 = aln[s].seq[j0], c1 = aln[s].seq[j1], c2 = aln[s].seq[j2];
            if (isGap(c0) || isGap(c1) || isGap(c2)) continue;
            const codon = c0 + c1 + c2;
            if (translateCodon(codon, STANDARD_CODE) === "*") score -= stopPenalty;
          }
        }
        return score;
      }

      const s0 = scoreDelta(0), s1 = scoreDelta(1), s2 = scoreDelta(2);
      return argmax3(s0, s1, s2);
    }

    function inferBackboneMaskNoRef(aln, opts) {
      // 2-state HMM over all columns: state 1 backbone, state 0 insertion
      // Uses emission + insertion run priors from priorModelOptions (old naming).
      const L = aln[0].seq.length;
      const emission = makeColumnEmissionFn(aln, opts);
      const prior = getPriorParams(opts);

      // Interpret the same prior as before in a no-ref setting:
      // insertionStartPenalty / insertionExtendPenalty are penalties for being in insertion state (0).
      // To preserve historical behavior loosely, we flip signs by using them as insertion rewards/penalties:
      // Here, we simply penalize transitions into insertion and staying in insertion using their negatives.
      const insStart = (opts.priorModelOptions && opts.priorModelOptions.insertionStartPenalty != null)
        ? opts.priorModelOptions.insertionStartPenalty
        : -2.0;
      const insExtend = (opts.priorModelOptions && opts.priorModelOptions.insertionExtendPenalty != null)
        ? opts.priorModelOptions.insertionExtendPenalty
        : -0.2;

      // dp[j][s]
      const NEG = -1e300;
      let dp0 = 0.0, dp1 = 0.0;
      const bt = new Int8Array(L * 2); // store prev state for each state at j
      for (let j = 0; j < L; j++) {
        const e0 = emission(j, 0);
        const e1 = emission(j, 1);

        // Transition penalties:
        // entering insertion (from backbone): insStart
        // staying insertion: insExtend
        // backbone has no penalty (simple)
        const ndp0_from0 = dp0 + insExtend + e0;
        const ndp0_from1 = dp1 + insStart + e0;
        const ndp1_from0 = dp0 + e1;
        const ndp1_from1 = dp1 + e1;

        let ndp0, ndp1, prev0, prev1;
        if (ndp0_from0 >= ndp0_from1) { ndp0 = ndp0_from0; prev0 = 0; } else { ndp0 = ndp0_from1; prev0 = 1; }
        if (ndp1_from0 >= ndp1_from1) { ndp1 = ndp1_from0; prev1 = 0; } else { ndp1 = ndp1_from1; prev1 = 1; }

        bt[j * 2 + 0] = prev0;
        bt[j * 2 + 1] = prev1;

        dp0 = ndp0;
        dp1 = ndp1;
      }

      // Backtrace
      const backboneMask = new Array(L).fill(false);
      let s = (dp1 >= dp0) ? 1 : 0;
      for (let j = L - 1; j >= 0; j--) {
        backboneMask[j] = (s === 1);
        s = bt[j * 2 + s];
      }
      return backboneMask;
    }

    function inferFrameNoReference(aln, opts) {
      const L = aln[0].seq.length;
      const backboneMask = inferBackboneMaskNoRef(aln, opts);
      const delta = chooseOffsetByStopRate(aln, backboneMask, opts);

      // Assign frameVec by counting backbone columns
      const frameVec = new Array(L).fill(null);
      let k = 0;
      for (let j = 0; j < L; j++) {
        if (!backboneMask[j]) continue;
        frameVec[j] = (k + delta) % 3;
        k++;
      }

      return {
        frameVec,
        backboneMask,
        meta: { mode: "noref", delta },
      };
    }

    /* ---------------------------
     * Public API: inferFrame
     * ------------------------- */

    /**
     * inferFrame(alignment, options)
     *
     * Options (backward-compatible shape; some fields are interpreted differently when reference is used):
     *
     *   options.reference (optional):
     *     { index: 0, mode: "hard"|"soft", frameOffset: 0 }
     *
     *     - index: which sequence in alignment is the trusted reference.
     *     - mode:
     *         "hard": between consecutive reference bases, the number of backbone columns in reference-gap runs
     *                 must be a multiple of 3 (i.e. ≡ 0 mod 3).
     *         "soft": residues 1/2 are allowed but penalized by priorModelOptions.insertionEndResidPenalty[r].
     *     - frameOffset: phase offset applied to reference ungapped index (0 means ref[0] is phase 0).
     *
     *   options.columnModelOptions:
     *     { pBackbone: 0.9, pInsertion: 0.05 }  // occupancy model
     *
     *   options.priorModelOptions:
     *     {
     *       insertionStartPenalty: -2.0,
     *       insertionExtendPenalty: -0.2,
     *       insertionEndResidPenalty: [0, -6, -6]
     *     }
     *
     *     With reference:
     *       - insertionStartPenalty / insertionExtendPenalty penalize selecting BACKBONE columns inside ref-gap segments
     *         (i.e. explaining ref gaps as ref deletions).
     *       - insertionEndResidPenalty is used only in soft mode as a penalty on (#backbone_in_gap mod 3).
     *
     *   options.codonModelOptions:
     *     { stopCodonPenalty: 8, aaConsensusWeight: 0.2 }
     *     (Used only in no-reference mode for a weak global offset heuristic.)
     *
     *   options.columnScoreFn (optional):
     *     ({colIndex, state, nonGapCount, N}) => number
     *     Add custom per-column additive score (state 0 insertion, 1 backbone).
     */
    function inferFrame(alignment, options) {
      const aln = normalizeAlignment(alignment);
      const opts = options || {};

      if (opts.reference && typeof opts.reference.index === "number") {
        return inferFrameWithReference(aln, opts);
      }
      return inferFrameNoReference(aln, opts);
    }

    /* ---------------------------
     * Frame-corrected alignment builder
     * ------------------------- */

    /**
     * buildCorrectedAlignment(alignment, inferred, options)
     *
     * Returns a new alignment whose column indices are consistent with the assigned frame:
     *   - Backbone columns are padded so that (outputColumnIndex % 3) == frameVec (0/1/2).
     *
     * Options:
     *   - mode: "pad" (default) or "remove"
     *       "pad": keep insertion columns, but insert padding all-gap columns when needed to keep backbone in-phase.
     *       "remove": drop all insertion/noise columns (backboneMask=false). This will change sequences.
     *
     *   - padToMultipleOf3: boolean (default false)
     *       If true, pads at the beginning so the first backbone column lands on its desired phase.
     *
     *   - keepInsertions: boolean (default true for "pad", ignored for "remove")
     *       If false, behaves like remove for insertion columns, but still pads backbone phases.
     */
    function buildCorrectedAlignment(alignment, inferred, options) {
      const aln = normalizeAlignment(alignment);
      const L = aln[0].seq.length;

      assert(inferred && Array.isArray(inferred.backboneMask) && Array.isArray(inferred.frameVec), "inferred must come from inferFrame()");
      assert(inferred.backboneMask.length === L && inferred.frameVec.length === L, "inferred vectors length mismatch");

      const opts = options || {};
      const mode = (opts.mode === "remove") ? "remove" : "pad";
      const padToMultipleOf3 = !!opts.padToMultipleOf3;
      const keepInsertions = (opts.keepInsertions != null) ? !!opts.keepInsertions : (mode === "pad");

      const N = aln.length;

      // Identify first backbone column (for optional leading padding)
      let firstBB = -1;
      for (let j = 0; j < L; j++) {
        if (inferred.backboneMask[j]) { firstBB = j; break; }
      }

      // Output builders
      const outChars = Array.from({ length: N }, () => []);
      const outFrame = [];

      let outIndex = 0;

      function pushPadColumn() {
        for (let i = 0; i < N; i++) outChars[i].push("-");
        outFrame.push(null);
        outIndex++;
      }

      function pushRealColumn(j, isBackbone) {
        for (let i = 0; i < N; i++) outChars[i].push(aln[i].seq[j]);
        outFrame.push(isBackbone ? (outIndex % 3) : null);
        outIndex++;
      }

      // Leading pad so first backbone phase matches (if requested)
      if (padToMultipleOf3 && firstBB >= 0) {
        const want = inferred.frameVec[firstBB];
        if (want != null) {
          while ((outIndex % 3) !== want) pushPadColumn();
        }
      }

      for (let j = 0; j < L; j++) {
        const isBB = !!inferred.backboneMask[j];

        if (!isBB) {
          if (mode === "remove") continue;
          if (!keepInsertions) continue;
          // Keep insertion column; it may shift outIndex, so we will re-pad before next backbone.
          pushRealColumn(j, false);
          continue;
        }

        const want = inferred.frameVec[j];
        assert(want === 0 || want === 1 || want === 2, "backbone column missing phase (frameVec)");

        while ((outIndex % 3) !== want) pushPadColumn();
        pushRealColumn(j, true);
      }

      const outAlignment = aln.map((s, i) => ({ id: s.id, seq: outChars[i].join("") }));

      return {
        alignment: outAlignment,
        frameVec: outFrame,
        meta: {
          mode,
          padToMultipleOf3,
          keptInsertions: keepInsertions && mode === "pad",
          inputLength: L,
          outputLength: outAlignment[0].seq.length,
        },
      };
    }

    /**
     * toFasta(alignment, wrap) -> string
     */
    function toFasta(aln, wrap = 80) {
      const lines = [];
      for (const s of aln) {
        lines.push(">" + s.id);
        const seq = s.seq;
        if (!wrap || wrap <= 0) lines.push(seq);
        else for (let i = 0; i < seq.length; i += wrap) lines.push(seq.slice(i, i + wrap));
      }
      return lines.join("\n") + "\n";
    }

    /**
     * cleanInFrame(alignment, inferred) -> { alignment, removedCols }
     *
     * Drops the minimum number of columns from each consecutive run of inferred insertion columns
     * so that the run length becomes a multiple of 3. Columns with lowest non-gap support are
     * preferred for removal.
     */
    function cleanInFrame(alignment, inferred) {
      const aln = normalizeAlignment(alignment);
      const L = aln[0].seq.length;
      assert(inferred && Array.isArray(inferred.backboneMask), "inferred.backboneMask is required");

      const isIns = inferred.backboneMask.map(b => !b);
      const toDrop = new Array(L).fill(false);

      let j = 0;
      while (j < L) {
        if (!isIns[j]) { j++; continue; }
        const start = j;
        while (j < L && isIns[j]) j++;
        const end = j;
        const runLen = end - start;
        const r = runLen % 3;
        if (r === 0) continue;

        const cols = [];
        for (let c = start; c < end; c++) {
          const sup = columnNonGapCount(aln, c);
          const edge = Math.min(c - start, (end - 1) - c);
          cols.push({ c, sup, edge });
        }
        cols.sort((a, b) => (a.sup !== b.sup) ? (a.sup - b.sup) : (a.edge - b.edge));

        for (let t = 0; t < r; t++) toDrop[cols[t].c] = true;
      }

      const cleaned = aln.map(s => ({ id: s.id, seq: "" }));
      for (let i = 0; i < aln.length; i++) {
        let out = "";
        for (let k = 0; k < L; k++) if (!toDrop[k]) out += aln[i].seq[k];
        cleaned[i].seq = out;
      }

      const removed = [];
      for (let k = 0; k < L; k++) if (toDrop[k]) removed.push(k);

      return { alignment: cleaned, removedCols: removed };
    }

    /* ---------------------------
     * Exports
     * ------------------------- */

    const FrameClean = {
      inferFrame,
      buildCorrectedAlignment,
      cleanInFrame,
      toFasta,
      models: {
        STANDARD_CODE,
        translateCodon,
      },
      // Internal hooks are intentionally minimal in this build.
      _internal: {
        inferFrameWithReference,
        inferFrameNoReference,
        computeReferencePhases,
        columnNonGapCount,
      },
    };

    return FrameClean;
  });
