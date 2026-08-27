using PORPIDpipeline

cases = [
    (["ACGTACGT", "ACGTTCGT"], ["ACGTACGT", "ACGTTCGT"]),
    (["AACCGGTT", "AACC-GTT"], ["TTAACCGGTTAA", "TTAACC-GTTAA"]),
    (["ACGTACGT", "ACGTACGT"], ["GGACGTTACGTCC", "GGACG-TACGTCC"]),
]

for (case_index, (panel, sample)) in enumerate(cases)
    extracted, scores = PORPIDpipeline.extract_and_score_misalignments(sample, seqs2profile(panel))
    for row in eachindex(extracted)
        println(case_index, '\t', row, '\t', scores[row], '\t', extracted[row])
    end
end
