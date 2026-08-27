using PORPID

counts = Dict{String,Int32}(
    "AACCGGT" => 3, "AACCGGTA" => 4, "AACCGGTT" => 100, "AACCGTT" => 2,
    "AAACCGGTT" => 2, "TTGGCCAA" => 20, "TTGGCCAT" => 1,
)
tag_to_index, index_to_tag = tag_index_mapping(keys(counts))
likelihoods = prob_observed_tags_given_reals(tag_to_index, PORPID.PacBioErrorModel(0.005), 1)
assignments = LDA(likelihoods, index_counts(counts, tag_to_index))
for observed in 1:length(assignments)
    tag = index_to_tag[observed]
    parent = index_to_tag[assignments[observed][1]]
    probability = assignments[observed][2]
    disposition = probability < 0.995 ? 3 : (length(tag) != 8 ? 4 : 0)
    println(tag, '\t', parent, '\t', counts[tag], '\t', probability, '\t', disposition)
end
