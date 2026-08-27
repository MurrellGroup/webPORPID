using PORPIDpipeline, PORPID, RobustAmpliconDenoising, BioSequences, FASTX

input_path = abspath(ARGS[1])
sample_name = "sample_1"
sample_config = Dict(
    "cDNA_primer" => "CCGCTacgtaaNNNNNNNNGTCA",
    "sec_str_primer" => "TAGG",
)

mktempdir() do directory
    demux = joinpath(directory, "demux"); mkpath(demux)
    counts = redirect_stdout(devnull) do
        chunked_filter_apply(input_path, demux, chunked_quality_demux; chunk_size=10000, f_kwargs=[
            :demux_dir => demux, :samples => Dict(sample_name => sample_config), :verbose => false,
            :error_rate => 0.05, :min_length => 20, :max_length => 300,
            :primer_tol => 0, :primer_window => 100, :primer_chop => 0,
            :label_prefix => "seq", :error_out => true,
        ])
    end
    println("counts\t", join(counts, '\t'))

    demux_path = joinpath(demux, sample_name * ".fastq.gz")
    id = uppercase(match(r"[a-z]+", sample_config["cDNA_primer"]).match)
    umi_range = findfirst(r"N+", sample_config["cDNA_primer"])
    suffix = replace(sample_config["cDNA_primer"][first(umi_range):end], "N" => "n") * "*"
    template_text = id * suffix
    cfg = Configuration(); cfg.files = [demux_path]; cfg.filetype = fastq
    cfg.start_inclusive = 0; cfg.end_inclusive = 6 + length(suffix) + 2; cfg.try_reverse_complement = false
    push!(cfg.templates, Template(sample_name, template_text))
    grouped = Dict{String,Vector{Any}}()
    callback(source_file_name, template, tag, output_sequence, score) = push!(get!(grouped, tag, Any[]), output_sequence)
    redirect_stdout(devnull) do
        extract_tags_from_file(demux_path, cfg, callback; print_every=10000, print_callback=(_)->nothing)
    end
    haskey(grouped, "REJECTS") && delete!(grouped, "REJECTS")
    tag_counts = Dict(tag => Int32(length(records)) for (tag, records) in grouped)
    tag_to_index, index_to_tag = tag_index_mapping(keys(tag_counts))
    likelihoods = prob_observed_tags_given_reals(tag_to_index, PORPID.PacBioErrorModel(0.005), 1)
    assignments = LDA(likelihoods, index_counts(tag_counts, tag_to_index))
    for observed in eachindex(assignments)
        tag = index_to_tag[observed]; println("family\t", tag, '\t', tag_counts[tag], '\t', assignments[observed][2])
    end

    tag = first(keys(grouped))
    reads = [String(FASTQ.sequence(record)) for record in grouped[tag]]
    draft = RobustAmpliconDenoising.consensus_seq(reads)
    final = RobustAmpliconDenoising.refine_ref(RobustAmpliconDenoising.refine_ref(draft, reads), reads)
    _, _, matches, _ = RobustAmpliconDenoising.get_matches(final, reads, 0)
    minimum_agreement = round(minimum(matches[2:end-1]); digits=2)
    sid_start = first(findfirst(r"[a-z]+", sample_config["cDNA_primer"]))
    trim_primer = uppercase(sample_config["cDNA_primer"][sid_start:end])
    sequence = PORPIDpipeline.reverse_complement(RobustAmpliconDenoising.primer_trim(final, trim_primer))
    println("consensus\t", tag, '\t', minimum_agreement, '\t', sequence)
end
