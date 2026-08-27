using SeededAlignment, BioSequences

const REFERENCE = "ATGCCTTGGGCCATCGGACCATATGTTTACGATGGGCAGCTGACTACCGACAACCGTCAATTCGTCTCAGAGAAGTAA"

function replace_at(sequence, position, base)
    chars = collect(sequence); chars[position] = base; join(chars)
end

cases = [
    ("exact", REFERENCE), ("synonymous", replace_at(REFERENCE, 12, 'T')),
    ("codon_deletion", REFERENCE[1:30] * REFERENCE[34:end]),
    ("codon_insertion", REFERENCE[1:33] * "GCT" * REFERENCE[34:end]),
    ("frameshift_deletion", REFERENCE[1:31] * REFERENCE[33:end]),
    ("ambiguous", replace_at(REFERENCE, 20, 'N')), ("late_start", "TTG" * REFERENCE[4:end]),
    ("early_stop", REFERENCE[1:end-3]), ("bad_match", "ATG" * repeat("GCT", 26) * "TAA"),
    ("flanking_sequence", "CC" * REFERENCE * "GG"),
]

function read_records(path)
    records = Tuple{String,String}[]; name = ""; sequence = IOBuffer()
    for line in eachline(path)
        if startswith(line, ">")
            !isempty(name) && push!(records, (name, String(take!(sequence))))
            name = line[2:end]
        else
            write(sequence, strip(line))
        end
    end
    !isempty(name) && push!(records, (name, String(take!(sequence))))
    records
end

mktempdir() do directory
    reference_path = joinpath(directory, "reference.fasta"); query_path = joinpath(directory, "queries.fasta")
    functional_path = joinpath(directory, "functionals.fasta"); rejected_path = joinpath(directory, "nonfunctionals.fasta")
    open(reference_path, "w") do io; println(io, ">reference\n", REFERENCE); end
    open(query_path, "w") do io
        for (name, sequence) in cases; println(io, ">", name, "\n", sequence); end
    end
    redirect_stdout(devnull) do
        SeededAlignment.filter_and_align(reference_path, query_path, functional_path, rejected_path; match_thresh=0.7)
    end
    functional = read_records(functional_path); rejected = isfile(rejected_path) ? read_records(rejected_path) : Tuple{String,String}[]
    for (name, _) in cases
        accepted = findfirst(record -> startswith(record[1], name * " ") || record[1] == name, functional)
        if !isnothing(accepted)
            nt = replace(functional[accepted][2], "-" => ""); aa = String(BioSequences.translate(LongDNA{4}(nt)))
            println(join((name, "true", "pass", nt, aa), '\t')); continue
        end
        failed = findfirst(record -> startswith(record[1], name * " ") || record[1] == name, rejected)
        descriptor = isnothing(failed) ? "reject" : rejected[failed][1]
        reason_match = match(r"(ambiguousSymbols-reject|frameshift-reject|lateStart-reject|earlyStop-reject|badMatch-reject)", descriptor)
        reason = isnothing(reason_match) ? "reject" : reason_match.match
        println(join((name, "false", reason, "", ""), '\t'))
    end
end
