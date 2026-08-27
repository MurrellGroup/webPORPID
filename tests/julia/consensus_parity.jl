import RobustAmpliconDenoising

mutable struct Generator
    state::UInt32
end

function next!(generator::Generator)
    generator.state = generator.state * UInt32(1664525) + UInt32(1013904223)
    generator.state
end

below!(generator::Generator, limit::Int) = Int((next!(generator) >> 8) % UInt32(limit))
base!(generator::Generator) = "ACGT"[below!(generator, 4) + 1]

function template_sequence(generator::Generator, length::Int)
    join(base!(generator) for _ in 1:length)
end

function noisy(generator::Generator, source::String, deletion::Int, insertion::Int, substitution::Int)
    output = IOBuffer()
    for value in source
        event = below!(generator, 10000)
        event < deletion && continue
        event < deletion + insertion && print(output, base!(generator))
        if event < deletion + insertion + substitution
            replacement = base!(generator)
            while replacement == value
                replacement = base!(generator)
            end
            print(output, replacement)
        else
            print(output, value)
        end
    end
    String(take!(output))
end

function emit(name::String, reads::Vector{String})
    draft = RobustAmpliconDenoising.consensus_seq(reads)
    draft2 = RobustAmpliconDenoising.refine_ref(draft, reads)
    final = RobustAmpliconDenoising.refine_ref(draft2, reads)
    _, _, matches, _ = RobustAmpliconDenoising.get_matches(final, reads, 0)
    minimum_agreement = round(minimum(matches[2:end-1]); digits=2)
    println(name, '\t', minimum_agreement, '\t', final)
end

generator = Generator(UInt32(0x6d2b79f5))
source = template_sequence(generator, 720)
emit("identical", fill(source, 5))
emit("substitution", [noisy(generator, source, 0, 0, 220) for _ in 1:9])
emit("mixed_indel", [noisy(generator, source, 120, 120, 280) for _ in 1:11])
emit("high_indel", [noisy(generator, source, 350, 350, 400) for _ in 1:13])

terminal = String[]
for index in 0:8
    read = noisy(generator, source, 100, 100, 180)
    index < 5 && (read = "ACG" * read)
    index % 3 == 0 && length(read) > 8 && (read = read[6:end])
    push!(terminal, read)
end
emit("terminal_overhang", terminal)
