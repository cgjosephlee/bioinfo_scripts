#!/usr/bin/env julia

using Printf
using Statistics
using Klib

const partitions = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9]

function get_Nxx(seq_lengths::Array{Int,1}, slen::Int)
    thresholds = partitions * slen
    # Nxx_len = zeros(Int, (length(partitions), 3))
    Nxx_len = zeros(Int, (3, length(partitions)))
    accu_len = 0
    n = 0
    for i in 1:1:length(partitions)
        while accu_len < thresholds[i]
            n += 1
            accu_len += seq_lengths[n]
        end
        # Nxx_len[i,:] = [n, seq_lengths[n], accu_len]
        Nxx_len[:,i] = [n, seq_lengths[n], accu_len]
    end
    return Nxx_len
end

function generate_QualDict()
    # phred 33
    code = Vector{UInt8}("""!"#\$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~""")
    @assert length(code) == 94
    d = Dict{UInt8, Float32}()
    for c in code
        d[c] = Float32(10) ^ ((c - UInt8(33)) / Float32(-10.0))
    end
    return d
end
const QD = generate_QualDict()

function cal_meadQ(qual::String)::Float32
    qualBytes = Vector{UInt8}(qual)
    sumProb = Float32(0)
    for i in qualBytes
        sumProb += QD[i]
    end
    avgProb = sumProb / Float32(length(qualBytes))
    Float32(-10.0) * log10(avgProb)
end

# function cal_meadQ(qual::String)::Float32
#     qualBytes = Vector{UInt8}(qual)
#     ## t1
#     # avgProb = sum(Float32(10.0) .^ ((qualBytes .- UInt8(33)) / Float32(-10.0))) / Float32(length(qualBytes))
#     ## t2
#     sumProb = Float32(0)
#     for i in qualBytes
#         sumProb += Float32(10) ^ ((i - UInt8(33)) / Float32(-10.0))
#     end
#     avgProb = sumProb / Float32(length(qualBytes))
#     Float32(-10.0) * log10(avgProb)
# end

function usage()
    println("Usage: fastq_stats.jl [-q] [-s] [--dump FILE] <in.fq.gz>")
    exit(1)
end

function main(args::Vector{String})
    ### handle opts ###

    readQ = true
    singleLine = false
    dumpToFile = false
    dumpArg = ""
    for (opt, arg) in Klib.getopt(args, "qs", ["dump="])
        if opt == "-q"
            readQ = true
        elseif opt == "-s"
            singleLine = true
        elseif opt == "--dump"
            dumpToFile = true
            dumpArg = arg
        else
            usage()
        end
    end
    if length(args) != 1
        usage()
    end

    if args[1] == "-"
        fn = stdin
    else
        fn = Klib.GzFile(args[1])
    end
    fx = Klib.FastxReader(fn)

    if dumpToFile
        if dumpArg == "-" || dumpArg == "auto"
            fnDump = open(args[1]*".dump", "w")
        else
            fnDump = open(dumpArg, "w")
        end
    end

    ### main ###

    n, slen = 0, 0
    seq_lengths = Int[]
    seq_meanQ = Float32[]
    while (r = read(fx)) != nothing
        n += 1
        push!(seq_lengths, sizeof(r.seq))
        slen += sizeof(r.seq)
        if readQ
            q = cal_meadQ(r.qual)
            push!(seq_meanQ, q)
            if dumpToFile
                write(fnDump, @sprintf "%s\t%s\t%.6f\n" r.name sizeof(r.seq) q)
            end
        end
    end
    if dumpToFile
        close(fnDump)
    end
    sort!(seq_lengths, rev=true)
    Nxx_len = get_Nxx(seq_lengths, slen)

    if singleLine
        @printf "%d\t%d\t%.1f\t%d\t%d\t%d\t%d\t\n" slen n slen/n seq_lengths[1] seq_lengths[end] Nxx_len[2,9] Nxx_len[2,17]
    else
        @printf "input file: %s\n" args[1]
        @printf "\n"
        @printf "minimum length: %d\n" seq_lengths[end]
        @printf "maximum length: %d\n" seq_lengths[1]
        @printf "total seqs:     %d\n" n
        @printf "total length:   %d (%.2f Gbp)\n" slen slen/1e9
        @printf "avg. length:    %.1f\n" slen/n
        @printf "median length:  %.1f\n" median(seq_lengths)
        @printf "N50:            %d\n" Nxx_len[2,9]
        @printf "N90:            %d\n" Nxx_len[2,17]
        @printf "\n"
        @printf "# total_base seq_num mean max min N50 N90\n"
        @printf "%d\t%d\t%.1f\t%d\t%d\t%d\t%d\t\n" slen n slen/n seq_lengths[1] seq_lengths[end] Nxx_len[2,9] Nxx_len[2,17]
        @printf "\n"
        @printf "# Nxx count Nxx_len accu_len\n"
        for (n,(a,b,c)) in zip(partitions, eachcol(Nxx_len))
            @printf "N%d\t%d\t%d\t%d\n" n*100 a b c
        end
        if readQ
            @printf "\n"
            @printf "# quality stats\n"
            @printf "Q>=20: %s\n" sum(seq_meanQ .>= 20)
            @printf "Q>=10: %s\n" sum(seq_meanQ .>= 10)
            @printf "Q>= 7: %s\n" sum(seq_meanQ .>= 7)
            @printf "Q>= 5: %s\n" sum(seq_meanQ .>= 5)
        end
    end

    # plot
    # TODO
end

main(ARGS)
