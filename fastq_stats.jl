#!/usr/bin/env julia

using Printf
using Statistics
using Klib

const partitions = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9]

function get_Nxx(seq_lengths::Array{Int,1}, slen::Integer)
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

function main(args)
    if length(args) == 0
        println("Usage: fastq_stats.jl <in.fq.gz>")
        return
    end

    fx = Klib.FastxReader(Klib.GzFile(args[1]))
    n, slen = 0, 0
    seq_lengths = Int[]
    while (r = read(fx)) != nothing
        n += 1
        push!(seq_lengths, sizeof(r.seq))
        slen += sizeof(r.seq)
    end
    sort!(seq_lengths, rev=true)
    Nxx_len = get_Nxx(seq_lengths, slen)

    # println(n, "\t", slen)
    # println(seq_lengths[1],"\t", seq_lengths[end])
    # println(Nxx_len)
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
    @printf "# total_base seq_num mean max min N50 N90\n"
    @printf "%d\t%d\t%.1f\t%d\t%d\t%d\t%d\t\n" slen n slen/n seq_lengths[1] seq_lengths[end] Nxx_len[2,9] Nxx_len[2,17]
    @printf "\n"
    @printf "# Nxx count Nxx_len accu_len\n"
    for (n,(a,b,c)) in zip(partitions, eachcol(Nxx_len))
        @printf "N%d\t%d\t%d\t%d\n" n*100 a b c
    end

    # plot
    # todo
end

main(ARGS)
