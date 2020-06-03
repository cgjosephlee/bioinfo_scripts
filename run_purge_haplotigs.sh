#!/bin/bash
set -e

input=$1
reads=$2
depth=${3:-200}
thread=16

if [[ $# == 2 || $# == 3 ]]; then
    test -s aligned.bam || minimap2 -t $thread -ax map-ont --secondary=no $input $reads | samtools sort -@ $thread -o aligned.bam
    purge_haplotigs hist -t $thread -b aligned.bam -g $input -d $depth
fi

# manual step
# check hist.png and decide coverage cutoff

if [[ $# == 5 ]]; then
    cutoff_L=$3
    cutoff_M=$4
    cutoff_H=$5
    purge_haplotigs cov -i aligned.bam.gencov -l $cutoff_L -m $cutoff_M -h $cutoff_H
    purge_haplotigs purge -t $thread -g $input -c coverage_stats.csv
fi
