#!/bin/bash
set -e

k=21
ploidy=2  # default=2
thread=8
mem=64  # GB

while getopts "k:p:t:" opt; do
    case $opt in
        k) k=$OPTARG ;;
        p) ploidy=$OPTARG ;;
        t) thread=$OPTARG ;;
    esac
done

shift $((OPTIND-1))
fof=$1
prefix=$2

mkdir -p tmp
kmc -k${k} -t${thread} -m${mem} -ci1 -cs10000 @${fof} ${prefix}_freq tmp/
kmc_tools transform ${prefix}_freq histogram ${prefix}_hist -cx10000

genomescope.R -k ${k} -i ${prefix}_hist -o . -p $ploidy -n ${prefix}_genomescope

L=$(smudgeplot.py cutoff ${prefix}_hist L)
U=$(smudgeplot.py cutoff ${prefix}_hist U)
echo 'cutoff' $L $U # these need to be sane values
# L should be like 20 - 200
# U should be like 500 - 3000
kmc_tools transform ${prefix}_freq -ci$L -cx$U dump -s ${prefix}_dump
smudgeplot.py hetkmers -o ${prefix}_pairs ${prefix}_dump
smudgeplot.py plot -k $k -o ${prefix} ${prefix}_pairs_coverages.tsv
