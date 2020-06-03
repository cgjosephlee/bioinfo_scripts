#!/bin/bash
set -e

fof=$1
prefix=$2
thread=8

mkdir -p tmp
kmc -k21 -t${thread} -m256 -ci1 -cs10000 @${fof} ${prefix}_freq tmp/
kmc_tools transform ${prefix}_freq histogram ${prefix}_hist -cx10000

genomescope.R -i ${prefix}_hist -o . -p 2 -n ${prefix}_genomescope

L=$(smudgeplot.py cutoff ${prefix}_hist L)
U=$(smudgeplot.py cutoff ${prefix}_hist U)
echo 'cutoff' $L $U # these need to be sane values
# L should be like 20 - 200
# U should be like 500 - 3000
kmc_tools transform ${prefix}_freq -ci$L -cx$U dump -s ${prefix}_dump
smudgeplot.py hetkmers -o ${prefix}_pairs ${prefix}_dump
smudgeplot.py plot -o ${prefix} ${prefix}_pairs_coverages.tsv
