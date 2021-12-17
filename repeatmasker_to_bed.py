#!/usr/bin/env python3
'''
Convert RepeatMasker output "ref.fa.out" to bed format.
Read "ref.fa.align" and output Kimura divergence if exist.
Output may need to be sorted.

Usage:
repeatmasker_to_bed.py ref.fa.out | sort -k1,1 -k2,2n > ref.fa.out.bed
'''

import sys

in_out = sys.argv[1]
in_align = in_out[:-4] + '.align'

try:
    f = open(in_align)
    has_align = True
except FileNotFoundError:
    has_align = False

# find kimura divergence value by (ID, SW_score)
kimuraDiv = {}
if has_align:
    for line in f:
        line = line.strip().split()
        if len(line) == 0:
            continue
        # key = (ID, score)
        k = (line[-1], line[0])

        # search Kimura
        line = next(f)
        while not line.startswith('Matrix'):
            line = next(f)
        line = next(f)
        if line.startswith('Kimura'):
            kimuraDiv[k] = line.strip()[26:]
            next(f)
        next(f)
    f.close()

with open(in_out) as f:
    # skip first 3 lines
    next(f)
    next(f)
    next(f)
    for line in f:
        line = line.strip().split()
        if line[10] in ('Simple_repeat', 'Low_complexity', 'tRNA', 'rRNA'):
            continue
        if line[8] == '+':
            s = '+'
        else:
            s = '-'
        try:
            c, f = line[10].split('/')
        except:
            c, f = line[10], '.'
        # mimic BED12 format, strandness at 6th column
        # chr start end ID score strand class family match %div (%Kimura)
        if has_align:
            k = kimuraDiv[(line[14], line[0])]
            print(f'{line[4]}\t{int(line[5])-1}\t{line[6]}\t{line[14]}\t{line[0]}\t{s}\t{c}\t{f}\t{line[9]}\t{line[1]}\t{k}')
        else:
            print(f'{line[4]}\t{int(line[5])-1}\t{line[6]}\t{line[14]}\t{line[0]}\t{s}\t{c}\t{f}\t{line[9]}\t{line[1]}')
