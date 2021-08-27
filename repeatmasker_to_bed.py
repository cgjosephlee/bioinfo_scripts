#!/usr/bin/env python3
'''
input: ref.fa.out
output may need to be sorted

repeatmasker_to_bed.py ref.fa.out | sort -k1,1 -k2,2n > ref.fa.out.bed
'''

import sys

with open(sys.argv[1]) as f:
    # skip first 3 lines
    next(f)
    next(f)
    next(f)
    for line in f:
        line = line.strip().split()
        if line[8] == '+':
            s = '+'
        else:
            s = '-'
        try:
            c, f = line[10].split('/')
        except:
            c, f = line[10], '.'
        # mimic BED12 format, strand on 6th column
        # chr start end ID score strand class family match
        print(f'{line[4]}\t{int(line[5])-1}\t{line[6]}\t{line[14]}\t{line[0]}\t{s}\t{c}\t{f}\t{line[9]}')
