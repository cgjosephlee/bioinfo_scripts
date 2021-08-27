#!/usr/bin/env python3
'''
run this script after merge

repeatmasker_to_bed.py ref.fa.out | sort -k1,1 -k2,2n > ref.fa.out.bed
bedtools merge -d -1 -c 5,7 -o collapse -i ref.fa.out.bed > ref.fa.out.merge.bed
repeatmasker_PerBaseBest.py ref.fa.out.merge.bed > ref.fa.out.merge.PerBaseBest.bed
'''

import sys

with open(sys.argv[1]) as f:
    for line in f:
        line = line.strip()
        cols = line.split()
        c = cols[4].split(',')  # classes
        if len(c) == 1:
            print(line)
            continue
        s = [int(x) for x in cols[3].split(',')]  # scores

        # i = s.index(max(s))  # only return first if multiple
        i = [i for i, j in enumerate(s) if j == max(s)]  # return all maximum
        if len(i) != 1 and len(set([c[x] for x in i])) != 1:
            print(f'[WARNING]\t{line}', file=sys.stderr)
        i = i[0]

        cols[3] = str(s[i])
        cols[4] = c[i]
        print('\t'.join(cols))
