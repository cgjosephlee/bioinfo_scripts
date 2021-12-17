#!/usr/bin/env python3
'''
Solve overlaps and get best record for each base by descending score.

Region      Score
-----------------
======          3
  ======        1
    ======      4
      ======    2
-----------------
333344444422

Usage:
repeatmasker_to_bed.py ref.fa.out | sort -k1,1 -k2,2n > ref.fa.out.bed
repeatmasker_PerBaseBest.py ref.fa.out.bed > ref.fa.out.PerBaseBest.bed

Note:
https://github.com/lh3/cgranges
https://github.com/biocore-ntnu/ncls
https://github.com/biocore-ntnu/pyranges
'''

import sys
from intervaltree import Interval, IntervalTree

def substract_region(Ast, Aed, Bst, Bed) -> list:
    # return A - B
    # A =====
    # B   =====
    if Ast < Bst and Aed <= Bed:
        return [(Ast, Bst)]
    # A   =====
    # B =====
    elif Ast >= Bst and Aed > Bed:
        return [(Bed, Aed)]
    # A =========
    # B   =====
    elif Ast < Bst and Aed > Bed:
        return [(Ast, Bst), (Bed, Aed)]
    # A   =====
    # B =========
    elif Ast >= Bst and Aed <= Bed:
        return []
    # A       =====
    # B =====  or   =====
    elif (Ast <= Bst and Aed <= Bst) or (Ast >= Bed and Aed >= Bed):
        return [(Ast, Aed)]
    else:
        raise ValueError((Ast, Aed, Bst, Bed))

def substract_many_regions(Ast: int, Aed: int, BRegions: list) -> list:
    ARegions = [(Ast, Aed)]
    for Bst, Bed in BRegions:
        tmp = []
        for Ast, Aed in ARegions:
            tmp += substract_region(Ast, Aed, Bst, Bed)
        ARegions = tmp
    return ARegions

def overlaps_to_regions(ovl: set) -> list:
    r = []
    while ovl:
        o = ovl.pop()  # o: Interval
        r.append((o.begin, o.end))
    return sorted(r)

def solve_overlaps(lines: list) -> IntervalTree:
    t = IntervalTree()
    # insert to tree by descending score
    lines.sort(key=lambda x: int(x[4]), reverse=True)
    for line in lines:
        st = int(line[1])
        ed = int(line[2])
        ovl = t.overlap(st, ed)
        if ovl:
            # remove overlapped regions and insert
            BRegions = overlaps_to_regions(ovl)
            ARegions = substract_many_regions(st, ed, BRegions)
            for Ast, Aed in ARegions:
                t.addi(Ast, Aed, line)
        else:
            t.addi(st, ed, line)
    return t

def generate_output(t: IntervalTree):
    for st, ed, line in sorted(t):
        info = "\t".join(line[3:])
        print(f'{line[0]}\t{st}\t{ed}\t{info}')

if __name__ == '__main__':
    with open(sys.argv[1]) as f:
        # initialize
        for line in f:
            line = line.strip().split()
            lines = [line]
            chrom = line[0]
            break
        # iterate
        for line in f:
            line = line.strip().split()
            if line[0] == chrom:
                lines.append(line)
            else:
                t = solve_overlaps(lines)
                generate_output(t)
                lines = [line]
                chrom = line[0]
        t = solve_overlaps(lines)
        generate_output(t)
