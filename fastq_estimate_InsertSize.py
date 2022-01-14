#!/usr/bin/env python3
'''
Note:
Mapping rate might not be accurate, due to more stringent criteria for filtering proper pairs and minimap2 itself is not designed for short reads.
But it is enough for fast estimating insert size.
'''

import sys
import argparse
import statistics
import mappy

def map_in_proper_pair(r1, r2, minQ=10):
    # map to same contig
    if not r1.ctg == r2.ctg:
        return False
    # on different strand
    if not r1.strand != r2.strand:
        return False
    # MAPQ
    if not (r1.mapq >= minQ and r2.mapq >= minQ):
        return False
    # inward orientation, FR
    if not ((r1.strand == 1 and r1.r_st <= r2.r_st) or \
            (r1.strand == -1 and r1.r_en >= r2.r_en)):
        return False
    return True

def get_insert_size(r1, r2, minQ=10, maxInsertSize=3000):
    if map_in_proper_pair(r1, r2, minQ):
        if r1.strand == 1:
            l = r2.r_en - r1.r_st
        else:
            l = r1.r_en - r2.r_st
        # if l <= 100:
        #     print(f'{r1.strand} {r1.r_st} {r1.r_en} {r2.strand} {r2.r_st} {r2.r_en}')
        if l <= maxInsertSize:
            return l
    return -1

parser = argparse.ArgumentParser(description='Estimate insert size using mappy.')
parser.add_argument('ref', type=str, help='reference, .fa or .mmi')
parser.add_argument('fq1', type=str, help='read1, .fq or .fa')
parser.add_argument('fq2', type=str, help='read2, .fq or .fa')
parser.add_argument('-n', type=int, metavar='INT', default=1000, help='find n proper pairs (%(default)s)')
parser.add_argument('-v', action='store_true', help='print values')
args = parser.parse_args()

ref = args.ref
refIdx = mappy.Aligner(ref, preset='sr', best_n=1)
if not refIdx: raise Exception('ERROR: failed to load/build index')

fq1 = mappy.fastx_read(args.fq1)
fq2 = mappy.fastx_read(args.fq2)

tot = 0
n = 0
values = []
while n < args.n:
    try:
        _, r1, _ = next(fq1)
        _, r2, _ = next(fq2)
        tot += 1
    except StopIteration:
        print('WARN: no enough reads', file=sys.stderr)
        break
    m = refIdx.map(r1, seq2=r2)
    try:
        l = get_insert_size(next(m), next(m))
    except StopIteration:
        continue
    if l != -1:
        values.append(l)
        n += 1
        print(f'{n/args.n:.1%}', end='\r', file=sys.stderr)

print(f'Total:  {tot}', file=sys.stderr)
print(f'Mapped: {n}', file=sys.stderr)
print(f'Rate:   {n/tot:.1%}', file=sys.stderr)  # might no be accurate
print(f'Min:    {min(values)}', file=sys.stderr)
print(f'Max:    {max(values)}', file=sys.stderr)
print(f'Mean:   {statistics.mean(values):.1f}', file=sys.stderr)
print(f'Median: {statistics.median(values)}', file=sys.stderr)

if args.v:
    for v in values:
        print(v)
