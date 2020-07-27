#!/usr/bin/env python3

import sys
from collections import defaultdict

def cal_Nxx(l, N):
    assert isinstance(l, list)
    assert isinstance(N, int)
    assert N > 0
    assert N < 100
    l.sort(reverse=True)
    threshold = sum(l) * N * 0.01
    accu_len = 0
    idx = 0
    while accu_len < threshold:
        accu_len += l[idx]
        idx += 1
    else:
        return l[idx-1]

infile = sys.argv[1]
minLen = int(sys.argv[2])

d = defaultdict(list)
tot_seq = 0
tot_base = 0
with open(infile) as f:
    h = next(f).strip().split()
    assert len(h) == 39
    assert h[20] == 'barcode_arrangement'

    for line in f:
        line = line.strip().split()
        if int(line[13]) >= minLen:
            d[line[20]].append(int(line[13]))  # length is after trimming barcode
            tot_seq += 1
            tot_base += int(line[13])

print('barcode\tseq\tbase\tseq_ratio\tbase_ratio\tmax\tmin\tN50')
for k in sorted(d.keys()):
    print('{}\t{}\t{}\t{:.2%}\t{:.2%}\t{}\t{}\t{}'.format(
        k,
        len(d[k]),
        sum(d[k]),
        len(d[k]) / tot_seq,
        sum(d[k]) / tot_base,
        max(d[k]),
        min(d[k]),
        cal_Nxx(d[k], 50)
    ))
print(f'total\t{tot_seq}\t{tot_base}')
