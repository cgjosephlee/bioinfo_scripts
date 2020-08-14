#!/usr/bin/env python3

import sys
import time
import numpy as np

fname = sys.argv[1]
prob_diff = float(sys.argv[2])  # abs(prob_1 - prob_0) > prob_diff
min_depth = 4

cutoff = (0.5 + prob_diff * 0.5, 0.5 - prob_diff * 0.5)
print('prob_cutoff = {}'.format(cutoff), file=sys.stderr)

print('[{}] Start reading file...'.format(time.strftime('%H:%M:%S', time.localtime())), file=sys.stderr)

fin = open(fname)

all_freqs = {}
for line in fin:
    line = line.strip().split()
    prob = float(line[7])
    if prob > cutoff[0]:
        methyl = True
    elif prob < cutoff[1]:
        methyl = False
    else:
        continue

    chr = line[0]
    pos = int(line[1])
    if line[2] == '-':
        pos -= 1  # position of C

    # [total, methaylated]
    if methyl:
        all_freqs[(chr, pos)] = all_freqs.setdefault((chr, pos), np.array([0, 0])) + [1, 1]
    else:
        all_freqs[(chr, pos)] = all_freqs.setdefault((chr, pos), np.array([0, 0])) + [1, 0]
fin.close()

print('[{}] Post-manipulation...'.format(time.strftime('%H:%M:%S', time.localtime())), file=sys.stderr)

all_freqs = all_freqs.items()
# filter depth
all_freqs = [x for x in all_freqs if x[1][0] >= min_depth]
# sort
all_freqs.sort(key=lambda x: (x[0][0], x[0][1]))

print('[{}] Writing output...'.format(time.strftime('%H:%M:%S', time.localtime())), file=sys.stderr)

for k, v in all_freqs:
    print('{}\t{}\t{}\t{}\t{:.3f}\t{}\t{}'
          .format(
              k[0],
              k[1],
              k[1] + 1,
              'C',
              v[1] / v[0],
              v[1],
              v[0]
          ))

print('[{}] Finish.'.format(time.strftime('%H:%M:%S', time.localtime())), file=sys.stderr)
