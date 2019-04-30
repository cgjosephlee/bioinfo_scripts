#!/usr/bin/env python3
'''
Remove duplicated entries in fastq.

Key insertion order is preserved in python 3.7+. Otherwise the order may be shuffled.
'''

import sys

try:
    fin = open(sys.argv[1])
except IOError:
    fin = sys.stdin

data = {}
total = 0
for line in fin:
    id = line.strip()[1:]

    line = next(fin)
    seq = line.strip()

    next(fin)
    line = next(fin)
    qual = line.strip()

    data[id] = (seq, qual)
    total += 1

for k, v in data.items():
    print('@{}\n{}\n+\n{}'.format(k, v[0], v[1]))

print('Done!', file=sys.stderr)
print('{}/{} entries were removed.'.format(total - len(data.keys()), total), file=sys.stderr)
