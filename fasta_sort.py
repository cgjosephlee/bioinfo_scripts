#!/usr/bin/env python3
'''
Sort fasta by length.
No rename. No wrapping or unwrapping.
'''

import sys

fin = sys.argv[1]
data = {}

with open(fin) as f:
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            h = line[1:]
            data[h] = {'seq': [], 'length': 0}
        else:
            data[h]['seq'].append(line)
            data[h]['length'] += len(line)

data = sorted(data.items(), key=lambda x: x[1]['length'], reverse=True)
# return list of tuples

for k, v in data:
    print('>{}\n{}'.format(k, '\n'.join(v['seq'])))

