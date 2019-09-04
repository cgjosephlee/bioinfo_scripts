#!/usr/bin/env python3
'''
Sort a fasta.
No rename. No wrapping or unwrapping.

NOTE:
A more advanced implementation: http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec:SeqIO-sort
'''

import argparse
import re

parser = argparse.ArgumentParser(description='Sort a fasta.')
parser.add_argument('fasta', type=str)
parser.add_argument('--by', type=str, choices=['l', 'n', 'N'], default='l',
                    help='l: descending length; n: alphabetically ordered ID; N: naturally ordered ID (%(default)s)')
args = parser.parse_args()

fin = args.fasta
key = args.by

seq = {}
seq_len = {}
with open(fin) as f:
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            h = line[1:]
            seq[h] = []
            seq_len[h] = 0
        else:
            seq[h].append(line)
            seq_len[h] += len(line)

# sort by descending length
if key == 'l':
    sorted_ID = [x for x,_ in sorted(seq_len.items(), reverse=True, key=lambda x: x[1])]
# lexical sort by ID
elif key == 'n':
    sorted_ID = sorted(seq.keys())
# natural sort by ID, case insensitive
elif key == 'N':
    # https://stackoverflow.com/a/16090640/7859425
    sorted_ID = sorted(seq.keys(),
                       key=lambda s: [int(t) if t.isdigit() else t.lower() for t in re.split('(\d+)', s)])

for i in sorted_ID:
    print('>{}\n{}'.format(i, '\n'.join(seq[i])))
