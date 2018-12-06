#!/usr/bin/env python3

'''
Usage: dna2sixframe.py fasta genetic_code

0-----------------x
 1-----------------
  2--------------xx
ATCGATCGATCGATCGGAG
x-----------------3
-----------------4
xx--------------5

- Trailing residual bases are trimmed.
'''

import sys
from Bio import SeqIO

IN = sys.argv[1]
code = sys.argv[2]

for rec in SeqIO.parse(IN, 'fasta'):
    F0 = rec.seq
    F1 = F0[1:]
    F2 = F0[2:]
    F3 = F0.reverse_complement()
    F4 = F3[1:]
    F5 = F3[2:]
    sixFrames = [F0, F1, F2, F3, F4, F5]
    for i in range(6):
        s = sixFrames[i]
        r = len(s) % 3
        if r != 0:
            s = s[:-r]
        print('>{}:{}\n{}'.format(rec.id, i, s.translate(table=code)))
