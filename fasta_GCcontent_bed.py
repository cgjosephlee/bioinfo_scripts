#!/usr/bin/env python3
'''
Usage: fasta_GCcontent_bed.py <fasta> <window size>
Output: bed formated GC percentage and counts binned by window size
'''

import sys
from Bio import SeqIO

fin = sys.argv[1]
window = int(sys.argv[2])

for rec in SeqIO.parse(fin, 'fasta'):
    length = len(rec.seq)
    st = 0
    while st + window < length:
        # took from biopython, Ns are counted as denomirator
        gc = sum(rec.seq[st:st+window].count(x) for x in ['G', 'C', 'g', 'c', 'S', 's'])
        try:
            p = gc / window
        except ZeroDivisionError:
            p = 0.0
        print('{}\t{}\t{}\t{:.4f}\t{}'.format(rec.id, st, st + window, p, gc))
        st += window
    else:
        gc = sum(rec.seq[st:].count(x) for x in ['G', 'C', 'g', 'c', 'S', 's'])
        try:
            p = gc / (length - st)
        except ZeroDivisionError:
            p = 0.0
        print('{}\t{}\t{}\t{:.4f}\t{}'.format(rec.id, st, length, p, gc))
