#!/usr/bin/env python3

import sys
from Bio import AlignIO

FIN = sys.argv[1]
# 0-indexed
breakpoints = sorted([int(x) for x in sys.argv[2].split(',')])

if FIN.endswith('.phy') or FIN.endswith('.phylip'):
    # aln_fmt = 'phylip'
    aln_fmt = 'phylip-relaxed'
    out_suffix = 'phy'
elif FIN.endswith('.fa') or FIN.endswith('.fasta'):
    aln_fmt = 'fasta'
    out_suffix = 'fa'
else:
    raise ValueError('Unknown format.')

in_aln = AlignIO.read(FIN, aln_fmt)
breakpoints = [0] + breakpoints + [in_aln.get_alignment_length()]

for i in range(1, len(breakpoints)):
    st = breakpoints[i-1]
    ed = breakpoints[i]
    FOUT = '{}.{}..{}.{}'.format(FIN, st, ed, out_suffix)
    print(FOUT, file=sys.stderr)
    out_aln = open(FOUT, 'w')
    AlignIO.write(in_aln[:,st:ed], out_aln, aln_fmt)
    out_aln.close()
