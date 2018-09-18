#!/usr/bin/env python3

### for nucleotide only! ###

import sys
import re
from Bio import AlignIO
import argparse
import pandas as pd

def usage():
    parser = argparse.ArgumentParser(description='generate pairwise sequence identity from multiple sequence alignment fasta (default: long table)')
    parser.add_argument('fasta', help='MSA fasta')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-tg', action='store_true', help='output wide table for global identities')
    group.add_argument('-tl', action='store_true', help='output wide table for identities')
    return parser.parse_args()

def pair_id(s1, s2):
    # trun N into gap, and ignore base consisting only gap in both seq
    # global: consider gaps in either one seq as mismatch
    # local: ingore all gaps
    s1_len = len(s1.seq.ungap('-'))
    s2_len = len(s2.seq.ungap('-'))
    aln_len = len(s1)
    glb_len = aln_len
    loc_len = aln_len
    match = 0
    for pos in range(aln_len):
        i, j = str(s1[pos]).upper(), str(s2[pos]).upper()
        if i == 'N': i = '-'
        if j == 'N': j = '-'
        if i == j == '-':
            glb_len -= 1
            loc_len -= 1
        elif i == '-' or j == '-':
            loc_len -= 1
        elif i == j:
            match += 1
    #glb_id = '{:.4f}'.format(match/glb_len)
    #loc_id = '{:.4f}'.format(match/loc_len)
    glb_id = round(match/glb_len, 4)
    loc_id = round(match/loc_len, 4)
    return [s1.id, s2.id, s1_len, s2_len, glb_len, glb_id, loc_len, loc_id]

args = usage()
IN = args.fasta
FA = AlignIO.read(IN, 'fasta')

aln_len = FA.get_alignment_length()
seq_num = len(FA)

outs = []
glb_val = []
loc_val = []
for n1 in range(seq_num):
    for n2 in range(n1, seq_num):
        o = pair_id(FA[n1], FA[n2])
        outs.append(o)
        if o[0] != o[1]:
            glb_val.append(o[5])
            loc_val.append(o[7])

glb_mean = sum(glb_val) / len(glb_val)
loc_mean = sum(loc_val) / len(loc_val)

if args.tg or args.tl:
    df = pd.DataFrame()
    if args.tg:
        for out in outs:
            df.at[out[1], out[0]] = out[5]
        s = 'global'
        i = glb_mean
    elif args.tl:
        for out in outs:
            df.at[out[1], out[0]] = out[7]
        s = 'local'
        i = loc_mean
    print('# number of sequences: {}\n'
          '# alignment length: {}\n'
          '# mean of {} identity: {:.4f}'.format(seq_num, aln_len, s, i))
    df.to_csv(sys.stdout, sep='\t')
else:
    print('# number of sequences: {}\n'
          '# alignment length: {}\n'
          '# mean of global identity: {:.4f}\n'
          '# mean of local identity: {:.4f}\n'
          '# seq1\tseq2\tlen1\tlen2\tglb_len\tglb_id\tloc_len\tloc_id'.format(seq_num, aln_len, glb_mean, loc_mean))
    for out in outs:
        if out[0] != out[1]:
            print('\t'.join(str(x) for x in out))
