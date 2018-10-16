#!/usr/bin/env python3

### for nucleotide only! ###

import sys
import io
# import re
from Bio import AlignIO, SeqIO
import argparse
import pandas as pd
import subprocess as sp
from threading import Thread
from queue import Queue

def usage():
    parser = argparse.ArgumentParser(description='generate pairwise sequence identity from fasta (nt only, require muscle)')
    parser.add_argument('fasta', help='fasta')
    parser.add_argument('-aln', action='store_true', help='input is aligned fasta')
    parser.add_argument('-t', type=int, default=1, help='threads for alignment (default: 1)')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-tg', action='store_true', help='output wide table for global identities (default: long table)')
    group.add_argument('-tl', action='store_true', help='output wide table for local identities (default: long table)')
    return parser.parse_args()

def pair_identity(s1, s2):
    # convert N into gap, and ignore positions consisting only gaps in both sequences
    # global: consider gap in either one of sequences as mismatch
    # local: ignore all gaps
    aln_len = len(s1)
    s1_len = len(s1.seq.ungap('-'))
    s2_len = len(s2.seq.ungap('-'))
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
    # glb_id = '{:.4f}'.format(match/glb_len)
    # loc_id = '{:.4f}'.format(match/loc_len)
    glb_id = round(match / glb_len, 4)
    loc_id = round(match / loc_len, 4)
    return [s1.id, s2.id, s1_len, s2_len, glb_len, glb_id, loc_len, loc_id]

def pairwise_alignment(s1, s2):
    records = [s1, s2]
    proc = sp.Popen(['muscle', '-quiet'], stdin=sp.PIPE, stdout=sp.PIPE, universal_newlines=True)
    SeqIO.write(records, proc.stdin, 'fasta')
    outs = io.StringIO(proc.communicate()[0])
    # assuming the order of muscle outputs would not change while only 2 sequences
    aln_fa = AlignIO.read(outs, 'fasta')
    out = pair_identity(aln_fa[0], aln_fa[1])
    return out

def worker():
    while not queue.empty():
        s1, s2, n = queue.get()
        outs[n] = pairwise_alignment(s1, s2)

args = usage()
IN = args.fasta
NUM_THREADS = args.t

# main process
if not args.aln:
    FA = [x for x in SeqIO.parse(IN, 'fasta')]
    seq_num = len(FA)

    job_n = 0
    queue = Queue()
    for n1 in range(seq_num):
        for n2 in range(n1, seq_num):
            queue.put([FA[n1], FA[n2], job_n])
            job_n += 1

    total_runs = queue.qsize()
    outs = [None] * queue.qsize()
    threads = [Thread(target=worker) for x in range(NUM_THREADS)]
    print('open {} threads, running {} jobs...'.format(NUM_THREADS, total_runs), file=sys.stderr)
    for th in threads:
        th.start()
    for th in threads:
        th.join()
    print('finish!', file=sys.stderr)
else:
    FA = AlignIO.read(IN, 'fasta')
    seq_num = len(FA)

    outs = []
    for n1 in range(seq_num):
        for n2 in range(n1, seq_num):
            outs.append(pair_identity(FA[n1], FA[n2]))

glb_val = []
loc_val = []
for o in outs:
    if o[0] != o[1]:
        glb_val.append(o[5])
        loc_val.append(o[7])
glb_mean = sum(glb_val) / len(glb_val)
loc_mean = sum(loc_val) / len(loc_val)

# output results
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
          '# mean of {} identity: {:.4f}'.format(seq_num, s, i))
    df.to_csv(sys.stdout, sep='\t')
else:
    print('# number of sequences: {}\n'
          '# mean of global identity: {:.4f}\n'
          '# mean of local identity: {:.4f}\n'
          '# seq1\tseq2\tlen1\tlen2\tglb_len\tglb_id\tloc_len\tloc_id'.format(seq_num, glb_mean, loc_mean))
    for out in outs:
        if out[0] != out[1]:
            print('\t'.join(str(x) for x in out))
