#!/usr/bin/env python3

import sys
import os
import re
import argparse
import subprocess as sp

parser = argparse.ArgumentParser(description='Remove self-contained sequences in assembly.')
parser.add_argument('fasta', type=str,
                    help='fasta file')
parser.add_argument('-t', type=int, default=8,
                    help='threads')
parser.add_argument('--prog', type=str, choices=['nucmer', 'minimap2'], default='minimap2',
                    help='program')
parser.add_argument('--min-identity', type=float, default=.98,
                    help='0-1 (%(default)s)')
parser.add_argument('--min-coverage', type=float, default=.98,
                    help='0-1 (%(default)s)')
parser.add_argument('--dry', action='store_true',
                    help='don\'t run')
parser.add_argument('--redo', action='store_true',
                    help='force redo')
args = parser.parse_args()

in_fa = args.fasta
threads = args.t
prog = args.prog
min_id = args.min_identity
min_cov = args.min_coverage

if prog == 'nucmer':
    if not os.path.isfile('{}.delta'.format(in_fa)) or args.redo:
        print('Running nucmer...', file=sys.stderr)
        cmd = 'nucmer --maxmatch -p {0} -t {1} {0} {0}'.format(in_fa, threads)  # very slow!
        sp.run(cmd, shell=True, check=True)
    if not os.path.isfile('{}.delta.coords'.format(in_fa)) or args.redo:
        print('Running show-coords...', file=sys.stderr)
        cmd = 'show-coords -HrcolT {}.delta > {}.delta.coords'.format(in_fa, threads)
        sp.run(cmd, shell=True, check=True)
    else:
        print('Use existing result...', file=sys.stderr)

    with open('{}.delta.coords'.format(in_fa)) as f:
        excluded_seqs = []
        for line in f:
            line = line.strip().split()
            id = float(line[6]) / 100
            cov = float(line[10]) / 100
            if len(line) == 14 and line[-1] == '[CONTAINS]':
                if id > min_id and cov > min_cov:
                    print('\t'.join(line[:-1]))
                    excluded_seqs.append(line[-2])
elif prog == 'minimap2':
    # https://lh3.github.io/minimap2/minimap2.html
    # https://github.com/lh3/minimap2/blob/master/cookbook.md#constructing-self-homology-map
    if not os.path.isfile('{}.paf'.format(in_fa)) or args.redo:
        print('Running minimap2...', file=sys.stderr)
        # asm20 = -k19 -w10 -A1 -B4 -O6,26 -E2,1 -s200 -z200 -N50 --min-occ-floor=100
        # cmd = 'minimap2 -t {1} -D -x asm20 --cs {0} {0} > {0}.paf'.format(in_fa, threads)
        cmd = 'minimap2 -t {1} -DP -k19 -w19 -m900 --cs {0} {0} > {0}.paf'.format(in_fa, threads)
        sp.run(cmd, shell=True, check=True)
    else:
        print('Use existing result...', file=sys.stderr)

    with open('{}.paf'.format(in_fa)) as f:
        excluded_seqs = []
        for line in f:
            line = line.strip().split()
            if line[0] == line[5]:
                continue
            alnlen = int(line[10])
            id = int(line[9]) / alnlen
            tlen = int(line[6])
            tcov = (int(line[8]) - int(line[7])) / tlen
            qlen = int(line[1])
            qcov = (int(line[3]) - int(line[2])) / qlen
            if id > min_id and qcov > min_cov:
                print(
                    '\t'.join(line[0:12]),
                    'id:f:{}\tqv:f:{}'.format(round(id, 5), round(qcov, 5)), sep='\t')
                excluded_seqs.append(line[0])

excluded_seqs = set(excluded_seqs)
# print(excluded_seqs)

if not args.dry:
    out_fa = re.sub('.fa$|.fasta$', '', in_fa) + '_rmdup.fasta'
    if len(excluded_seqs) > 0:
        print('Removing {} sequences.'.format(len(excluded_seqs)), file=sys.stderr)
        print(excluded_seqs, file=sys.stderr)
        with open(in_fa) as f, open(out_fa, 'wt') as fout:
            for line in f:
                if line.startswith('>'):
                    p = True
                    if line.strip().split()[0][1:] in excluded_seqs:
                        # print('Remove {}.'.format(line.strip()[1:]), file=sys.stderr)
                        p = False
                if p:
                    print(line, end='', file=fout)
    else:
        print('Nothing to remove, make a soft link.', file=sys.stderr)
        in_fa_base = os.path.basename(in_fa)
        sp.run(['ln', '-snf', in_fa_base, out_fa], check=True)
else:
    print('Removing {} sequences.'.format(len(excluded_seqs)), file=sys.stderr)
    print(excluded_seqs, file=sys.stderr)
