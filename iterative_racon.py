#!/usr/bin/env python3

'''
Run racon iterativily. Will detect finished runs.

Require:
    racon
    minimap2
'''

import sys
import os
import re
import argparse
import subprocess as sp

parser = argparse.ArgumentParser(description='Run racon iterativily with nanopore reads. Will detect finished runs.')
parser.add_argument('fasta')
parser.add_argument('fq')
parser.add_argument('-p', type=str, help='output prefix (input basename)', default=None)
parser.add_argument('-i', type=int, help='number of iterations (3)', default=3)
parser.add_argument('-t', type=int, help='threads (20)', default=20)
parser.add_argument('--racon', type=str, help='racon executable if not in $PATH', default='racon')
parser.add_argument('--minimap2', type=str, help='minimap2 executable if not in $PATH', default='minimap2')
args = parser.parse_args()

RACON = args.racon
MINIMAP = args.minimap2
mapper = 'map-ont'

fa_base = args.fasta
fq = args.fq
if not args.p:
    prefix = re.sub('.fa$|.fasta$', '', os.path.basename(fa_base)) + '_racon'
else:
    prefix = args.p
iter_n = args.i
threads = args.t

fa_out = fa_base
for n in range(iter_n):
    fa_in = fa_out
    fa_out = '{}_{}'.format(prefix, n + 1)
    if os.path.isfile(fa_out + '.fasta') and os.path.getsize(fa_out + '.fasta') > 0:
        print('Found iteration {} result, skip.'.format(n + 1), file=sys.stderr)
        fa_out = fa_out + '.fasta'
        continue

    aln = 'reads{}.paf'.format(n + 1)
    print('Iteration {} start...'.format(n + 1), file=sys.stderr)
    if not os.path.isfile(aln):
        sp.run([MINIMAP, '-x', mapper, '-t', threads, fa_in, fq], stdout=aln, stderr='minimap.{}.err'.format(n + 1), check=True)

    sp.run([RACON, '-t', threads, aln, fa_in], stdout=fa_out, stderr='racon.{}.err'.format(n + 1), check=True)
    print('Iteration {} finish!'.format(n + 1), file=sys.stderr)
    fa_out = fa_out + '.fasta'
