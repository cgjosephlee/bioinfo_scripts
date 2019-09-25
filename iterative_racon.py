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
import time
import subprocess as sp

parser = argparse.ArgumentParser(description='Run racon iterativily with nanopore reads. Will detect finished runs. Use uncompressed fastq will fasten racon significantly.')
parser.add_argument('fasta')
parser.add_argument('fq')
parser.add_argument('-p', type=str, help='output prefix (input basename)', default=None)
parser.add_argument('-i', type=int, help='number of iterations (%(default)s)', default=4)
parser.add_argument('-t', type=int, help='threads (%(default)s)', default=20)
parser.add_argument('--racon', type=str, default='racon',
                    help='racon executable if not in $PATH')
parser.add_argument('--racon-opts', type=str, default='',
                    help='additional racon options')
parser.add_argument('--minimap2', type=str, default='minimap2',
                    help='minimap2 executable if not in $PATH')
parser.add_argument('--minimap2-format', type=str, default='sam', choices=['sam', 'paf'],
                    help='minimap2 output format (%(default)s)')
parser.add_argument('--minimap2-preset', type=str, default='map-ont',
                    help='minimap2 config preset (%(default)s)')
parser.add_argument('--minimap2-opts', type=str, default='',
                    help='additional minimap2 options')
args = parser.parse_args()

RACON = args.racon
RACON_opts = args.racon_opts.split()
MINIMAP = args.minimap2
MINIMAP_format = args.minimap2_format
MINIMAP_preset = args.minimap2_preset
MINIMAP_opts = args.minimap2_opts.split()

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
    fa_out = '{}_{}'.format(prefix, n + 1) + '.fasta'
    if os.path.isfile(fa_out) and os.path.getsize(fa_out) > 0:
        print('Found iteration {} result, skip.'.format(n + 1), file=sys.stderr)
        continue

    startTime = time.time()
    print('[{}] Iteration {} start...'.format(time.strftime("%H:%M:%S", time.localtime(startTime)), n + 1), file=sys.stderr)

    # '-c' or '-a' is crucial in minimap2 mapping approach!
    # https://github.com/lh3/minimap2/blob/master/FAQ.md#1-alignment-different-with-option--a-or--c
    if MINIMAP_format == 'sam':
        aln = 'reads{}.sam'.format(n + 1)
        format_opt = '-a'
    elif MINIMAP_format == 'paf':
        aln = 'reads{}.paf'.format(n + 1)
        format_opt = '-c'

    with open('racon.{}.log'.format(n + 1), 'w') as ERR:
        if not os.path.isfile(aln) or os.path.getsize(aln) == 0:
            with open(aln, 'w') as OUT:
                cmd = [MINIMAP, format_opt, '-x', MINIMAP_preset, '-t', str(threads)] + MINIMAP_opts + [fa_in, fq]
                sp.run(cmd, stdout=OUT, stderr=ERR, check=True)

        with open(fa_out, 'w') as OUT:
            cmd = [RACON, '-t', str(threads)] + RACON_opts + [fq, aln, fa_in]
            sp.run(cmd, stdout=OUT, stderr=ERR, check=True)

    print('[{}] Iteration {} finish! (Elapsed: {} sec)'.format(
        time.strftime("%H:%M:%S", time.localtime()),
        n + 1,
        int(time.time() - startTime)), file=sys.stderr)
