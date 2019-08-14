#!/usr/bin/env python3

'''
Run pilon iterativily. Will detect finished runs.

Require:
    java
    pilon
    bwa/smalt
    samtools
'''

import sys
import os
import re
import argparse
import subprocess as sp

parser = argparse.ArgumentParser(description='Run pilon iterativily. Will detect finished runs.')
parser.add_argument('fasta')
parser.add_argument('fq1')
parser.add_argument('fq2')
parser.add_argument('-j', type=str, help='pilon jar', required=True)
parser.add_argument('-m', type=int, help='max mem for java in Gb (64)', default=64)
parser.add_argument('-i', type=int, help='number of iterations (5)', default=5)
parser.add_argument('-p', type=str, help='output prefix (input basename)', default=None)
parser.add_argument('-t', type=int, help='threads (20)', default=20)
parser.add_argument('--aligner', type=str, help='bwa or smalt (bwa)', default='bwa')
args = parser.parse_args()

fa_base = args.fasta
fq1 = args.fq1
fq2 = args.fq2
pilon_jar = args.j
aligner = args.aligner
maxMem = args.m
if not args.p:
    prefix = re.sub('.fa$|.fasta$', '', os.path.basename(fa_base)) + '_pilon'
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

    print('Iteration {} start...'.format(n + 1), file=sys.stderr)
    if not os.path.isfile(fa_out + '.bam'):
        if aligner == 'bwa':
            sp.run(['bwa', 'index', fa_in], check=True)
            sp.run('bwa mem -v 1 -t {} {} {} {} | samtools sort -@ {} -O bam -o {}.bam -'.format(threads, fa_in, fq1, fq2, threads, fa_out), shell=True, check=True)
            sp.run('rm {0}.amb {0}.ann {0}.bwt {0}.pac {0}.sa'.format(fa_in), shell=True, check=True)
        elif aligner == 'smalt':
            sp.run('smalt index {0} {0}'.format(fa_in), shell=True, check=True)
            sp.run('smalt map -n {} {} {} {} | samtools sort -@ {} -O bam -o {}.bam -'.format(threads, fa_in, fq1, fq2, threads, fa_out), shell=True, check=True)
            sp.run('rm {0}.sma {0}.smi'.format(fa_in), shell=True, check=True)
        sp.run('samtools index {}.bam'.format(fa_out), shell=True, check=True)

    sp.run('java -Xms512m -Xmx{}G -jar {} --genome {} --frags {}.bam --output {} --changes --threads {}'.format(maxMem, pilon_jar, fa_in, fa_out, fa_out, threads), shell=True, check=True)
    print('Iteration {} finish!'.format(n + 1), file=sys.stderr)
    fa_out = fa_out + '.fasta'
