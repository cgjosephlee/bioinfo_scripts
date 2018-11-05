#!/usr/bin/env python3

import sys
import os
import argparse
import subprocess as sp

parser = argparse.ArgumentParser(description='run iterative pilon')
parser.add_argument('fasta')
parser.add_argument('fq1')
parser.add_argument('fq2')
parser.add_argument('-j', type=str, help='pilon jar', required=True)
parser.add_argument('-p', type=str, help='output prefix', default='')
parser.add_argument('-i', type=int, help='number of iterations', default=1)
parser.add_argument('-t', type=int, help='threads', default=1)
args = parser.parse_args()

fa_base = args.fasta
fq1 = args.fq1
fq2 = args.fq2
pilon_jar = args.j
if not args.p:
    prefix = fa_base.rstrip('.fa|.fasta') + '_pilon'
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
    sp.run(['bwa', 'index', fa_in], check=True)
    sp.run('bwa mem -v 1 -t {} {} {} {} | samtools view -h -F 4 - | samtools sort -@ {} -O bam -o {}.bam -'.format(threads, fa_in, fq1, fq2, threads, fa_out), shell=True, check=True)
    sp.run('samtools index {}.bam'.format(fa_out), shell=True, check=True)
    # sp.run(['java', '-Xmx64G', '-jar', pilon_jar, '--genome', fa_in, '--frags', fa_out + '.bam', '--output', fa_out, '--changes', '--threads', threads], check=True)
    sp.run('java -Xmx64G -jar {} --genome {} --frags {}.bam --output {} --changes --threads {}'.format(pilon_jar, fa_in, fa_out, fa_out, threads), shell=True, check=True)
    sp.run('rm {0}.amb {0}.ann {0}.bwt {0}.pac {0}.sa'.format(fa_in), shell=True, check=True)
    print('Iteration {} finish!'.format(n + 1), file=sys.stderr)
    fa_out = fa_out + '.fasta'
