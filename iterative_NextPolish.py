#!/usr/bin/env python3
'''
Run NextPolish for short reads iterativily. Will detect finished runs.

Require:
    bwa
    samblaster
    samtools
'''

import sys
import os
import re
import argparse
import time
import subprocess as sp

parser = argparse.ArgumentParser(description='Run NextPolish iterativily. Will detect finished runs.')
parser.add_argument('fasta')
parser.add_argument('fq1')
parser.add_argument('fq2')
parser.add_argument('-i', type=int, help='number of iterations (%(default)s)', default=2)
parser.add_argument('-p', type=str, help='output prefix (input basename)', default=None)
parser.add_argument('-t', type=int, help='threads (%(default)s)', default=20)
parser.add_argument('--dir', type=str, help='path to NextPolish dir', required=True)
parser.add_argument('--aligner', type=str, help='bwa, bwa2 (%(default)s)', default='bwa2')
parser.add_argument('--exec', type=str, help='path of aligner executable')
args = parser.parse_args()

fa_base = args.fasta
fq1 = os.path.realpath(args.fq1)
fq2 = os.path.realpath(args.fq2)
iter_n = args.i
prefix = args.p
threads = args.t
aligner = args.aligner

if not fa_base == os.path.basename(fa_base):
    os.chdir(os.path.dirname(fa_base))
    fa_base = os.path.basename(fa_base)

if not prefix:
    prefix = re.sub('.fa$|.fasta$', '', os.path.basename(fa_base)) + '_NextPolish'
else:
    prefix = args.p

ALIGNER = args.exec
if not ALIGNER:
    # assuming in $path
    if aligner == 'bwa':
        ALIGNER = 'bwa'
    elif aligner == 'bwa2':
        ALIGNER = 'bwa-mem2'

NP1 = os.path.join(args.dir, 'lib/nextpolish1.py')

# check program
def check_prog(cmd):
    p = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.STDOUT)
    o = p.communicate()[0]
    # 127: command not found; 2: python find no file
    if p.returncode == 127 or p.returncode == 2:
        print(o.decode(), file=sys.stderr)
        raise FileNotFoundError()
    else:
        return o.decode()

o = check_prog(['python', NP1, '-h'])
o = check_prog([ALIGNER])
o = check_prog(['samblaster', '-h'])
o = check_prog(['samtools'])

fa_out = fa_base
for n in range(iter_n):
    fa_in = fa_out
    fa_out = '{}_{}'.format(prefix, n + 1)
    if os.path.isfile(fa_out + '.fasta') and os.path.getsize(fa_out + '.fasta') > 0:
        print('Found iteration {} result, skip.'.format(n + 1), file=sys.stderr)
        fa_out = fa_out + '.fasta'
        continue

    startTime = time.time()
    print('[{}] Iteration {} start...'.format(time.strftime("%H:%M:%S", time.localtime(startTime)), n + 1), file=sys.stderr)

    with open('NextPolish.{}.log'.format(fa_out), 'w') as ERR:
        ### step 1
        if aligner == 'bwa':
            sp.run([ALIGNER, 'index', '-p', fa_in + '.bwa1', fa_in], check=True, stdout=ERR, stderr=sp.STDOUT)
            sp.run(f'{ALIGNER} mem -v 1 -t {threads} {fa_in}.bwa1 {fq1} {fq2} | samblaster | samtools sort -@ {threads} -O bam -o {fa_out}.1.bam -', shell=True, check=True, stderr=ERR)
            sp.run(f'rm {fa_in}.bwa1*', shell=True, check=True)
        elif aligner == 'bwa2':
            sp.run([ALIGNER, 'index', '-p', fa_in + '.bwa2', fa_in], check=True, stdout=ERR, stderr=sp.STDOUT)
            sp.run(f'{ALIGNER} mem -v 1 -t {threads} {fa_in}.bwa2 {fq1} {fq2} | samblaster | samtools sort -@ {threads} -O bam -o {fa_out}.1.bam -', shell=True, check=True, stderr=ERR)
            sp.run(f'rm {fa_in}.bwa2*', shell=True, check=True)

        sp.run(f'samtools index {fa_out}.1.bam', shell=True, check=True, stderr=ERR)
        sp.run(f'samtools faidx {fa_in}', shell=True, check=True, stderr=ERR)

        with open(f'{fa_out}.1.fasta', 'w') as SO:
            sp.run(f'python {NP1} -t 1 -g {fa_in} -p {threads} -s {fa_out}.1.bam', shell=True, check=True, stdout=SO, stderr=ERR)
        fa_in = fa_out + '.1.fasta'

        ### step 2
        if aligner == 'bwa':
            sp.run([ALIGNER, 'index', '-p', fa_in + '.bwa1', fa_in], check=True, stdout=ERR, stderr=sp.STDOUT)
            sp.run(f'{ALIGNER} mem -v 1 -t {threads} {fa_in}.bwa1 {fq1} {fq2} | samblaster | samtools sort -@ {threads} -O bam -o {fa_out}.2.bam -', shell=True, check=True, stderr=ERR)
            sp.run(f'rm {fa_in}.bwa1*', shell=True, check=True)
        elif aligner == 'bwa2':
            sp.run([ALIGNER, 'index', '-p', fa_in + '.bwa2', fa_in], check=True, stdout=ERR, stderr=sp.STDOUT)
            sp.run(f'{ALIGNER} mem -v 1 -t {threads} {fa_in}.bwa2 {fq1} {fq2} | samblaster | samtools sort -@ {threads} -O bam -o {fa_out}.2.bam -', shell=True, check=True, stderr=ERR)
            sp.run(f'rm {fa_in}.bwa2*', shell=True, check=True)

        sp.run(f'samtools index {fa_out}.2.bam', shell=True, check=True, stderr=ERR)
        sp.run(f'samtools faidx {fa_in}', shell=True, check=True, stderr=ERR)

        with open(f'{fa_out}.fasta', 'w') as SO:
            sp.run(f'python {NP1} -t 2 -g {fa_in} -p {threads} -s {fa_out}.2.bam', shell=True, check=True, stdout=SO, stderr=ERR)
        fa_out = fa_out + '.fasta'

    print('[{}] Iteration {} finish! (Elapsed: {} sec)'.format(
        time.strftime("%H:%M:%S", time.localtime()),
        n + 1,
        int(time.time() - startTime)), file=sys.stderr)
