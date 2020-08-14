#!/usr/bin/env python3
'''
Run pilon iterativily. Will detect finished runs.

Require:
    java 8
    pilon
    bwa/smalt
    samblaster
    samtools
'''

import sys
import os
import re
import argparse
import time
import subprocess as sp

parser = argparse.ArgumentParser(description='Run pilon iterativily. Will detect finished runs.')
parser.add_argument('fasta')
parser.add_argument('fq1')
parser.add_argument('fq2')
parser.add_argument('-j', type=str, help='pilon jar', required=True)
parser.add_argument('-m', type=int, help='max mem for java in Gb (%(default)s)', default=64)
parser.add_argument('-i', type=int, help='number of iterations (%(default)s)', default=5)
parser.add_argument('-p', type=str, help='output prefix (input basename)', default=None)
parser.add_argument('-t', type=int, help='threads (%(default)s)', default=20)
parser.add_argument('--aligner', type=str, help='bwa, bwa2 or smalt (%(default)s)', default='bwa')
parser.add_argument('--exec', type=str, help='path of aligner executable')
args = parser.parse_args()

fa_base = args.fasta
fq1 = os.path.realpath(args.fq1)
fq2 = os.path.realpath(args.fq2)
pilon_jar = args.j
maxMem = args.m
iter_n = args.i
prefix = args.p
threads = args.t
aligner = args.aligner
exec = args.exec

if not fa_base == os.path.basename(fa_base):
    os.chdir(os.path.dirname(fa_base))
    fa_base = os.path.basename(fa_base)
if not exec:
    # assuming in $path
    if aligner == 'bwa':
        exec = 'bwa'
    elif aligner == 'bwa2':
        exec = 'bwa-mem2'
    elif aligner == 'smalt':
        exec = 'smalt'
if not prefix:
    prefix = re.sub('.fa$|.fasta$', '', os.path.basename(fa_base)) + '_pilon'
else:
    prefix = args.p

# check program
def check_prog(cmd):
    p = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.STDOUT)
    o = p.communicate()[0]
    if p.returncode == 127:
        print(o.decode(), file=sys.stderr)
        raise FileNotFoundError()
    else:
        return o.decode()
# java? pilon?
o = check_prog([exec])
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

    with open('pilon.{}.log'.format(fa_out), 'w') as ERR:
        if not os.path.isfile(fa_out + '.bam'):
            if aligner == 'bwa':
                sp.run([exec, 'index', fa_in], check=True, stderr=ERR)
                sp.run('{} mem -v 1 -t {} {} {} {} | samblaster | samtools sort -@ {} -O bam -o {}.bam -'.format(exec, threads, fa_in, fq1, fq2, threads, fa_out), shell=True, check=True, stderr=ERR)
                sp.run('rm {0}.amb {0}.ann {0}.bwt {0}.pac {0}.sa'.format(fa_in), shell=True, check=True)
            elif aligner == 'bwa2':
                sp.run([exec, 'index', fa_in], check=True, stdout=ERR, stderr=sp.STDOUT)
                sp.run('{} mem -v 1 -t {} {} {} {} | samblaster | samtools sort -@ {} -O bam -o {}.bam -'.format(exec, threads, fa_in, fq1, fq2, threads, fa_out), shell=True, check=True, stderr=ERR)
                sp.run('rm {0}.0123 {0}.amb {0}.ann {0}.bwt.2bit.64 {0}.bwt.8bit.32 {0}.pac'.format(fa_in), shell=True, check=True)
            elif aligner == 'smalt':
                sp.run('{0} index {1} {1}'.format(exec, fa_in), shell=True, check=True, stderr=ERR)
                sp.run('{} map -n {} {} {} {} | samblaster | samtools sort -@ {} -O bam -o {}.bam -'.format(exec, threads, fa_in, fq1, fq2, threads, fa_out), shell=True, check=True, stderr=ERR)
                sp.run('rm {0}.sma {0}.smi'.format(fa_in), shell=True, check=True)
            sp.run('samtools index {}.bam'.format(fa_out), shell=True, check=True, stderr=ERR)

        sp.run('java -Xms512m -Xmx{}G -jar {} --genome {} --frags {}.bam --output {} --changes --threads {}'.format(maxMem, pilon_jar, fa_in, fa_out, fa_out, threads), shell=True, check=True, stdout=ERR)
        fa_out = fa_out + '.fasta'

    print('[{}] Iteration {} finish! (Elapsed: {} sec)'.format(
        time.strftime("%H:%M:%S", time.localtime()),
        n + 1,
        int(time.time() - startTime)), file=sys.stderr)
