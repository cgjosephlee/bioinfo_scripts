#!/usr/bin/env python3

import sys
import os
import re
import argparse
import subprocess as sp

parser = argparse.ArgumentParser(description='')
parser.add_argument('fasta', type=str,
                    help='fasta file')
parser.add_argument('-t', type=int, default=8,
                    help='threads for nucmer')
args = parser.parse_args()

in_fa = args.fasta
threads = args.t

if not os.path.isfile('{}.delta'.format(in_fa)):
    print('Running nucmer...', file=sys.stderr)
    cmd = 'nucmer -p {0} -t {1} {0} {0}'.format(in_fa, threads)
    sp.run(cmd, shell=True, check=True)

cmd = ['show-coords', '-HrcoT', in_fa + '.delta']
proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True)
out, err = proc.communicate()
# print(out, err)
if proc.returncode == 0:
    excluded_seqs = []
    for line in out:
        # print(line)
        line = line.strip().split()
        if len(line) == 12 and line[11] == '[CONTAINS]':
            excluded_seqs.append(line[10])
else:
    raise sp.CalledProcessError(err)

out_fa = re.sub('.fa$|.fasta$', '', in_fa) + '_rmdup.fasta'
excluded_seqs = set(excluded_seqs)
if len(excluded_seqs) > 0:
    print('Removing {} sequences.'.format(len(excluded_seqs)), file=sys.stderr)
    # print(excluded_seqs)
    with open(in_fa) as f, open(out_fa, 'wt') as fout:
        for line in f:
            if line.startswith('>'):
                p = True
                if line.strip()[1:] in excluded_seqs:
                    print('Remove {}.'.format(line.strip()[1:]), file=sys.stderr)
                    p = False
            if p:
                print(line, end='', file=fout)
else:
    print('Nothing to remove, make a soft link.', file=sys.stderr)
    sp.run(['ln', '-snf', in_fa, out_fa], check=True)
