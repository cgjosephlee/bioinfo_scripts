#!/usr/bin/env python3
'''
Parse trf output and find potential telomeric repeats.
Require trf and samtools.
'''

import sys
import os
import re
import argparse
import subprocess

def rc_seq(seq):
    # handle upper case only
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def SlidingWindowREObj(seq):
    l = len(seq)
    seq = seq * 2
    sw = []
    for i in range(l):
        sw.append(seq[i:i+l])
        sw.append(rc_seq(seq[i:i+l]))
    return re.compile('(?:{})'.format('|'.join(sw)))

parser = argparse.ArgumentParser(description='Find telomeric sequences.')
parser.add_argument('fasta',
                    help='fasta')
parser.add_argument('-f', action='store_true',
                    help='force redo')
parser.add_argument('-r', type=str, metavar='repeat',
                    help='user provided sequence to match telomeric sequences (only return matched)')
parser.add_argument('--window', type=int, default=5000, metavar='',
                    help='size of terminal window (5000)')
parser.add_argument('--min_copy_number', type=int, default=5, metavar='',
                    help='minimal copy numbers (5)')
parser.add_argument('--min_repeat_len', type=int, default=4, metavar='',
                    help='minimal repeat length (4)')

args = parser.parse_args()
in_fa = args.fasta
repeat_seq = args.r
search_window = args.window
min_copy_number = args.min_copy_number
min_repeat_len = args.min_repeat_len
trf_opt = [2, 7, 7, 80, 10, 50, 500]

if repeat_seq:
    match_mode = True
    repeat_seq = SlidingWindowREObj(repeat_seq.upper())
else:
    match_mode = False

# parse fasta index
fa_idx = in_fa + '.fai'
if args.f or not os.path.isfile(fa_idx):
    print('Building index...', file=sys.stderr)
    cmd = ['samtools', 'faidx', in_fa]
    subprocess.run(cmd, check=True)

ctg_len = {}
with open(fa_idx) as f:
    for line in f:
        line = line.strip().split()
        ctg_len[line[0]] = int(line[1])

# run trf
trf_opt = [str(x) for x in trf_opt]
trf_out = '{}.{}.dat'.format(in_fa, '.'.join(trf_opt))
if args.f or not os.path.isfile(trf_out):
    print('Running trf...', file=sys.stderr)
    cmd = ['trf', in_fa] + trf_opt + ['-h', '-ngs']
    with open(trf_out, 'w') as f:
        subprocess.run(cmd, stdout=f, stderr=subprocess.DEVNULL, check=True)

# parse trf output
with open(trf_out) as f:
    FoundInSTART = set()
    FoundInEND = set()
    FoundInMID = set()
    for line in f:
        if line.startswith('@'):
            ctg = line[1:].strip().split()[0]
        else:
            line = line.strip()
            # skip blank line
            if len(line) == 0:
                continue
            line = line.split()
            # filter copy number and repeat length
            if float(line[3]) < min_copy_number or int(line[2]) < min_repeat_len:
                continue

            if match_mode and not repeat_seq.match(line[13]):
                continue

            # START
            if int(line[1]) < search_window:
                FoundInSTART.add(ctg)
                print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(ctg, ctg_len[ctg], 'START',
                                                              line[0], line[1], line[2], line[3], line[13]))
            # END
            elif ctg_len[ctg] - int(line[0]) < search_window:
                FoundInEND.add(ctg)
                print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(ctg, ctg_len[ctg], 'END',
                                                              line[0], line[1], line[2], line[3], line[13]))
            elif match_mode:
                FoundInMID.add(ctg)
                print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(ctg, ctg_len[ctg], 'MID',
                                                              line[0], line[1], line[2], line[3], line[13]))

if match_mode:
    FoundInBoth = FoundInSTART.intersection(FoundInEND)
    FoundInSingle = FoundInSTART.union(FoundInEND).difference(FoundInBoth)
    print()
    print('Found in both ends:  {}\n'
          'Found in single end: {}\n'
          'Found in interval:   {}'.format(len(FoundInBoth), len(FoundInSingle), len(FoundInMID)))
    print(FoundInBoth, FoundInSingle, FoundInMID, sep='\n')
