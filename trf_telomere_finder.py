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
    assert isinstance(seq, str)
    seq = seq.upper()
    l = len(seq)
    seq = seq * 2
    sw = []
    for i in range(l):
        sw.append(seq[i:i+l] + '$')
        sw.append(rc_seq(seq[i:i+l]) + '$')
    return re.compile('(?:{})'.format('|'.join(sw)))

def parse_faidx(handle):
    ctg_len = {}
    for line in handle:
        line = line.strip().split()
        ctg_len[line[0]] = int(line[1])
    return ctg_len

'''
cols 1+2:  Indices of the repeat relative to the start of the sequence
col 3:     Period size of the repeat
col 4:     Number of copies aligned with the consensus pattern
col 5:     Size of consensus pattern (may differ slightly from the period size)
col 6:     Percent of matches between adjacent copies overall
col 7:     Percent of indels between adjacent copies overall
col 8:     Alignment score
cols 9-12: Percent composition for each of the four nucleotides
col 13:    Entropy measure based on percent composition
col 14:    Consensus sequence
col 15:    Repeat sequence
'''
def parse_trf_output(handle, ctg_len,
                     min_copy_number=5, min_repeat_len=4, search_window=5000,
                     repeat_seq=None, silent=False):
    FoundInSTART = set()
    FoundInEND = set()
    FoundInMID = set()
    if repeat_seq:
        repeat_seq = SlidingWindowREObj(repeat_seq)
    for line in handle:
        if line.startswith('@'):
            ctg = line[1:].strip().split()[0]
        else:
            line = line.strip()
            # skip blank line
            if len(line) == 0:
                continue
            line = line.split()
            # filter copy number and repeat consensus size
            if float(line[3]) < min_copy_number or int(line[2]) < min_repeat_len:
                continue

            if repeat_seq and not repeat_seq.match(line[13]):
                continue

            line[0] = int(line[0])
            line[1] = int(line[1])
            # cols: ctg ctg_len label pos_start pos_end length repeat_size repeat_copies repeat_seq
            # START
            if line[1] < search_window:
                FoundInSTART.add(ctg)
                if not silent:
                    print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                        ctg, ctg_len[ctg], 'START',
                        line[0], line[1], line[1]-line[0]+1, line[2], line[3], line[13]))
            # END
            elif ctg_len[ctg] - line[0] < search_window:
                FoundInEND.add(ctg)
                if not silent:
                    print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                        ctg, ctg_len[ctg], 'END',
                        line[0], line[1], line[1]-line[0]+1, line[2], line[3], line[13]))
            # MID
            elif repeat_seq:
                FoundInMID.add(ctg)
                if not silent:
                    print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                        ctg, ctg_len[ctg], 'MID',
                        line[0], line[1], line[1]-line[0]+1, line[2], line[3], line[13]))
    return FoundInSTART, FoundInEND, FoundInMID

def parse_args():
    parser = argparse.ArgumentParser(description='Find telomeric sequences.')
    parser.add_argument('fasta',
                        help='fasta')
    parser.add_argument('-f', action='store_true',
                        help='force redo')
    parser.add_argument('-r', type=str, metavar='repeat', default=None,
                        help='user provided sequence to match telomeric sequences (only return matched)')
    parser.add_argument('--window', type=int, default=5000, metavar='',
                        help='size of terminal window (5000)')
    parser.add_argument('--min_copy_number', type=int, default=5, metavar='',
                        help='minimal copy numbers (5)')
    parser.add_argument('--min_repeat_len', type=int, default=4, metavar='',
                        help='minimal repeat length (4)')
    return parser.parse_args()

def main():
    args = parse_args()
    in_fa = args.fasta
    repeat_seq = args.r
    search_window = args.window
    min_copy_number = args.min_copy_number
    min_repeat_len = args.min_repeat_len
    trf_opt = [2, 7, 7, 80, 10, 50, 500]  # default

    if repeat_seq:
        match_mode = True
    else:
        match_mode = False

    # run samtools faidx
    fa_idx = in_fa + '.fai'
    if args.f or not os.path.isfile(fa_idx):
        print('Building index...', file=sys.stderr)
        cmd = ['samtools', 'faidx', in_fa]
        subprocess.run(cmd, check=True)

    # run trf
    trf_opt = [str(x) for x in trf_opt]
    trf_out = '{}.{}.dat'.format(in_fa, '.'.join(trf_opt))
    if args.f or not os.path.isfile(trf_out):
        print('Running trf...', file=sys.stderr)
        cmd = ['trf', in_fa] + trf_opt + ['-h', '-ngs']
        with open(trf_out, 'w') as f:
            subprocess.run(cmd, stdout=f, stderr=subprocess.DEVNULL, check=True)

    # parse fasta index
    with open(fa_idx) as f:
        ctg_len = parse_faidx(f)

    # parse trf output
    with open(trf_out) as f:
        FoundInSTART, FoundInEND, FoundInMID = parse_trf_output(f, ctg_len,
                                                                min_copy_number, min_repeat_len, search_window,
                                                                repeat_seq)

    if match_mode:
        FoundInBoth = FoundInSTART.intersection(FoundInEND)
        FoundInSingle = FoundInSTART.union(FoundInEND).difference(FoundInBoth)
        print()
        print('Found in both ends : {}\n'
              'Found in single end: {}\n'
              'Found in interval  : {}'.format(len(FoundInBoth), len(FoundInSingle), len(FoundInMID)))
        print(sorted(FoundInBoth))
        print(sorted(FoundInSingle))
        print(sorted(FoundInMID))

if __name__ == '__main__':
    main()
