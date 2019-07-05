#!/usr/bin/env python3

'''
Discard sequences shorter than cutoff.
'''

import sys
import gzip
import io
import os
import argparse
import re

# def print_wrapped(s, wrap=80, **kwargs):
#     s = s.strip()
#     l = len(s)
#     st = 0
#     while st + wrap < l:
#         print(s[st:st+wrap], **kwargs)
#         st += wrap
#     else:
#         print(s[st:], **kwargs)

# def check_format(handle, beginner):
#     if next(handle).startswith(beginner):
#         handle.seek(0)
#         return
#     else:
#         raise ValueError('Invalid file format.')

def fasta_filter(fin, fout, cutoff):
    rec_tot = 0
    rec_flt = 0
    base_tot = 0
    base_flt = 0
    seq_length = 0
    seq_stored = []
    # first record
    for line in fin:
        if line.startswith(b'>'):
            rec_tot += 1
            header = line
            break
        else:
            raise ValueError('Invalid file format.')
    for line in fin:
        if line.startswith(b'>'):
            if seq_length >= cutoff:
                fout.write(header)
                for s in seq_stored:
                    fout.write(s)
            else:
                rec_flt += 1
                base_flt += seq_length
            base_tot += seq_length
            rec_tot += 1
            header = line
            seq_length = 0
            seq_stored = []
        elif line != b'\n':  # skip blank line
            seq_length += (len(line.strip()))  # \n not trimmed
            seq_stored.append(line)
    # last record
    base_tot += seq_length
    if seq_length >= cutoff:
        fout.write(header)
        for s in seq_stored:
            fout.write(s)
    else:
        rec_flt += 1
        base_flt += seq_length
    return rec_tot, rec_flt, base_tot, base_flt

def fastq_filter(fin, fout, cutoff):
    rec_tot = 0
    rec_flt = 0
    base_tot = 0
    base_flt = 0
    try:
        # first record
        for l1 in fin:
            if not l1.startswith(b'@'):
                raise ValueError('Invalid file format.')
            l2 = next(fin)
            rec_tot += 1
            base_tot += (len(l2) - 1)
            if len(l2) - 1 >= cutoff:  # \n not trimmed
                fout.write(l1)
                fout.write(l2)
                fout.write(next(fin))
                fout.write(next(fin))
            else:
                rec_flt += 1
                base_flt += (len(l2) - 1)
                next(fin)  # l3
                next(fin)  # l4
            break
        for l1 in fin:
            l2 = next(fin)
            rec_tot += 1
            base_tot += (len(l2) - 1)
            if len(l2) - 1 >= cutoff:  # \n not trimmed
                fout.write(l1)
                fout.write(l2)
                fout.write(next(fin))
                fout.write(next(fin))
            else:
                rec_flt += 1
                base_flt += (len(l2) - 1)
                next(fin)  # l3
                next(fin)  # l4
    except StopIteration:
        raise ValueError('Invalid file format.')
    return rec_tot, rec_flt, base_tot, base_flt

def parse_args():
    parser = argparse.ArgumentParser(description='Discard sequences shorter than cutoff.')
    parser.add_argument('fastx', type=str, nargs='?',
                        help='fasta(q), can be gzipped, leave blank to read from STDIN')
    parser.add_argument('-c', metavar='cutoff', type=int, default=1000,
                        help='length cutoff (1000)')
    parser.add_argument('-o', metavar='file', type=str, default=None,
                        help='output file (AUTO)')
    group1 = parser.add_mutually_exclusive_group(required=True)
    group1.add_argument('--fa', action='store_true',
                        help='fasta format')
    group1.add_argument('--fq', action='store_true',
                        help='fastq format')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    fastx = args.fastx
    foutname = args.o
    cutoff = args.c
    stdin = False if fastx else True
    gz = True if fastx and fastx.endswith('.gz') else False
    if args.fa:
        mode = 'fa'
    elif args.fq:
        mode = 'fq'

    # use binary stream
    if stdin:
        # fin = sys.stdin.buffer
        fin = os.fdopen(sys.stdin.fileno(), 'rb', closefd=False)
    elif gz:
        fin = io.BufferedReader(gzip.open(fastx))
    else:
        fin = open(fastx, 'rb')

    if foutname:
        fout = open(foutname, 'wb')
    elif stdin:
        fout = os.fdopen(sys.stdout.fileno(), 'wb', closefd=False)
    else:
        foutname = re.sub('.fa$|.fa.gz$|.fasta$|.fasta.gz$|.fq$|.fq.gz$|.fastq$|.fastq.gz$', '', fastx)
        if mode == 'fa':
            foutname = '{}.min{}.fa'.format(foutname, cutoff)
        elif mode == 'fq':
            foutname = '{}.min{}.fq'.format(foutname, cutoff)
        fout = open(foutname, 'wb')

    if mode == 'fa':
        print('Input is fasta. Cutoff is {}.'.format(cutoff), file=sys.stderr)
        rec_tot, rec_flt, base_tot, base_flt = fasta_filter(fin, fout, cutoff)
    elif mode == 'fq':
        print('Input is fastq. Cutoff is {}.'.format(cutoff), file=sys.stderr)
        rec_tot, rec_flt, base_tot, base_flt = fastq_filter(fin, fout, cutoff)

    # finish!
    fin.close()
    fout.close()
    print('Finish!', file=sys.stderr)
    if not stdin:
        print('{} generated.'.format(foutname), file=sys.stderr)
    print('{}/{} ({:.2%}) records are discarded.'.format(rec_flt, rec_tot, rec_flt / rec_tot), file=sys.stderr)
    print('{}/{} ({:.2%}) bases are discarded.'.format(base_flt, base_tot, base_flt / base_tot), file=sys.stderr)
