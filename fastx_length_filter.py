#!/usr/bin/env python3

'''
Discard sequences shorter than cutoff.

TODO:
    compressed output? gzip vs pigz?
'''

import sys
import gzip
import io
import argparse
import re

def print_wrapped(s, wrap=80, **kwargs):
    s = s.strip()
    l = len(s)
    st = 0
    while st + wrap < l:
        print(s[st:st+wrap], **kwargs)
        st += wrap
    else:
        print(s[st:], **kwargs)

def fasta_filter(fin, fout, cutoff, wrap=80):
    rec_tot = 0
    rec_flt = 0
    base_tot = 0
    base_flt = 0
    s = None
    try:
        for line in fin:
            line = line.strip()
            if line.startswith('>'):
                if s and len(s) >= cutoff:
                    base_tot += len(s)
                    print(h, file=fout)
                    print_wrapped(s, file=fout, wrap=wrap)
                elif s:
                    rec_flt += 1
                    base_tot += len(s)
                    base_flt += len(s)
                rec_tot += 1
                h = line
                s = ''
            else:
                s = s + line
        # last record
        base_tot += len(s)
        if len(s) >= cutoff:
            print(h, file=fout)
            print_wrapped(s, file=fout, wrap=wrap)
        else:
            rec_flt += 1
            base_flt += len(s)
    except Exception:
        print('Maybe not fasta?', file=sys.stderr)
        raise
    return rec_tot, rec_flt, base_tot, base_flt

def fastq_filter(fin, fout, cutoff):
    rec_tot = 0
    rec_flt = 0
    base_tot = 0
    base_flt = 0
    try:
        for line in fin:
            l1 = line
            l2 = next(fin)
            rec_tot += 1
            base_tot += (len(l2) - 1)
            if len(l2) - 1 >= cutoff:  # \n not trimmed
                print(l1, end='', file=fout)
                print(l2, end='', file=fout)
                print(next(fin), end='', file=fout)  # l3
                print(next(fin), end='', file=fout)  # l4
            else:
                rec_flt += 1
                base_flt += (len(l2) - 1)
                next(fin)  # l3
                next(fin)  # l4
    except Exception:
        print('Maybe not fastq?', file=sys.stderr)
        raise
    return rec_tot, rec_flt, base_tot, base_flt

def parse_args():
    parser = argparse.ArgumentParser(description='Discard sequences shorter than cutoff.')
    parser.add_argument('fastx', type=str, nargs='?',
                        help='fasta or fastq, can be gzipped')
    parser.add_argument('-c', metavar='cutoff', type=int, default=1000,
                        help='length cutoff (1000)')
    parser.add_argument('-o', metavar='file', type=str, default=None,
                        help='output file ({base}.min{cutoff}.fx)')
    group1 = parser.add_mutually_exclusive_group()
    group1.add_argument('--fa', action='store_true',
                        help='fasta from STDIN (default when no input file)')
    group1.add_argument('--fq', action='store_true',
                        help='fastq from STDIN')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    fastx = args.fastx
    cutoff = args.c
    gz = False
    stdin = True
    if fastx:
        stdin = False
        if fastx.endswith(('.fa', '.fasta')):
            mode = 'fa'
        elif fastx.endswith(('.fa.gz', '.fasta.gz')):
            mode = 'fa'
            gz = True
        elif fastx.endswith(('.fq', '.fastq')):
            mode = 'fq'
        elif fastx.endswith(('.fq.gz', '.fastq.gz')):
            mode = 'fq'
            gz = True
        else:
            raise IOError('Unknown input file format ({})'.format(fastx))
    else:
        if args.fq:
            mode = 'fq'
        else:
            mode = 'fa'

    if args.o:
        fout = open(args.o, 'w')
    elif stdin:
        fout = sys.stdout
    else:
        fout = re.sub('.fa$|.fa.gz$|.fasta$|.fasta.gz$|.fq$|.fq.gz$|.fastq$|.fastq.gz$', '', fastx)
        if mode == 'fa':
            fout = '{}.min{}.fa'.format(fout, cutoff)
        elif mode == 'fq':
            fout = '{}.min{}.fq'.format(fout, cutoff)
        fout = open(fout, 'w')

    if stdin:
        fin = sys.stdin
    elif gz:
        fin = io.TextIOWrapper(io.BufferedReader(gzip.open(fastx)))
    else:
        fin = open(fastx)

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
    print('{}/{} ({:.2%}) records are discarded.'.format(rec_flt, rec_tot, rec_flt / rec_tot), file=sys.stderr)
    print('{}/{} ({:.2%}) bases are discarded.'.format(base_flt, base_tot, base_flt / base_tot), file=sys.stderr)
