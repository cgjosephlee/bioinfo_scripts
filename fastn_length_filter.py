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

def fastq_get_read_length(fin):
    print('Reading through file...', file=sys.stderr)
    # rec_tot = 0
    # base_tot = 0
    lengths = []
    try:
        # first record
        for l1 in fin:
            if not l1.startswith(b'@'):
                raise ValueError('Invalid file format.')
            lengths.append(len(next(fin)) - 1)  # l2
            next(fin)  # l3
            next(fin)  # l4
            break
        for l1 in fin:
            lengths.append(len(next(fin)) - 1)  # l2
            next(fin)  # l3
            next(fin)  # l4
    except StopIteration:
        raise ValueError('Invalid file format.')
    return sorted(lengths, reverse=True)

def cal_new_cutoff(lengths, cutoff, targetBase):
    print('Total {} records, {} bases.'.format(len(lengths), sum(lengths)), file=sys.stderr)
    print('Target is {} bases.'.format(targetBase), file=sys.stderr)
    if targetBase > sum(lengths):
        print('You probably need to sequence more.', file=sys.stderr)
        sys.exit(1)
    baseSum = 0
    for i in lengths:
        if baseSum < targetBase:
            baseSum += i
        else:
            break
    if i > cutoff:
        print('New cutoff is {}, yielding {} bases.'.format(i, baseSum), file=sys.stderr)
        return i
    else:
        print('New cutoff is {} < {}, use {} still.'.format(i, cutoff, cutoff), file=sys.stderr)
        return cutoff

def fasta_filter(fin, fout, cutoff):
    rec_tot = 0
    rec_flt = 0
    base_tot = 0
    base_flt = 0
    lengths = []

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
                lengths.append(seq_length)
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
        lengths.append(seq_length)
        for s in seq_stored:
            fout.write(s)
    else:
        rec_flt += 1
        base_flt += seq_length
    return rec_tot, rec_flt, base_tot, base_flt, sorted(lengths, reverse=True)

def fastq_filter(fin, fout, cutoff):
    rec_tot = 0
    rec_flt = 0
    base_tot = 0
    base_flt = 0
    lengths = []
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
                lengths.append(len(l2) - 1)
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
                lengths.append(len(l2) - 1)
            else:
                rec_flt += 1
                base_flt += (len(l2) - 1)
                next(fin)  # l3
                next(fin)  # l4
    except StopIteration:
        raise ValueError('Invalid file format.')
    return rec_tot, rec_flt, base_tot, base_flt, sorted(lengths, reverse=True)

def parse_args():
    parser = argparse.ArgumentParser(description='Discard sequences shorter than cutoff.')
    parser.add_argument('fastn', type=str,
                        help='fasta(q), can be gzipped, \'-\' to read from STDIN')
    parser.add_argument('-c', metavar='cutoff', type=int, default=1000,
                        help='length cutoff (1000)')
    parser.add_argument('-o', metavar='file', type=str, default=None,
                        help='output file (AUTO)')
    parser.add_argument('-t', metavar='base', type=int, default=None,
                        help=argparse.SUPPRESS)
    # 'descending read lengths and keep this many bases only (require input from file)'
    group1 = parser.add_mutually_exclusive_group(required=True)
    group1.add_argument('--fa', action='store_true',
                        help='fasta format')
    group1.add_argument('--fq', action='store_true',
                        help='fastq format')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    FIN = args.fastn
    foutname = args.o
    cutoff = args.c
    targetBase = args.t
    is_stdin = True if FIN == '-' else False
    gz = True if FIN and FIN.endswith('.gz') else False
    if args.fa:
        mode = 'fa'
    elif args.fq:
        mode = 'fq'

    # use binary stream
    if is_stdin:
        if targetBase:
            raise IOError('Target base method requires input from file.')
        print('Reading from STDIN.', file=sys.stderr)
        # handle = sys.stdin.buffer
        handle = os.fdopen(sys.stdin.fileno(), 'rb', closefd=False)
    elif gz:
        handle = io.BufferedReader(gzip.open(FIN))
    else:
        handle = open(FIN, 'rb')

    if mode == 'fa':
        print('Input is fasta. Cutoff is {}.'.format(cutoff), file=sys.stderr)
    elif mode == 'fq':
        print('Input is fastq. Cutoff is {}.'.format(cutoff), file=sys.stderr)

    if targetBase:
        if mode == 'fa':
            # todo
            sys.exit(99)
        elif mode == 'fq':
            cutoff = cal_new_cutoff(fastq_get_read_length(handle), cutoff, targetBase)
        handle.seek(0)

    if foutname:
        fout = open(foutname, 'wb')
    elif is_stdin:
        fout = os.fdopen(sys.stdout.fileno(), 'wb', closefd=False)
    else:
        foutname = re.sub('.fa$|.fa.gz$|.fasta$|.fasta.gz$|.fq$|.fq.gz$|.fastq$|.fastq.gz$', '', FIN)
        if mode == 'fa':
            foutname = '{}.min{}.fa'.format(foutname, cutoff)
        elif mode == 'fq':
            foutname = '{}.min{}.fq'.format(foutname, cutoff)
        fout = open(foutname, 'wb')

    if mode == 'fa':
        rec_tot, rec_flt, base_tot, base_flt, lengths = fasta_filter(handle, fout, cutoff)
    elif mode == 'fq':
        rec_tot, rec_flt, base_tot, base_flt, lengths = fastq_filter(handle, fout, cutoff)
    handle.close()
    fout.close()

    # left reads
    total_len = sum(lengths)
    total_len_50 = total_len * 0.5
    total_len_90 = total_len * 0.9
    read_N50 = None
    read_N90 = None
    accu_len = 0

    for n, v in enumerate(lengths, 1):
        accu_len += v
        if not read_N50 and accu_len > total_len_50:
            read_N50 = (n, v)
        if not read_N90 and accu_len > total_len_90:
            read_N90 = (n, v)
            break

    print('Finish!', file=sys.stderr)
    if not is_stdin:
        print('{} generated.'.format(foutname), file=sys.stderr)
    print('{}/{} ({:.2%}) records are discarded.'.format(rec_flt, rec_tot, rec_flt / rec_tot), file=sys.stderr)
    print('{}/{} ({:.2%}) bases are discarded.'.format(base_flt, base_tot, base_flt / base_tot), file=sys.stderr)
    print('# total_base seq_num mean max min N50 N90', file=sys.stderr)
    print(
        total_len,
        len(lengths),
        round(total_len / len(lengths), 1),
        lengths[0],
        lengths[-1],
        read_N50[1],
        read_N90[1],
        sep='\t',
        file=sys.stderr
    )
