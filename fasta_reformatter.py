#!/usr/bin/env python3

import sys
import os
import argparse

def FastaParser(handle):
    """
    Iterate over Fasta records as string tuples.
    Took from Biopython FastaIO.
    """
    # Skip any text before the first record (e.g. blank lines, comments)
    for line in handle:
        if line[0] == '>':
            title = line[1:].rstrip()
            break
    else:   # no break encountered
        return  # Premature end of file, or just empty?

    lines = []
    for line in handle:
        if line[0] == '>':
            yield title, ''.join(lines).replace(" ", "").replace("\r", "")
            lines = []
            title = line[1:].rstrip()
            continue
        lines.append(line.rstrip())
    yield title, ''.join(lines).replace(" ", "").replace("\r", "")

def print_wrapped(s, wrap=80, **kwargs):
    if wrap > 0:
        s = s.strip()
        l = len(s)
        st = 0
        while st + wrap < l:
            print(s[st:st+wrap], **kwargs)
            st += wrap
        else:
            print(s[st:], **kwargs)
    else:
        print(s, **kwargs)

def parse_args():
    parser = argparse.ArgumentParser(description='wrap or unwrap fasta (to stdout)')
    parser.add_argument('fasta', nargs='?', type=str, default=sys.stdin,
                        help='fasta file (default=stdin)')
    parser.add_argument('-n', type=int, default=80,
                        help='wrap by N characters (default=80, set 0 to unwrap)')
    parser.add_argument('-i', action='store_true',
                        help='overwrite file in-place (use with caution!)')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    n = args.n
    FIN = args.fasta

    if FIN != sys.stdin and args.i:
        os.rename(FIN, FIN + '.bak')
        print('backup file {}.bak created'.format(FIN), file=sys.stderr)
        FOUT = open(FIN, 'w')
        FIN = FIN + '.bak'
    else:
        FOUT = sys.stdout

    if FIN != sys.stdin:
        FIN = open(FIN)

    for rec in FastaParser(FIN):
        print('>' + rec[0], file=FOUT)
        print_wrapped(rec[1], wrap=n, file=FOUT)

    FIN.close()
    FOUT.close()
