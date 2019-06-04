#!/usr/bin/env python3
'''
support output .nex, including source info
'''

import sys
import argparse

def FastaParser(handle):
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

def parse_all_files(files, sep=None):
    seqs = {}
    for file in files:
        with open(file) as f:
            for title, seq in FastaParser(f):
                title = title.split(sep)[0]
                seqs.setdefault(title, []).append(seq)
    return seqs

# def print_wrapped(s, wrap=80, **kwargs):
#     if wrap > 0:
#         s = s.strip()
#         l = len(s)
#         st = 0
#         while st + wrap < l:
#             print(s[st:st+wrap], **kwargs)
#             st += wrap
#         else:
#             print(s[st:], **kwargs)
#     else:
#         print(s, **kwargs)

def parse_args():
    parser = argparse.ArgumentParser(description='concatenate fasta files by sequence IDs')
    parser.add_argument('fasta', nargs='+', type=str,
                        help='fasta files')
    # parser.add_argument('-n', type=int, default=80,
    #                     help='wrap by N characters (default=80, set 0 to unwrap)')
    parser.add_argument('--sep', type=str, default=None,
                        help='separator (default="\\s")')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    files = args.fasta
    sep = args.sep

    seqs = parse_all_files(files, sep)
    for k, v in seqs.items():
        if len(v) != len(files):
            print('{} does not exist in all fasta files. Skip.'.format(k), file=sys.stderr)
            continue
        print('>{}\n{}'.format(k, ''.join(v)))
