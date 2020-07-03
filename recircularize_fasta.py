#!/usr/bin/env python3

import argparse
from Bio import SeqIO

def parse_arg():
    parser = argparse.ArgumentParser(description='recircularize a circular sequence to a specific position (single sequence only)')
    parser.add_argument('fasta', type=str, help='fasta file')
    parser.add_argument('-n', metavar='POS', type=int, default=1,
                        help='start position (1-base)')
    parser.add_argument('-r', action='store_true',
                        help='reverse complement')
    parser.add_argument('-t', metavar='TITLE', type=str,
                        help='provide new title')
    return parser.parse_args()

# main
args = parse_arg()

with open(args.fasta, 'r') as f:
    FA = SeqIO.read(f, 'fasta')

if not args.r:
    title = FA.id + ' [start:%d]' % (args.n)
    seq = FA.seq[args.n - 1:] + FA.seq[:args.n - 1]
else:
    title = FA.id + ' [start:%d,reverse_complement]' % (args.n)
    seq = FA.seq[args.n:] + FA.seq[:args.n]
    seq = seq.reverse_complement()

if args.t:
    print('>' + args.t)
else:
    print('>' + title)
print(seq)
