#!/usr/bin/env python3

import sys
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='wrap or unwrap fasta')
parser.add_argument('fasta', nargs='?', type=str, default=sys.stdin,
                    help='fasta file (default=stdin)')
parser.add_argument('-n', type=int, default=0,
                    help='wrap by N characters (default=0, unwrap)')
parser.add_argument('-i', action='store_true',
                    help='overwrite file in-place (use with caution!)')
args = parser.parse_args()

FIN = args.fasta
n = args.n
recs = [x for x in SeqIO.parse(FIN, 'fasta')]

FOUT = open(args.fasta, 'w') if args.i else sys.stdout
SeqIO.FastaIO.FastaWriter(FOUT, wrap=n).write_file(recs)
FOUT.close()
