#!/usr/bin/env python3
'''
Basic statistics with pandas.
'''
import sys
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Calculate basic statstics from a file.')
parser.add_argument('file', type=str, default='-',
                    help='file, "-" to read from stdin')
parser.add_argument('-d', metavar='str', type=str, default='\t',
                    help='delimiter (\\t)')
parser.add_argument('--header', action='store_true',
                    help='has header')
parser.add_argument('--cols', metavar='list', type=str,
                    help='target columns, 0-base index or colnames (with -h), comma separated (all)')
parser.add_argument('--skip', metavar='int', type=int,
                    help='skip top n rows')
args = parser.parse_args()

infile = args.file if args.file != '-' else sys.stdin
delim = args.d
has_header = 0 if args.header else None
skip_rows = args.skip
try:
    cols = [int(x) for x in args.cols.split(',')] if args.cols else None
except ValueError:
    cols = [x for x in args.cols.split(',')] if args.cols else None

df = pd.read_csv(infile,
                 sep=delim,
                 header=has_header,
                 usecols=cols,
                 skiprows=skip_rows)

print(df.describe().T)
