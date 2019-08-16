#!/usr/bin/env python3
import sys
import argparse
from collections import defaultdict
import statistics as stat
# import numpy as np

parser = argparse.ArgumentParser(description='Easily calculate basic statstics from specific columns of a file.')
parser.add_argument('file', type=str, default='-',
                    help='file, "-" to read from stdin')
parser.add_argument('-d', metavar='str', type=str, default='\t',
                    help='delimiter (\\t)')
parser.add_argument('-c', metavar='int', type=int, default=[1], nargs='+',
                    help='target columns, 1-base index, separated by space (1)')
parser.add_argument('-s', metavar='int', type=int, default=0,
                    help='skip top n lines (0)')
args = parser.parse_args()

infile = args.file
delim = args.d
cols = args.c
skip = args.s

if infile == '-':
    handle = sys.stdin
else:
    handle = open(infile)

total_rows = 0
arrays = defaultdict(list)

for line in handle:
    for _ in range(skip):
        next(handle)
    total_rows += 1
    line = line.strip().split(delim)
    for i in cols:
        try:
            arrays[i].append(float(line[i - 1]))
        except ValueError:
            # Inf, NaN, etc.
            continue

if infile != '-':
    handle.close()

for n, i in enumerate(cols):
    if n > 0:
        print('')

    array = arrays[i]
    print('''\
Column : {}
Count  : {}
Invalid: {}
Max    : {}
Min    : {}
Mean   : {:.4f}
Stdev  : {:.4f}
Median : {}\
'''.format(
        i,
        len(array),
        total_rows - len(array),
        max(array),
        min(array),
        stat.mean(array),
        stat.stdev(array),
        stat.median(array)
    ))
