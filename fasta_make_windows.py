#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description='make windows')
parser.add_argument('faidx',
                    help='fasta index file')
parser.add_argument('-w', type=int,
                    help='window size')
parser.add_argument('-o', type=int, default=0,
                    help='window overlaps')
parser.add_argument('-m', type=int, default=0,
                    help='min window size, region is combined with previous one if shorter than this number')
parser.add_argument('-b', choices=[0, 1], type=int, default=0,
                    help='0-base or 1-base')
parser.add_argument('-f', choices=['tab', 'samtools'], default='tab',
                    help='format')
args = parser.parse_args()

w = args.w
o = args.o
s = args.w - args.o
m = args.m
assert w > o
assert w > m

handle = open(args.faidx)
for line in handle:
    line = line.strip().split()
    ID = line[0]
    LEN = int(line[1])

    starts = []
    i = 0
    while i < LEN:
        # print(i)
        starts.append(i)
        i += s

    if m and LEN - starts[-1] < m:
        starts.pop()

    for st in starts[:-1]:
        ed = st + w
        if args.b == 0:
            if args.f == 'tab':
                print(f'{ID}\t{st}\t{ed}')
            elif args.f == 'samtools':
                print(f'{ID}:{st}-{ed}')
        elif args.b == 1:
            if args.f == 'tab':
                print(f'{ID}\t{st+1}\t{ed}')
            elif args.f == 'samtools':
                print(f'{ID}:{st+1}-{ed}')
    st = starts[-1]
    ed = LEN
    if args.b == 0:
        if args.f == 'tab':
            print(f'{ID}\t{st}\t{ed}')
        elif args.f == 'samtools':
            print(f'{ID}:{st}-{ed}')
    elif args.b == 1:
        if args.f == 'tab':
            print(f'{ID}\t{st+1}\t{ed}')
        elif args.f == 'samtools':
            print(f'{ID}:{st+1}-{ed}')

handle.close()
