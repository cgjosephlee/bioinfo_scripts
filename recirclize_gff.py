#!/usr/bin/env python3

import sys
import argparse


def parse_arg():
    parser = argparse.ArgumentParser(description='restate a gff from a specific coordinate')
    parser.add_argument('gff', type=str, help='gff file')
    parser.add_argument('-n', metavar='POS', type=int, default=1,
                        help='start from this position')
    parser.add_argument('-r', action='store_true',
                        help='reverse complement')
    parser.add_argument('-l', metavar="LEN", type=int, default=float('Inf'),  # required=True,
                        help='sequence full length (default: Inf)')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args()

# main
args = parse_arg()

with open(args.gff, "r") as f:
    for line in f.readlines():
        if not line.startswith('#'):
            line = line.strip()
            line = line.split('\t')
            if not args.r:
                start = int(line[3]) - args.n + 1
                len = int(line[4]) - int(line[3]) + 1
                if start < 0:
                    start = start + args.l
                end = start + len - 1
                print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(line[0], line[1], line[2], start, end, line[5], line[6], line[7], line[8]))
            elif args.r and args.l:
                start = args.n - int(line[4]) + 1
                len = int(line[4]) - int(line[3]) + 1
                if start < 0:
                    start = start + args.l
                end = start + len - 1
                if line[6] == '-':
                    orientation = '+'
                else:
                    orientation = '-'
                print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(line[0], line[1], line[2], start, end, line[5], orientation, line[7], line[8]))
