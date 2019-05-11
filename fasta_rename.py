#!/usr/bin/env python3

import sys
import gzip
import argparse

def fa2dict_sort(handle):
    data = {}
    for line in handle:
        line = line.strip()
        if line.startswith('>'):
            h = line[1:]
            data[h] = {'seq': [], 'length': 0}
        else:
            data[h]['seq'].append(line)
            data[h]['length'] += len(line)
    data = sorted(data.items(), key=lambda x: x[1]['length'], reverse=True)
    return data

parser = argparse.ArgumentParser(description='rename fasta IDs serially')
parser.add_argument('fasta', type=str,
                    help='fasta file ("-" to read stdin)')
parser.add_argument('-s', action='store_true',
                    help='sort sequences by length')
parser.add_argument('-p', type=str, default='scaff',
                    help='ID prefix (scaff)')
parser.add_argument('-n', type=int, default=4,
                    help='digits of serial number (4)')
args = parser.parse_args()
prefix = args.p
digits = args.n

if args.fasta == '-':
    FIN = sys.stdin
elif args.fasta.endswith('.gz'):
    FIN = gzip.open(args.fasta, 'rt')
else:
    FIN = open(args.fasta)

n = 1
if not args.s:
    for line in FIN:
        if line.startswith('>'):
            print('>{}{:0>{}d}'.format(prefix, n, digits))
            n += 1
        else:
            print(line.strip())
else:
    FA = fa2dict_sort(FIN)
    for k, v in FA:
        print('>{}{:0>{}d}'.format(prefix, n, digits))
        print('{}'.format('\n'.join(v['seq'])))
        n += 1

if FIN != sys.stdin:
    FIN.close()
