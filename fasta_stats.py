#!/usr/bin/env python3

import gzip
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='fasta statistics')
    parser.add_argument('fasta',
                        help='fasta file')
    parser.add_argument('cutoff', nargs='?', type=int,
                        help='minimun length cutoff (0)')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    FIN = args.fasta
    FIN = gzip.open(FIN) if FIN.endswith('.gz') else open(FIN, 'rb')
    cutoff = args.cutoff if args.cutoff else 0

    lengths = {}
    total_len = 0
    total_N = 0
    total_seqs = 0
    for line in FIN:
        line = line.strip()
        if line.startswith(b'>'):
            title = line[1:].decode()
            lengths[title] = 0
            total_seqs += 1
        else:
            lengths[title] += len(line)
            total_len += len(line)
            total_N += (line.count(b'N') + line.count(b'n'))

    lengths = sorted(lengths.items(), key=lambda x: x[1], reverse=True)
    total_len_50 = total_len * 0.5
    total_len_90 = total_len * 0.9
    ctg_N50 = None
    ctg_N90 = None
    accu_len = 0
    above_len = [1e2, 1e3, 1e4, 1e5, 1e6]
    above_counts = [0] * len(above_len)

    for n, (k, v) in enumerate(lengths, 1):
        accu_len += v
        if not ctg_N50 and accu_len > total_len_50:
            ctg_N50 = (n, k, v)
        if not ctg_N90 and accu_len > total_len_90:
            ctg_N90 = (n, k, v)
        for i in range(len(above_len)):
            if v > above_len[i]:
                above_counts[i] += 1
    FIN.close()

    print('''\
input fasta file: {}
length cutoff: {}

minimum length: {} ({})
maximum length: {} ({})
2nd long seq:   {} ({})
3rd long seq:   {} ({})
total seqs:   {}
total length: {}
avg. length:  {:.3f}
number of Ns: {} ({:.3%})
L90: {} ({})
N90: {}
L50: {} ({})
N50: {}
'''.format(args.fasta,
           cutoff,
           lengths[-1][1], lengths[-1][0],
           lengths[0][1], lengths[0][0],
           lengths[1][1], lengths[1][0],
           lengths[2][1], lengths[2][0],
           total_seqs,
           total_len,
           total_len / total_seqs,
           total_N, total_N / total_len,
           ctg_N90[0], ctg_N90[1],
           ctg_N90[2],
           ctg_N50[0], ctg_N50[1],
           ctg_N50[2]
           ))

    print('''\
seq counts:
> 100 bp:  {}
>   1 Kbp: {}
>  10 Kbp: {}
> 100 Kbp: {}
>   1 Mbp: {}\
'''.format(above_counts[0],
           above_counts[1],
           above_counts[2],
           above_counts[3],
           above_counts[4]
           ))
