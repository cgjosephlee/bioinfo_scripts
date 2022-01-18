#!/usr/bin/env python3

import gzip
import io
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='fasta statistics')
    parser.add_argument('fasta',
                        help='fasta file')
    parser.add_argument('-l', metavar='len', type=int, default=0,
                        help='minimun sequence length (0)')
    parser.add_argument('-s', action='store_true',
                        help='print single line format')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    FIN = args.fasta
    FIN = io.BufferedReader(gzip.open(FIN)) if FIN.endswith('.gz') else open(FIN, 'rb')
    cutoff = args.l

    lengths = {}
    Ns = {}
    for line in FIN:
        line = line.strip()
        if line.startswith(b'>'):
            title = line[1:].decode()
            lengths[title] = 0
            Ns[title] = 0
        else:
            lengths[title] += len(line)
            Ns[title] += (line.count(b'N') + line.count(b'n'))
    FIN.close()

    lengths = sorted(filter(lambda x: x[1] > cutoff, lengths.items()), key=lambda x: x[1], reverse=True)
    total_len = sum([x[1] for x in lengths])
    titles_pass = [x[0] for x in lengths]
    total_N = sum([x[1] for x in filter(lambda x: x[0] in titles_pass, Ns.items())])
    total_seqs = len(lengths)

    auN = sum([x[1]**2 for x in lengths]) / total_len
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

    if args.s:
        # total_base seq_num mean(kb) max(kb) N50(kb) L50 N90(kb) L90 N
        print(
            total_len,
            total_seqs,
            '{:.1f}'.format(total_len / total_seqs / 1e3),
            '{:.1f}'.format(lengths[0][1] / 1e3),
            '{:.1f}'.format(ctg_N50[2] / 1e3),
            ctg_N50[0],
            '{:.1f}'.format(ctg_N90[2] / 1e3),
            ctg_N90[0],
            total_N,
            sep='\t'
        )
    elif total_seqs < 3:
        print('''\
input fasta file: {}
length cutoff: {}

minimum length: {} ({})
maximum length: {} ({})
total seqs:   {}
total length: {}
avg. length:  {:.3f}
number of Ns: {} ({:.3%})
'''.format(args.fasta,
           cutoff,
           lengths[-1][1], lengths[-1][0],
           lengths[0][1], lengths[0][0],
           total_seqs,
           total_len,
           total_len / total_seqs,
           total_N, total_N / total_len,
           ))
    else:
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
L50: {} ({})
N50: {}
L90: {} ({})
N90: {}
auN: {:.3f}
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
           ctg_N50[0], ctg_N50[1],
           ctg_N50[2],
           ctg_N90[0], ctg_N90[1],
           ctg_N90[2],
           auN
           ))

        print('''\
seq counts:
> 100 bp:  {}
>   1 Kbp: {}
>  10 Kbp: {}
> 100 Kbp: {}
>   1 Mbp: {}
'''.format(above_counts[0],
           above_counts[1],
           above_counts[2],
           above_counts[3],
           above_counts[4]
           ))

        print('# total_base seq_num mean(kb) max(kb) N50(kb) L50 N90(kb) L90 N')
        print(
            total_len,
            total_seqs,
            '{:.1f}'.format(total_len / total_seqs / 1e3),
            '{:.1f}'.format(lengths[0][1] / 1e3),
            '{:.1f}'.format(ctg_N50[2] / 1e3),
            ctg_N50[0],
            '{:.1f}'.format(ctg_N90[2] / 1e3),
            ctg_N90[0],
            total_N,
            sep='\t'
        )
