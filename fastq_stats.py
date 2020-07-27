#!/usr/bin/env python3

import sys
import io
import gzip
import argparse
from statistics import median

# try new method, but slower QQ
# def get_read_length(handle):
#     lengths = []
#     readLen = 0
#     last = None
#     while True:
#         if not last:
#             for line in handle:
#                 if line[0] in b'>@':
#                     last = line
#                     break
#         if not last: break
#         last = None
#         for line in handle:  read sequence
#             if line.startswith((b'@', b'+', b'>')):
#                 last = line
#                 break
#             readLen += (len(line)-1)
#         if not last or not last.startswith(b'+'):  this is fasta
#             lengths.append(readLen)
#             readLen = 0
#             if not last: break
#         else:  this is fastq
#             try:
#                 next(handle)  ignore qual line
#             except StopIteration:  EOF
#                 print('what?')
#                 raise
#     return sorted(lengths, reverse=True)

def get_read_length(handle):
    lengths = []
    type = None
    for line in handle:
        if line.startswith(b'@'):
            lengths.append(len(next(handle)) - 1)
            next(handle)
            next(handle)
            type = 'fq'
            break
        elif line.startswith(b'>'):
            lengths.append(0)
            type = 'fa'
            break
    if type == 'fq':
        try:
            for line in handle:  # l1
                lengths.append(len(next(handle)) - 1)  # l2
                next(handle)  # l3
                next(handle)  # l4
        except StopIteration:
            print(line)
            raise StopIteration('Broken file?')
    elif type == 'fa':
        for line in handle:
            if line.startswith(b'>'):
                lengths.append(0)
            else:
                lengths[-1] += len(line.rstrip())
    else:
        raise OSError('Empty file?')
    return sorted(lengths, reverse=True)

def parse_args():
    parser = argparse.ArgumentParser(description='fastq statistics for long reads')
    parser.add_argument('fastq',
                        help='fastq file, \'-\' to read from STDIN')
    parser.add_argument('-c', metavar='cutoff', type=int, default=0,
                        help='minimun length cutoff (0)')
    parser.add_argument('-p', action='store_true',
                        help='plot histogram (require matplotlib)')
    parser.add_argument('-s', action='store_true',
                        help='print single line format')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    FIN = args.fastq
    cutoff = args.c

    is_stdin = True if FIN == '-' else False
    is_gzipped = True if FIN.endswith('.gz') else False
    if is_stdin:
        handle = sys.stdin.buffer
    elif is_gzipped:
        handle = io.BufferedReader(gzip.open(FIN))
    else:
        handle = open(FIN, 'rb')

    lengths = get_read_length(handle)
    handle.close()

    if cutoff != 0:
        lengths = [x for x in lengths if x >= cutoff]

    # calculate N10-N90
    total_len = sum(lengths)
    partitions = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9]
    thresholds = [total_len * x for x in partitions]
    read_Nxx = []
    accu_len = 0
    n = 0
    for current_threshold in thresholds:
        while accu_len < current_threshold:
            n += 1  # count
            accu_len += lengths[n-1]  # index
        read_Nxx.append((n, lengths[n-1], accu_len))

    read_N50 = read_Nxx[8]
    read_N90 = read_Nxx[16]

    if args.s:
        # total_base seq_num mean max min L50 L90
        print(
            total_len,
            len(lengths),
            round(total_len / len(lengths), 1),
            lengths[0],
            lengths[-1],
            read_N50[1],
            read_N90[1],
            sep='\t'
        )
    else:
        print('''\
input fastq file: {}

minimum length: {}
maximum length: {}
total seqs:     {}
total length:   {} ({:.2f} Gbp)
avg. length:    {:.1f}
median length:  {:.1f}
N50:            {}
N90:            {} '''.format(FIN, lengths[-1], lengths[0], len(lengths), total_len, total_len / 10 ** 9,
           total_len / len(lengths),
           median(lengths),
           read_N50[1],
           read_N90[1]
           ))
        print('# total_base seq_num mean max min N50 N90')
        print(
            total_len,
            len(lengths),
            round(total_len / len(lengths), 1),
            lengths[0],
            lengths[-1],
            read_N50[1],
            read_N90[1],
            sep='\t'
        )
        print('# Nxx count Nxx_len accu_len')
        for n, p in enumerate(partitions):
            print('N{}\t{}\t{}\t{}'.format(int(p*100), read_Nxx[n][0], read_Nxx[n][1], read_Nxx[n][2]))

    if args.p:
        try:
            import matplotlib as mpl
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError('Please install matplotlib to plot histogram.')
        print('\nPlotting histogram...', file=sys.stderr)
        out_png = FIN + '.hist.png'
        plt.hist(lengths, weights=lengths, bins=10**2)
        plt.axvline(read_N50[1], color='b', linewidth=.5)
        plt.axvline(read_N90[1], color='b', linewidth=.5)
        ticks_y = mpl.ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x / 10**6))
        plt.gca().yaxis.set_major_formatter(ticks_y)
        ticks_x = mpl.ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x / 10**3))
        plt.gca().xaxis.set_major_formatter(ticks_x)
        plt.title(FIN)
        plt.ylabel('Total bases (Mbp)')
        plt.xlabel('Read length (Kbp)')
        t = '''\
        Total yield = {:.2f} Gbp
        Max. length = {} bp
        N50 = {} bp
        N90 = {} bp'''.format(total_len / 10**9, lengths[0], read_N50[1], read_N90[1])
        plt.text(0.9, 0.9, t,
                 transform=plt.gca().transAxes,
                 verticalalignment='top',
                 horizontalalignment='right', ma='left')

        plt.savefig(out_png, dpi=300)
        print('{} generated.'.format(out_png), file=sys.stderr)
