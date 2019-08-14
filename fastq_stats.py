#!/usr/bin/env python3

import sys
import io
import gzip
import argparse
from statistics import median

def check_format(handle, beginner):
    if next(handle).startswith(beginner):
        handle.seek(0)
        return
    else:
        raise ValueError('Invalid file format.')

parser = argparse.ArgumentParser(description='fastq statistics for long reads')
parser.add_argument('fastq',
                    help='fastq file')
parser.add_argument('cutoff', nargs='?', type=int,
                    help='minimun length cutoff (0)')
parser.add_argument('-p', action='store_true',
                    help='plot histogram (require matplotlib)')
parser.add_argument('--fa', action='store_true',
                    help='support fasta formatted read file')
args = parser.parse_args()

FQ = args.fastq
cutoff = args.cutoff if args.cutoff else 0

is_gzipped = True if FQ.endswith('.gz') else False
if is_gzipped:
    handle = io.BufferedReader(gzip.open(FQ))
else:
    handle = open(FQ, 'rb')

lengths = []
if args.fa:  # fasta
    beginner = b'>'
    check_format(handle, beginner)

    for line in handle:
        if line.startswith(beginner):
            lengths.append(0)
        else:
            lengths[-1] += len(line.strip())
else:  # fastq
    beginner = b'@'
    check_format(handle, beginner)

    for line in handle:
        lengths.append(len(next(handle).strip()))
        next(handle)
        next(handle)
handle.close()

if cutoff != 0:
    lengths = [x for x in lengths if x > cutoff]
lengths.sort(reverse=True)
total_len = sum(lengths)
total_len_50 = total_len * 0.5
total_len_90 = total_len * 0.9
read_N50 = None
read_N90 = None
accu_len = 0

for n, v in enumerate(lengths, 1):
    accu_len += v
    if not read_N50 and accu_len > total_len_50:
        read_N50 = (n, v)
    if not read_N90 and accu_len > total_len_90:
        read_N90 = (n, v)
        break

print('''\
input fastq file: {}

minimum length: {}
maximum length: {}
total seqs:     {}
total length:   {} ({:.2f} Gbp)
avg. length:    {:.1f}
median length:  {:.1f}
N90:            {}
N50:            {}\
'''.format(FQ,
           lengths[-1],
           lengths[0],
           len(lengths),
           total_len, total_len / 10 ** 9,
           total_len / len(lengths),
           median(lengths),
           read_N90[1],
           read_N50[1]
           ))

if args.p:
    try:
        import matplotlib as mpl
        import matplotlib.pyplot as plt
    except ImportError:
        raise ImportError('Please install matplotlib to plot histogram.')
    print('\nPlotting histogram...', file=sys.stderr)
    out_png = FQ + '.hist.png'
    plt.hist(lengths, weights=lengths, bins=10**2)
    plt.axvline(read_N50[1], color='b', linewidth=.5)
    plt.axvline(read_N90[1], color='b', linewidth=.5)
    ticks_y = mpl.ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x / 10**6))
    plt.gca().yaxis.set_major_formatter(ticks_y)
    ticks_x = mpl.ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x / 10**3))
    plt.gca().xaxis.set_major_formatter(ticks_x)
    plt.title(FQ)
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
