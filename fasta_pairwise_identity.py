#!/usr/bin/env python3

'''
Calculate pairwise identity of all sequence combinations in the fasta.

NOTE:
- Muscle v3.8.31 has a bug in aligning long sequences. Upgrade to v3.8.341 or higher to fix this.
TODO:
- Support protein sequences.
'''

import sys
import os
import io
import argparse
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
import pandas as pd
import subprocess as sp
from threading import Thread
from queue import Queue
import time
import tempfile

def handle_args():
    class AddOptAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            values = values.split()
            setattr(namespace, self.dest, values)
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description='''\
Calculate pairwise sequence identity from alignment file, output in phylip
distance matrix format.

Use extarnal msa file:
muscle -in input.fasta -out alignment.fasta
{prog} -f f alignment.fasta

Use this script to perform all-vs-all pairwise alignment:
{prog} --aligner muscle -o pairwise_alignment.phy input.fasta
{prog} -f p pairwise_alignment.phy'''.format(prog=os.path.basename(sys.argv[0])))
    parser.add_argument('IN', metavar='fasta', help='fasta file')
    parser.add_argument('-f', metavar='FORMAT', type=str, choices=['f', 'p'], default='f',
                        help='input file format [f: pre-aligned fasta; p: multiple alignment phylip] (default: f)')
    parser.add_argument('-g', metavar='MODE', type=int, choices=[0, 1, 2], default=0,
                        help='gap manipulation [0: BLAST identity; 1: gap-excluded identity; 2: gap-compressed identity] (default: 0)')
    parser.add_argument('-o', metavar='FILE', type=argparse.FileType('w'), default=sys.stdout,
                        help='output file (default: stdout)')
    parser.add_argument('--aligner', metavar='PROG', choices=['muscle', 'mafft', 'needle'],
                        help='if this is given, the script performs pariwise alignment on input file and generate a multiple alignment file in relaxed-phylip format, which can be further proceeded by this script [needle; muscle (external); mafft (external)]')
    parser.add_argument('--thread', type=int, default=1,
                        help='threads for alignment')
    # argparse do not allow value starts with '-', a workaround is to assign equal sign `--opts='--var'`, or add a leadng space `--opts ' --var'`
    parser.add_argument('--opts', action=AddOptAction, default=[],
                        help='additional options sent to aligner (space separated list, an equal sign "--opts=" is required)')
    # parser.add_argument('--type', metavar='type', choices=['nuc', 'prot'], default='nuc',
    #                     help='type of input sequences [nuc, prot] (default: nuc)')
    return parser.parse_args()

def check_program(prog):
    def check_muscle():
        proc = sp.Popen(['muscle', '-version'], stdout=sp.PIPE, stderr=sp.PIPE)
        out, err = proc.communicate()
        if proc.returncode == 0:
            ver = out.decode().split()[1]
            print('MUSCLE {} is available.'.format(ver), file=sys.stderr)
        else:
            print('MUSCLE is not available in your path!', file=sys.stderr)
            sys.exit(1)

    def check_mafft():
        # mafft has no -h or -v argument, print help to stderr
        proc = sp.Popen(['mafft', '-h'], stdout=sp.PIPE, stderr=sp.PIPE)
        out, err = proc.communicate()
        m = re.search(r'MAFFT (v[\d.]+)', err.decode())
        if m:
            ver = m.group(1)
            print('MAFFT {} is available.'.format(ver), file=sys.stderr)
        else:
            print('MAFFT is not available in your path!', file=sys.stderr)
            sys.exit(1)

    if prog == 'muscle':
        check_muscle()
    elif prog == 'mafft':
        check_mafft()

def seq_identity(s1, s2, gap_mode=0):
    '''
    https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
    gap_mode=0, BLAST identity, matches over alignment length.
    gap_mode=1, gap-excluded identity, all gaps are ignored.
    gap_mode=2, gap-compressed identity, gaps with any length are counted as one difference.
    '''
    assert isinstance(s1, str)
    assert isinstance(s2, str)
    if len(s1) != len(s2):
        raise ValueError('Unequal length in sequences!')

    s1 = s1.upper()
    s2 = s2.upper()
    DNA_CHAR = set('ATCG')
    GAP_CHAR = '-'
    Match = 0
    Mismatch = 0
    Gap = 0
    GapCount = 0
    LastIsGap = False

    for b1, b2 in zip(s1, s2):
        if b1 in DNA_CHAR and b2 in DNA_CHAR:
            if b1 == b2:
                Match += 1
            elif b1 != b2:
                Mismatch += 1
            LastIsGap = False
        elif b1 == GAP_CHAR or b2 == GAP_CHAR:
            # ignore column containing only gap
            if b1 == b2:
                continue
            Gap += 1
            if not LastIsGap:
                GapCount += 1
            LastIsGap = True
        else:
            # ambiguous bases are counted as mismatch
            Mismatch += 1

    if gap_mode == 0:
        return Match / (Match + Mismatch + Gap)
    elif gap_mode == 1:
        return Match / (Match + Mismatch)
    elif gap_mode == 2:
        return Match / (Match + Mismatch + GapCount)

def pairwise_alignment(s1, s2, prog='muscle', opts=[]):
    '''
    s1 and s2 are SeqRrecord objects.
    Perform alignment and return list of aligned sequences.
    '''
    if prog == 'needle':
        # EMBOSS needle default: DNAfull matrix (5, -4), gap penality (-10, -0.5)
        alignment = pairwise2.align.globalms(s1.seq, s2.seq, 5, -4, -10, -0.5, one_alignment_only=True)[0]
        # blastn
        # alignment = pairwise2.align.globalms(s1.seq, s2.seq, 2, -3, -5, -2, one_alignment_only=True)[0]
        return [SeqRecord(alignment[0], id=s1.id), SeqRecord(alignment[1], id=s2.id)]
    elif prog == 'muscle':
        cmd = ['muscle', '-quiet'] + opts
        proc = sp.Popen(cmd, stdin=sp.PIPE, stdout=sp.PIPE, universal_newlines=True)
        proc_out, proc_err = proc.communicate(input='>{}\n{}\n>{}\n{}'.format(s1.id, s1.seq, s2.id, s2.seq))
        if proc.returncode == 0:
            return [x for x in SeqIO.parse(io.StringIO(proc_out), 'fasta')]
        else:
            raise sp.CalledProcessError(proc.returncode, cmd, proc_err)
    elif prog == 'mafft':
        cmd = ['mafft', '--quiet'] + opts
        # mafft don't support stdin, so write a temporary fasta
        tmp_fa = tempfile.mkstemp(suffix='.fa', text=True)[1]
        cmd.append(tmp_fa)
        with open(tmp_fa, 'w') as f:
            print('>{}\n{}\n>{}\n{}'.format(s1.id, s1.seq, s2.id, s2.seq), file=f)
        proc = sp.Popen(cmd, stdout=sp.PIPE, universal_newlines=True)
        proc_out, proc_err = proc.communicate()
        os.remove(tmp_fa)
        if proc.returncode == 0:
            return [x for x in SeqIO.parse(io.StringIO(proc_out), 'fasta')]
        else:
            raise sp.CalledProcessError(proc.returncode, cmd, proc_err)

def pairwise_alignment_through_fasta(fastas, aligner, opts, thread, handle):
    '''
    input list of SeqRecord objects
    output in phylip format
    '''
    def worker():
        while not queue.empty():
            n1, n2, n = queue.get()
            outs[n] = pairwise_alignment(fastas[n1], fastas[n2], aligner, opts)

    check_program(aligner)

    # put jobs in queue
    seq_num = len(fastas)
    total_jobs = 0
    queue = Queue()
    for n1 in range(seq_num):
        for n2 in range(n1 + 1, seq_num):
            queue.put([n1, n2, total_jobs])
            total_jobs += 1
    outs = [None] * total_jobs
    # multi-threading
    threads = [Thread(target=worker) for _ in range(thread)]
    print('Open {} threads, running {} jobs...'.format(thread, total_jobs), file=sys.stderr)
    for th in threads:
        th.start()
    for th in threads:
        # th.join()
        while th.is_alive():
            print('{} jobs are queueing...'.format(queue.qsize()), end='\r', file=sys.stderr)
            time.sleep(0.1)
    print('', file=sys.stderr)

    # printing alignments in phylip format
    for s1, s2 in outs:
        maxLen = len(s2.id) if len(s2.id) > len(s1.id) else len(s1.id)
        print('2 {}'.format(len(s1.seq)), file=handle)
        print('{: <{}} {}'.format(s1.id, maxLen, s1.seq), file=handle)
        print('{: <{}} {}'.format(s2.id, maxLen, s2.seq), file=handle)
    print('Finish!', file=sys.stderr)

def df_to_phylip_matrix(df, handle):
    df.sort_index(axis=0, inplace=True)  # row
    df.sort_index(axis=1, inplace=True)  # column
    # right padding for index
    df.reset_index(inplace=True)
    maxLen = max(len(x) for x in df['index'])
    df['index'] = df['index'].apply(lambda x: x + ' ' * (maxLen - len(x)))
    print(df.shape[0], file=handle)  # nrows
    for row in df.itertuples(index=False):
        print(' '.join(row), file=handle)

def RelaxedPhylipParser(handle):
    '''
    3 lines per unit
    '''
    for line in handle:
        yield [next(handle).strip().split(), next(handle).strip().split()]

def main():
    args = handle_args()

    if not args.aligner:
        df = pd.DataFrame()
        if args.f == 'f':
            FA = [x for x in SeqIO.parse(args.IN, 'fasta')]
            seq_num = len(FA)
            for n1 in range(seq_num):
                for n2 in range(n1 + 1, seq_num):
                    df.at[FA[n1].id, FA[n2].id] = df.at[FA[n2].id, FA[n1].id] = '{:.8f}'.format(seq_identity(str(FA[n1].seq), str(FA[n2].seq), args.g))
        elif args.f == 'p':
            with open(args.IN) as f:
                for alns in RelaxedPhylipParser(f):
                    df.at[alns[0][0], alns[1][0]] = df.at[alns[1][0], alns[0][0]] = '{:.8f}'.format(seq_identity(alns[0][1], alns[1][1], args.g))

        for i in df.columns:
            df.at[i, i] = '{:.8f}'.format(1.0)
        df_to_phylip_matrix(df, args.o)
    # perform pairwise alignment
    else:
        FA = [x for x in SeqIO.parse(args.IN, 'fasta')]
        pairwise_alignment_through_fasta(FA, args.aligner, args.opts, args.thread, args.o)

    if args.o != sys.stdout:
        args.o.close()

if __name__ == '__main__':
    main()
