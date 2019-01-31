#!/usr/bin/env python3

'''
Calculate pairwise identity of all sequence combinations in the fasta.

NOTE:
- For protein sequences, straightforward identity was calculated.
- Muscle v3.8.31 has a bug in aligning long sequences. Upgrade to v3.8.341 or higher to fix this.
TODO:
- Support protein sequences.
'''

import sys
import os
import io
import argparse
import re
from Bio import AlignIO, SeqIO
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
    parser = argparse.ArgumentParser(description='generate pairwise sequence identity from fasta')
    parser.add_argument('IN', metavar='fasta', help='fasta file')
    parser.add_argument('-m', metavar='mode', dest='mode', type=int, choices=[1, 2, 3], default=1,
                        help='run mode [1: do alignment, 2: input is aligned, 3: parse previous output only] (default: 1)')
    parser.add_argument('-t', metavar='thread', dest='thread', type=int, default=1,
                        help='threads for alignment (default: 1)')
    parser.add_argument('-f', metavar='format', dest='outfmt', type=int, choices=[1, 2, 3], default=1,
                        help='output format [1: long table, 2: wide table for global identity, 3: wide table for local (gapless) identity] (default: 1)')
    parser.add_argument('-o', metavar='file', dest='OUT', type=argparse.FileType('w'), default=sys.stdout,
                        help='output file (default: stdout)')
    parser.add_argument('--aligner', metavar='prog', choices=['muscle', 'mafft'], default='muscle',
                        help='alignment program [muscle, mafft] (default: muscle)')
    # argparse do not allow value starts with '-', a workaround is to assign equal sign `--opts='--var'`, or add a leadng space `--opts ' --var'`
    parser.add_argument('--opts', metavar='options', action=AddOptAction, default=[],
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

def pairwise_identity(s1, s2, type='nuc'):
    '''
    s1 and s2 are strings! Calculate identity only.
    Ignore positions consisting only gaps in both sequences.
    N to N would be a mismatch.
    global: consider gap in either one of sequences as mismatch
    local: ignore all gaps (gapless)
    '''
    try:
        if len(s1) != len(s2):
            raise ValueError('Unequal length in sequences!')
        aln_len = glb_len = loc_len = len(s1)
        s1_len = len(s1.replace('-', ''))
        s2_len = len(s2.replace('-', ''))

        s1 = s1.upper()
        s2 = s2.upper()
        match = 0
        if type == 'nuc':
            for pos in range(aln_len):
                i, j = s1[pos], s2[pos]
                if i == j == '-':
                    glb_len -= 1
                    loc_len -= 1
                elif i == '-' or j == '-':
                    loc_len -= 1
                elif i == j and i != 'N':
                    match += 1
        glb_id = round(match / glb_len, 4)
        loc_id = round(match / loc_len, 4)
    except ZeroDivisionError as e:
        glb_id = loc_id = float(0)
    return [s1_len, s2_len, glb_len, glb_id, loc_len, loc_id]

def pairwise_alignment(s1, s2, prog='muscle', opts=[]):
    '''
    s1 and s2 are SeqRrecord objects. Perform alignment and calculate identity.
    '''
    records = [s1, s2]
    try:
        if prog == 'muscle':
            cmd = ['muscle', '-quiet'] + opts
            proc = sp.Popen(cmd, stdin=sp.PIPE, stdout=sp.PIPE, universal_newlines=True)
            SeqIO.write(records, proc.stdin, 'fasta')
            proc_out, proc_err = proc.communicate()
            proc_out = io.StringIO(proc_out)
            if proc.returncode == 0:
                # assuming the order of muscle outputs would not change while only 2 sequences
                aln_fa = AlignIO.read(proc_out, 'fasta')
                out = [aln_fa[0].id, aln_fa[1].id] + pairwise_identity(str(aln_fa[0].seq), str(aln_fa[1].seq))
            else:
                raise ValueError
        elif prog == 'mafft':
            cmd = ['mafft', '--quiet'] + opts
            # mafft don't support stdin, so write a temporary fasta
            tmp_fa = tempfile.mkstemp(suffix='.fa', text=True)[1]
            SeqIO.write(records, tmp_fa, 'fasta')
            cmd.append(tmp_fa)
            proc = sp.Popen(cmd, stdout=sp.PIPE, universal_newlines=True)
            proc_out, proc_err = proc.communicate()
            proc_out = io.StringIO(proc_out)
            os.remove(tmp_fa)
            if proc.returncode == 0:
                aln_fa = AlignIO.read(proc_out, 'fasta')
                out = [aln_fa[0].id, aln_fa[1].id] + pairwise_identity(str(aln_fa[0].seq), str(aln_fa[1].seq))
            else:
                raise ValueError
        return out
    except ValueError as e:
        print(e, file=sys.stderr)
        print('Sequence 1: {}\n'
              'Sequence 2: {}'.format(s1.id, s2.id), file=sys.stderr)
        print(proc_err, file=sys.stderr)
        sys.exit(1)
    except UnicodeDecodeError as e:
        print(e, file=sys.stderr)
        print('Sequence 1: {}\n'
              'Sequence 2: {}'.format(s1.id, s2.id), file=sys.stderr)
        sys.exit(1)

def main():
    args = handle_args()

    ### main process ###
    # do alignment
    if args.mode == 1:
        check_program(args.aligner)

        FA = [x for x in SeqIO.parse(args.IN, 'fasta')]
        seq_num = len(FA)

        # multi-threading
        def worker():
            while not queue.empty():
                s1, s2, n = queue.get()
                outs[n] = pairwise_alignment(s1, s2, args.aligner, args.opts)

        job_n = 0
        queue = Queue()
        for n1 in range(seq_num):
            for n2 in range(n1, seq_num):
                queue.put([FA[n1], FA[n2], job_n])
                job_n += 1

        total_runs = queue.qsize()
        outs = [None] * queue.qsize()
        threads = [Thread(target=worker) for x in range(args.thread)]
        print('Open {} threads, running {} jobs...'.format(args.thread, total_runs), file=sys.stderr)
        for th in threads:
            th.start()
        for th in threads:
            # th.join()
            while th.is_alive():
                queueing = queue.qsize()
                print('{} jobs are queueing...'.format(queueing), end='\r', file=sys.stderr)
                time.sleep(0.5)
        print('', file=sys.stderr)
        print('Finish!', file=sys.stderr)
    # no alignment
    elif args.mode == 2:
        FA = AlignIO.read(args.IN, 'fasta')
        seq_num = len(FA)
        total_runs = queueing = int(seq_num * (seq_num + 1) / 2)
        print('Found {} sequences, {} jobs in total.'.format(seq_num, total_runs), file=sys.stderr)
        outs = []
        for n1 in range(seq_num):
            for n2 in range(n1, seq_num):
                outs.append([FA[n1].id, FA[n2].id] + pairwise_identity(str(FA[n1].seq), str(FA[n2].seq)))
                queueing -= 1
                print('{} jobs are queueing...'.format(queueing), end='\r', file=sys.stderr)
        print('', file=sys.stderr)
        print('Finish!', file=sys.stderr)
    # parse only
    elif args.mode == 3:
        if args.outfmt == 1:
            print('In parse only mode, please specify other output format.', file=sys.stderr)
            sys.exit(1)
        with open(args.IN) as f:
            outs = []
            seq_list = []
            for line in f.readlines():
                if line.startswith('#'):
                    continue
                else:
                    line = line.strip().split()
                    line[5] = float(line[5])
                    line[7] = float(line[7])
                    if line[0] not in seq_list:
                        seq_list.append(line[0])
                        line[2] = int(line[2])
                        outs.append([line[0], line[0], line[2], line[2], line[2], 1.0, line[2], 1.0])
                    outs.append(line)
            if line[1] not in seq_list:
                seq_list.append(line[1])
                line[3] = int(line[3])
                outs.append([line[1], line[1], line[3], line[3], line[3], 1.0, line[3], 1.0])
        seq_num = len(seq_list)

    glb_val = []
    loc_val = []
    for o in outs:
        if o[0] != o[1]:
            glb_val.append(o[5])
            loc_val.append(o[7])
    glb_mean = sum(glb_val) / len(glb_val)
    loc_mean = sum(loc_val) / len(loc_val)

    ### output results ###
    if args.outfmt == 1:
        print('# number of sequences: {}\n'
              '# mean of global identity: {:.4f}\n'
              '# mean of local identity: {:.4f}\n'
              '# seq1\tseq2\tlen1\tlen2\tglb_len\tglb_id\tloc_len\tloc_id'.format(seq_num, glb_mean, loc_mean), file=args.OUT)
        for out in outs:
            if out[0] != out[1]:
                print('\t'.join(str(x) for x in out), file=args.OUT)
    elif args.outfmt == 2 or args.outfmt == 3:
        df = pd.DataFrame()
        if args.outfmt == 2:
            for out in outs:
                df.at[out[1], out[0]] = out[5]
            s = 'global'
            i = glb_mean
        elif args.outfmt == 3:
            for out in outs:
                df.at[out[1], out[0]] = out[7]
            s = 'local'
            i = loc_mean
        print('# number of sequences: {}\n'
              '# mean of {} identity: {:.4f}'.format(seq_num, s, i), file=args.OUT)
        df.to_csv(args.OUT, sep='\t')

    # avoide to terminate ipython console
    if args.OUT != sys.stdout:
        args.OUT.close()

if __name__ == '__main__':
    main()
