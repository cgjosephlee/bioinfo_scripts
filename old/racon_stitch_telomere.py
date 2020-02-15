#!/usr/bin/env python3
'''
run `trf_telomere_finder.py` beforehand to generate required files and pick proper telomeric sequence
require minimap2

Ideology:
    1. parse trf output of before-racon assembly to find sequences with telomeres
    2. run minimap2 to find the overhangs between before- and after-racon assemblies
    3. stitch back the overhangs with telomere to post-racon assembly
'''

import sys
import argparse
import subprocess as sp
from collections import defaultdict
from Bio import SeqIO

# My_script
from trf_telomere_finder import parse_faidx, parse_trf_output

parser = argparse.ArgumentParser(description='stitch telomeres back for post-racon assembly')
parser.add_argument('fa1', type=str, help='raw assembly')
parser.add_argument('fa2', type=str, help='post-racon assembly')
parser.add_argument('seq', type=str, help='telomeric sequence')
parser.add_argument('-t', type=int, default=8, help='threads for minimap2')
args = parser.parse_args()

fa1 = args.fa1  # raw
fa2 = args.fa2  # post-racon
telomeric_seq = args.seq
threads = args.t

# parse trf output to get original scaffolds with telomeres
try:
    with open(fa1 + '.fai') as f:
        fa1_ctg_len = parse_faidx(f)

    trf_opt = [2, 7, 7, 80, 10, 50, 500]  # default
    trf_out = '{}.{}.dat'.format(fa1, '.'.join([str(x) for x in trf_opt]))
    with open(trf_out) as f:
        FoundInSTART, FoundInEND, FoundInMID = parse_trf_output(f, fa1_ctg_len,
                                                                repeat_seq=telomeric_seq,
                                                                silent=True)
except FileNotFoundError:
    print('Run `trf_telomere_finder.py` first.', file=sys.stderr)
    raise

FoundInBoth = FoundInSTART.intersection(FoundInEND)
FoundInEither = FoundInSTART.union(FoundInEND)

# run minimap2
# ideally racon do not change contig ID
minimap_out = {}
cmd = ['minimap2', '-x', 'asm5', '-t', str(threads), fa1, fa2]
proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.DEVNULL)
for line in proc.stdout:
    line = line.decode().strip().split()
    que = line[0]
    ref = line[5]
    if ref in FoundInEither and ref == que:
        # que_pos = (int(line[2]), int(line[3]))
        # ref_pos = (int(line[7]), int(line[8]))
        minimap_out[ref] = (int(line[2]), int(line[3]), int(line[7]), int(line[8]))
if proc.returncode:
    raise sp.CalledProcessError(proc.returncode, proc.args)

# parse raw fasta and store telomeric sequences
telomere_storage = defaultdict(dict)
for rec in SeqIO.parse(fa1, 'fasta'):
    if rec.id in FoundInSTART:
        telomere_storage[rec.id]['START'] = rec.seq[:minimap_out[rec.id][2]]
    if rec.id in FoundInEND:
        telomere_storage[rec.id]['END'] = rec.seq[minimap_out[rec.id][3]:]

# parse post-racon fasta and stitch telomere
for rec in SeqIO.parse(fa2, 'fasta'):
    if rec.id in FoundInEither:
        map_coord = minimap_out[rec.id]
        if rec.id in FoundInBoth:
            print('>{}\n{}'.format(
                rec.id,
                telomere_storage[rec.id]['START'] + rec.seq[map_coord[0]:map_coord[1]] + telomere_storage[rec.id]['END']
            ))
        elif rec.id in FoundInSTART:
            print('>{}\n{}'.format(
                rec.id,
                telomere_storage[rec.id]['START'] + rec.seq[map_coord[0]:]
            ))
        elif rec.id in FoundInEND:
            print('>{}\n{}'.format(
                rec.id,
                rec.seq[:map_coord[1]] + telomere_storage[rec.id]['END']
            ))
        else:
            raise OSError('something wrong!')
        print('Seq: {} stitched.'.format(rec.id), file=sys.stdout)
    else:
        print('>{}\n{}'.format(rec.id, rec.seq))
