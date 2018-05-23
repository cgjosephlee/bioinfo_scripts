#!/usr/bin/env python3

import sys
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='split fasta by usearch cluster file')
parser.add_argument('fasta', type=str, help='input fasta')
parser.add_argument('uc_cluster', type=str, help='usearch cluster file')
parser.add_argument('-p', metavar='prefix', type=str, default='out', help='output prefix (out)')
parser.add_argument('-n', metavar='NUM', type=int, default=1, help='filter clusters <= n sequences (1)')
args = parser.parse_args()

in_fa = args.fasta
in_uc = args.uc_cluster
prefix = args.p
seq_num = args.n

with open(in_uc) as f:
    clusters = {}
    seq_ids = {}
    for line in f:
        if not line.startswith('C'):
            line = line.strip()
            line = line.split('\t')
            cluster = int(line[1])
            id = line[8]
            clusters.setdefault(cluster, []).append(id)
            seq_ids[id] = cluster

outfiles = {}
for key in clusters.keys():
    if len(clusters[key]) > seq_num:
        outfiles[key] = open('{}_cluster{:0>2d}.fa'.format(prefix, int(key)), 'w')
    else:
        print('cluster {} omitted (size: {})'.format(int(key), len(clusters[key])))

for rec in SeqIO.parse(in_fa, 'fasta'):
    if len(clusters[seq_ids[rec.id]]) > seq_num:
        #print('>{}\n{}'.format(rec.id, rec.seq), file = outfiles[seq_ids[rec.id]])
        SeqIO.write(rec, outfiles[seq_ids[rec.id]], 'fasta')

for f in outfiles.values():
    f.close()

