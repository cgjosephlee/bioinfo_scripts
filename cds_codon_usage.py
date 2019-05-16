#!/usr/bin/env python3
'''
Calculate codon usage for cds sequences.
'''

import sys
import argparse
from Bio import SeqIO
from Bio.Data import CodonTable

aa_3to1 = {'Ala': 'A',
           'Arg': 'R',
           'Asn': 'N',
           'Asp': 'D',
           'Cys': 'C',
           'Gln': 'Q',
           'Glu': 'E',
           'Gly': 'G',
           'His': 'H',
           'Ile': 'I',
           'Leu': 'L',
           'Lys': 'K',
           'Met': 'M',
           'Phe': 'F',
           'Pro': 'P',
           'Ser': 'S',
           'Thr': 'T',
           'Trp': 'W',
           'Tyr': 'Y',
           'Val': 'V',
           '*':   '*'}

parser = argparse.ArgumentParser(description='Calculate codon usage.')
parser.add_argument('fasta',
                    help='fasta file')
parser.add_argument('CodonTable', nargs='?', type=int,
                    help='codon table code according to NCBI (4)')
args = parser.parse_args()

in_fasta = args.fasta
if args.CodonTable:
    table_id = args.CodonTable
else:
    table_id = 4

# a dict, stop codon not included
table = CodonTable.unambiguous_dna_by_id[table_id].forward_table
start_codons = CodonTable.unambiguous_dna_by_id[table_id].start_codons
stop_codons = CodonTable.unambiguous_dna_by_id[table_id].stop_codons
for c in stop_codons:
    table[c] = '*'

aaUsage = {}
for k in table.keys():
    aaUsage.setdefault(table[k], {}).update({k: 0})

aaFreq = {}
for rec in SeqIO.parse(in_fasta, 'fasta'):
    seq = str(rec.seq).upper()
    if len(seq) % 3 != 0 or seq[0:3] not in start_codons:
        print("{} is not a valid CDS, pass.".format(rec.id), file=sys.stderr)
        continue
    else:
        # ignore start codon since always Met
        for k in range(3, len(seq), 3):
            c = seq[k:k+3]
            try:
                aaFreq[table[c]] = aaFreq.setdefault(table[c], 0) + 1
                aaUsage[table[c]][c] += 1
            except KeyError:
                print('Unknow codon {} in {}.'.format(c, rec.id), file=sys.stderr)

for aa3 in aa_3to1.keys():
    aa1 = aa_3to1[aa3]
    dt = sorted(aaUsage[aa1].items(), key=lambda x: x[0])
    for k, v in dt:
        print('{}\t{}\t{}\t{:.2f}\t{}'.format(aa3, aa1, k, v / aaFreq[aa1], v))
