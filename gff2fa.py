#!/usr/bin/env python3

import sys
import argparse
import re
from Bio import SeqIO
# from pprint import pprint


def parse_arg():
    parser = argparse.ArgumentParser(description='generate cds and protein sequences from gff and fasta')
    parser.add_argument('-f', metavar='fasta', type=str, required=True, help='fasta file')
    parser.add_argument('-g', metavar='gff', type=str, required=True, help='gff3 file')
    parser.add_argument('-n', metavar='genetic table', default=4,
                        help='genetic table for translation (4, for mito)')
    parser.add_argument('-p', metavar='prefix', default='out',
                        help='output file name prefix (out)')
    # parser.add_argument('-t', action='store_true',
    #                     help='output translated protein sequences')
    return parser.parse_args()


def get_strand(s):
    if s == '-':
        return 0
    else:
        return 1


def get_qualifier(s):
    q = {}
    s = s.split(';')
    for f in s:
        f = f.split('=')
        q[f[0]] = f[1]
    return q


def parse_gff(file):
    n = -1
    gene_lt = []
    with open(file) as GFF:
        for line in GFF.readlines():
            line = line.strip()
            if not line.startswith('#'):
                line = line.split('\t')
                qualifiers = get_qualifier(line[8])
                if line[2] == 'gene':
                    n += 1
                    gene_lt.append({})
                    gene_lt[n]['chr'] = line[0]
                    gene_lt[n]['ID'] = qualifiers['ID']
                    gene_lt[n]['Name'] = qualifiers['Name']
                    gene_lt[n]['strand'] = get_strand(line[6])
                elif line[2] in ['mRNA', 'tRNA', 'rRNA']:
                    gene_lt[n]['type'] = line[2]
                    gene_lt[n]['product'] = qualifiers.get('product', gene_lt[n]['ID'])
                elif line[2] == 'exon':
                    exon_co = [int(line[3]) - 1, int(line[4])]
                    gene_lt[n].setdefault('location', []).append(exon_co)
    return gene_lt


def parse_fa(file):
    with open(file) as FA:
        FA_dict = SeqIO.to_dict(SeqIO.parse(FA, 'fasta'))
    return FA_dict


args = parse_arg()
in_seq_file = args.f
in_gff_file = args.g
prefix = args.p
table = args.n

FA_dict = parse_fa(in_seq_file)
genes = parse_gff(in_gff_file)
genes = filter(lambda x: not re.match(r'trn|orf', x['ID']), genes)

with open(prefix + '.cds.fa', 'w') as nt_out,\
        open(prefix + '.protein.faa', 'w') as prot_out:
    for rec in genes:
        nt_seq = ''
        prot_seq = ''
        try:
            for start, end in rec['location']:
                nt_seq += FA_dict[rec['chr']].seq[start: end]
        except KeyError:
            print('WARNING: gene has no exon ({})'.format(rec['ID']), file=sys.stderr)
            continue
        if rec['strand'] == 0:
            nt_seq = nt_seq.reverse_complement()
        print('>{}\n{}'.format(rec['ID'], nt_seq), file=nt_out)
        if rec['type'] == 'mRNA':
            if len(nt_seq) % 3 != 0:
                print('WARNING: cds length cannot divided by 3 ({})'.format(rec['ID']), file=sys.stderr)
                continue
            prot_seq = nt_seq.translate(table=table, to_stop=True)
            print('>{}\n{}'.format(rec['ID'], prot_seq), file=prot_out)
