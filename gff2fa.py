#!/usr/bin/env python3

import sys
import argparse
import re
from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError
# from pprint import pprint

# This script assumes a single progeny mRNA feature for a gene
#
# --CDS mode
# gene features have no sibling mRNA/exon, but only CDS (e.g. AF402141.1)


def parse_arg():
    parser = argparse.ArgumentParser(description='generate cds and protein sequences from gff and fasta')
    parser.add_argument('-f', metavar='fasta', type=str, required=True, help='fasta file')
    parser.add_argument('-g', metavar='gff', type=str, required=True, help='gff3 file')
    parser.add_argument('-o', metavar='output prefix', default='',
                        help='prefix for output filename (basename of fasta)')
    parser.add_argument('-p', metavar='title prefix', default='',
                        help='prefix for sequence titles (None)')
    parser.add_argument('--field', metavar='field', default='Name',
                        help='feature field for naming sequences (ID, Name, product, etc.) (Name)')
    parser.add_argument('--table', metavar='code', type=int, default=4,
                        help='genetic code table for translation (4, for mito)')
    parser.add_argument('--CDS', action='store_true',
                        help='[debug] for some exceptions (see descriptions in script)')
    parser.add_argument('--no-filter', action='store_true',
                        help='[debug] disable filter')
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


def parse_gff(file, CDS_mode):
    n = -1
    gene_lt = []
    with open(file) as GFF:
        for line in GFF.readlines():
            line = line.strip()
            if not line.startswith('#') and not line == '':
                line = line.split('\t')
                qualifiers = get_qualifier(line[8])
                if line[2] == 'gene':
                    n += 1
                    gene_lt.append({})
                    gene_lt[n]['chr'] = line[0]
                    gene_lt[n]['strand'] = get_strand(line[6])
                    gene_lt[n]['ID'] = qualifiers.get('ID', '')
                    gene_lt[n]['Name'] = qualifiers.get('Name', '')
                elif line[2] in ['mRNA', 'tRNA', 'rRNA']:
                    gene_lt[n]['type'] = line[2]
                    gene_lt[n]['product'] = qualifiers.get('product', '')
                elif line[2] == 'exon':
                    exon_co = [int(line[3]) - 1, int(line[4])]
                    gene_lt[n].setdefault('location', []).append(exon_co)
                # for some exceptions that genes have no mRNA and exon featurs, but only CDS
                elif line[2] == 'CDS' and CDS_mode:
                    gene_lt[n]['type'] = 'mRNA'
                    gene_lt[n]['product'] = qualifiers.get('product', '')
                    exon_co = [int(line[3]) - 1, int(line[4])]
                    gene_lt[n].setdefault('location', []).append(exon_co)
    return gene_lt


def parse_fa(file):
    with open(file) as FA:
        FA_dict = SeqIO.to_dict(SeqIO.parse(FA, 'fasta'))
    return FA_dict


# handle arguments
args = parse_arg()
in_seq_file = args.f
in_gff_file = args.g
table = args.table
field = args.field
if args.o == '':
    prefix = '.'.join(in_seq_file.split('.')[0:-1])
else:
    prefix = args.o
if args.p == '':
    title_prefix = args.p
else:
    title_prefix = args.p + '_'

# input files
FA_dict = parse_fa(in_seq_file)
genes = parse_gff(in_gff_file, args.CDS)

# filter (temporary)
if not args.no_filter:
    genes = filter(lambda x: not re.match(r'trn[A-Z]|tRNA|orf|ORF', x[field]), genes)

# main
with open(prefix + '.cds.fa', 'w') as nt_out,\
        open(prefix + '.protein.fa', 'w') as prot_out:
    for rec in genes:
        nt_seq = ''
        prot_seq = ''
        try:
            for start, end in rec['location']:
                nt_seq += FA_dict[rec['chr']].seq[start: end]
        except KeyError:
            print('WARNING: gene has no exon ({})'.format(rec[field]), file=sys.stderr)
            continue
        if rec['strand'] == 0:
            nt_seq = nt_seq.reverse_complement()
        print('>{}{}\n{}'.format(title_prefix, rec[field], nt_seq), file=nt_out)
        if rec['type'] == 'mRNA':
            try:
                prot_seq = nt_seq.translate(table=table, to_stop=True, cds=True)
                print('>{}{}\n{}'.format(title_prefix, rec[field], prot_seq), file=prot_out)
            except TranslationError as e:
                print('WARNING: {} ({})'.format(e, rec[field]), file=sys.stderr)
