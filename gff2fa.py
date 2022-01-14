#!/usr/bin/env python3
'''
Generate cds and protein sequences from gff and fasta file.
This script assumes a single mRNA for a gene.
'''

import sys
import argparse
import re
from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError

def get_strand(s):
    if s == '-':
        return 0
    else:
        return 1

def get_flags(s):
    q = {}
    s = s.split(';')
    for f in s:
        if f:
            f = f.split('=')
            q[f[0]] = f[1]
    return q

def read_gff(handle):
    # yield lines of a full gene record
    # assuming gene informations are followed by a gene feature line, not looking parent-child structure
    saved = []
    for line in handle:
        if line.startswith('#'):
            continue
        line = line.strip().split()
        if line[2] in ('region'):  # line to ignore, maybe more
            continue
        if line[2] == 'gene':
            if len(saved) > 0:
                yield saved
            saved = []
            saved.append(line)
        else:
            saved.append(line)
    if len(saved) > 0:
        yield saved

def parse_gff_gene_lines(lines):
    assert lines[0][2] == 'gene'  # start with gene feature
    gene = {
        'chr':    lines[0][0],
        'strand': get_strand(lines[0][6]),
        # 'start':  int(lines[0][3]),
        # 'end':    int(lines[0][4]),
        'INFO':   get_flags(lines[0][8]),  # or called flags, attributes
        'exon':   [],
        'CDS':    []
    }
    if lines[1][2] in ('mRNA'):
        gene['type'] = lines[1][2]
        gene['is_coding'] = True
    elif lines[1][2] in ('tRNA', 'rRNA'):
        gene['type'] = lines[1][2]
        gene['is_coding'] = False
    else:
        for line in lines[2:]:
            if line[2] in ('mRNA', 'tRNA', 'rRNA'):
                # mRNA feature not followed by gene?
                raise ValueError('Gff is not sorted?')
        # unknown gene type, just try to translate
        print('WARNING: No gene type (e.g. mRNA) line. ({})'.format(' '.join(lines[0])), file=sys.stderr)
        gene['type'] = None
        gene['is_coding'] = True
    # what we actually need for a gene
    num_gene_type = 0
    for line in lines[1:]:
        if line[2] == 'exon':
            gene['exon'].append((int(line[3]) - 1, int(line[4])))
        elif line[2] == 'CDS':
            gene['CDS'].append((int(line[3]) - 1, int(line[4])))
        elif line[2] in ('mRNA', 'tRNA', 'rRNA'):
            num_gene_type += 1
    if len(gene['CDS']) == 0 and len(gene['exon']) == 0:
        # raise ValueError('No exon nor CDS?')
        print('ERROR: No exon nor CDS? ({})'.format(' '.join(lines[0])), file=sys.stderr)
    if num_gene_type > 1:
        raise NotImplementedError('Not support alternative splicing gene.')
    return gene

def read_fa(file):
    with open(file) as FA:
        FA_dict = SeqIO.to_dict(SeqIO.parse(FA, 'fasta'))
    return FA_dict

def parse_arg():
    parser = argparse.ArgumentParser(description='generate cds and protein sequences from gff and fasta')
    parser.add_argument('-f', metavar='fasta', type=str, required=True, help='fasta file')
    parser.add_argument('-g', metavar='gff', type=str, required=True, help='gff3 file')
    parser.add_argument('-o', metavar='output prefix', default='',
                        help='prefix for output filename (basename of fasta)')
    parser.add_argument('-p', metavar='title prefix', default='',
                        help='prefix for sequence titles (None)')
    parser.add_argument('--field', metavar='field', default='Name',
                        help='field for naming sequences (ID, Name, etc.) (Name)')
    parser.add_argument('--table', metavar='code', type=int, default=4,
                        help='genetic code table for translation (4, for mito)')
    # https://biopython.org/docs/1.75/api/Bio.Seq.html#Bio.Seq.translate
    parser.add_argument('--no-tostop', action='store_false',
                        help='[debug] to_stop=False')
    parser.add_argument('--no-cds', action='store_false',
                        help='[debug] cds=False')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_arg()
    in_seq_file = args.f
    in_gff_file = args.g
    field = args.field
    if args.o == '':
        prefix = re.sub(r'.fa$|.fasta$', '', in_seq_file)
    else:
        prefix = args.o
    if args.p == '':
        title_prefix = args.p
    else:
        title_prefix = args.p + '_'

    # input
    GFF_handle = open(in_gff_file)
    FA_dict = read_fa(in_seq_file)

    # output
    nt_out = open(prefix + '.cds.fa', 'w')
    prot_out = open(prefix + '.prot.fa', 'w')

    num_genes = 0
    num_proteins = 0
    for glines in read_gff(GFF_handle):
        gene = parse_gff_gene_lines(glines)
        try:
            name = gene['INFO'][field]
        except KeyError:
            # todo
            raise
        header = '{}{}'.format(title_prefix, name)
        nt_seq = ''
        prot_seq = ''
        # try CDS first
        if len(gene['CDS']) > 0:
            if not gene['is_coding']:
                print('WARNING: Not sure if this is a coding gene. ({})'.format(name), file=sys.stderr)
            for s, e in gene['CDS']:
                nt_seq += FA_dict[gene['chr']].seq[s:e]
        elif len(gene['exon']) > 0:
            for s, e in gene['exon']:
                nt_seq += FA_dict[gene['chr']].seq[s:e]
        else:
            # raise ValueError('No exon nor CDS?')
            continue
        if gene['strand'] == 0:
            nt_seq = nt_seq.reverse_complement()
        print('>{}\n{}'.format(header, nt_seq), file=nt_out)
        if gene['is_coding']:
            try:
                prot_seq = nt_seq.translate(table=args.table, to_stop=args.no_tostop, cds=args.no_cds)
                print('>{}\n{}'.format(header, prot_seq), file=prot_out)
            except TranslationError as e:
                print('WARNING: {} ({})'.format(e, name), file=sys.stderr)

    GFF_handle.close()
    nt_out.close()
    prot_out.close()
