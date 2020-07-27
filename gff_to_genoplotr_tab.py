#!/usr/bin/env python3
'''
convert gff to tsv for genoplotr::read_dna_seg_from_tab(tab) or genoplotr::dna_seq(data.frame(tab)).
'''

import sys
import argparse

def get_qualifier(s):
    q = {}
    s = s.split(';')
    for f in s:
        f = f.split('=')
        q[f[0]] = f[1]
    return q

parser = argparse.ArgumentParser(description='convert gff file for genoplotR')
parser.add_argument('gff', type=str,
                    help='gff file')
parser.add_argument('-f', type=str, default='ID',
                    help='use this flag in gene feature line as name (ID)')
parser.add_argument('-l', type=int, choices=[1,2,3], default=1,
                    help='detail level (1: gene (default), 2: exon, 3: intron)')
parser.add_argument('--CDS', action='store_true',
                    help='use CDS prior to exon if both exist')
args = parser.parse_args()

in_gff = open(args.gff)

gff_info_to_name = args.f
level = args.l

grob_gene      = 'arrows'
grob_head_exon = 'head_exons'  # a hacked grob
grob_exon      = 'exons'
grob_intron    = 'introns'

default_col    = 'black'
default_fill   = 'red'

# keep chr info and may be used for filtering
# col names should start with name, start, end, strand
gpr_cols = [
    'name',
    'start',
    'end',
    'strand',
    'chr',
    'col',
    'fill',
    'gene_type'
]

def output_to_df(exons, level):
    outs = []
    if strand == '+':
        for n in range(len(exons)):
            e = exons[n]
            if len(exons) == 1:  # only one
                outs.append([gene_id, e[0], e[1], 1, chr, default_col, default_fill, grob_head_exon])
            elif n == 0:  # first one
                outs.append([gene_id, e[0], e[1], 1, chr, default_col, default_fill, grob_exon])
            elif n == len(exons) - 1:  # last one
                e0 = exons[n-1]
                if level == 3:
                    outs.append([gene_id, e0[1]+1, e[0]-1, 1, chr, default_col, default_fill, grob_intron])
                outs.append([gene_id, e[0],    e[1],   1, chr, default_col, default_fill, grob_head_exon])
            else:
                e0 = exons[n-1]
                if level == 3:
                    outs.append([gene_id, e0[1]+1, e[0]-1, 1, chr, default_col, default_fill, grob_intron])
                outs.append([gene_id, e[0],    e[1],   1, chr, default_col, default_fill, grob_exon])
    elif strand == '-':
        for n in range(len(exons)):
            e = exons[n]
            if n == 0:
                outs.append([gene_id, e[0], e[1], -1, chr, default_col, default_fill, grob_head_exon])
            else:
                e0 = exons[n-1]
                if level == 3:
                    outs.append([gene_id, e0[1]+1, e[0]-1, -1, chr, default_col, default_fill, grob_intron])
                outs.append([gene_id, e[0],    e[1],   -1, chr, default_col, default_fill, grob_exon])
    for out in outs:
        print('\t'.join([str(x) for x in out]))

print('\t'.join(gpr_cols))
gene_id = ''
if level == 1:
    for line in in_gff:
        if line.startswith('#'):
            continue
        line = line.strip().split()
        if line[2] == 'gene':
            chr    = line[0]
            strand = line[6]
            st     = int(line[3])
            ed     = int(line[4])
            try:
                gene_id = get_qualifier(line[8])[gff_info_to_name]
            except KeyError:
                gene_id = f'NA({st}-{ed})'
                print(f'Flag ({gff_info_to_name}) not found! Use {gene_id} .', file=sys.stderr)
            if strand == '+':
                out = [gene_id, st, ed, 1, chr, default_col, default_fill, grob_gene]
            else:
                out = [gene_id, st, ed, -1, chr, default_col, default_fill, grob_gene]
            print('\t'.join([str(x) for x in out]))
else:
    for line in in_gff:
        if line.startswith('#'):
            continue
        line = line.strip().split()
        if line[2] == 'gene':
            if gene_id:
                if len(CDS) == 0 and len(exons) == 0:
                    raise ValueError('No exon or CDS feature?')
                if (len(CDS) > 0 and args.CDS) or len(exons) == 0:
                    output_to_df(CDS, level)
                else:
                    output_to_df(exons, level)

            chr     = line[0]
            strand  = line[6]
            st      = int(line[3])
            ed      = int(line[4])
            exons   = []
            CDS     = []
            try:
                gene_id = get_qualifier(line[8])[gff_info_to_name]
            except KeyError:
                gene_id = f'NA({st}-{ed})'
                print(f'Flag ({gff_info_to_name}) not found! Use {gene_id} .', file=sys.stderr)
        elif line[2] == 'exon':
            exons.append([int(line[3]), int(line[4])])
        elif line[2] == 'CDS':
            CDS.append([int(line[3]), int(line[4])])
    # last one
    if len(CDS) == 0 and len(exons) == 0:
        raise ValueError('No exon or CDS feature?')
    if (len(CDS) > 0 and args.CDS) or len(exons) == 0:
        output_to_df(CDS, level)
    else:
        output_to_df(exons, level)
in_gff.close()
