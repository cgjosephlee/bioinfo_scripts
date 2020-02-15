#!/usr/bin/env python3
'''
convert gff to tsv for genoplotr::read_dna_seg_from_tab(tab) or genoplotr::dna_seq(data.frame(tab)).
'''

import sys

def get_qualifier(s):
    q = {}
    s = s.split(';')
    for f in s:
        f = f.split('=')
        q[f[0]] = f[1]
    return q

in_gff = sys.argv[1]

gff_exon_feature = 'CDS'
gff_info_to_name = 'ID'  # in gene feature line
# gff_info_to_name = 'locus_tag'

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

def output_to_df():
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
                outs.append([gene_id, e0[1]+1, e[0]-1, 1, chr, default_col, default_fill, grob_intron])
                outs.append([gene_id, e[0],    e[1],   1, chr, default_col, default_fill, grob_head_exon])
            else:
                e0 = exons[n-1]
                outs.append([gene_id, e0[1]+1, e[0]-1, 1, chr, default_col, default_fill, grob_intron])
                outs.append([gene_id, e[0],    e[1],   1, chr, default_col, default_fill, grob_exon])
    elif strand == '-':
        for n in range(len(exons)):
            e = exons[n]
            if n == 0:
                outs.append([gene_id, e[0], e[1], -1, chr, default_col, default_fill, grob_head_exon])
            else:
                e0 = exons[n-1]
                outs.append([gene_id, e0[1]+1, e[0]-1, -1, chr, default_col, default_fill, grob_intron])
                outs.append([gene_id, e[0],    e[1],   -1, chr, default_col, default_fill, grob_exon])
    for out in outs:
        print('\t'.join([str(x) for x in out]))

print('\t'.join(gpr_cols))
gene_id = ''
with open(in_gff) as GFF:
    for line in GFF:
        if line.startswith('#'):
            pass
        line = line.strip().split()
        if line[2] == 'gene':
            if gene_id:
                output_to_df()
            chr     = line[0]
            gene_id = get_qualifier(line[8])[gff_info_to_name]
            strand  = line[6]
            exons   = []
        elif line[2] == gff_exon_feature:
            exons.append([int(line[3]), int(line[4])])
output_to_df()
