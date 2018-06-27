#!/usr/bin/env python3

import sys
import argparse
import re
from Bio import AlignIO

def parse_args():
    parser = argparse.ArgumentParser(description='Extract introns from aligned fasta.\n'
                                    'Take given reference sequence (CDS) to find exon boundaries, and output ungapped intron sequences for each record.',
                                    formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('fasta', help='aligned fasta file')
    parser.add_argument('-l', action='store_true', help='output length information')
    parser.add_argument('-n', type=int, default=0, help='#seq as reference (0)')
    return parser.parse_args()

def get_boundaries():
    exon = re.compile(r'^(\w+)')
    intron = re.compile(r'^(-+)')
    boundaries = [0]
    is_exon = []
    k = 0
    while k < FA_aln.get_alignment_length():
        if exon.match(refREC[k:]):
            s = exon.match(refREC[k:]).group(0)
            k += len(s)
            boundaries.append(k)
            is_exon.append(True)
        else:
            s = intron.match(refREC[k:]).group(0)
            k += len(s)
            boundaries.append(k)
            is_exon.append(False)
    return boundaries, is_exon

# ordered by original sequence order
def ordered_by_input():
    for rec in FA_aln:
        if rec.id != refID:
            for j in range(len(is_exon)):
                if is_exon[j]:
                    pass
                elif not is_exon[j]:
                    n = (j+1)/2
                    print('>{}_intron{}\n{}'.format(rec.id, int(n), rec.seq[boundaries[j]:boundaries[j+1]].ungap('-')))

# ordered by intron number
def ordered_by_feature():
    for j in range(len(is_exon)):
        if is_exon[j]:
            pass
        elif not is_exon[j]:
            local_aln = FA_aln[:,boundaries[j]:boundaries[j+1]]
            for rec in local_aln:
                if rec.id != refID:
                    n = (j+1)/2
                    print('>{}_intron{}\n{}'.format(rec.id, int(n), rec.seq.ungap('-')))

def print_length():
    for j in range(len(is_exon)):
        if is_exon[j]:
            n = (j+2)/2
            print("exon  {:<2} len:{:<5} {{{}..{}}}".format(int(n), boundaries[j+1]-boundaries[j], boundaries[j], boundaries[j+1]-1))
        elif not is_exon[j]:
            n = (j+1)/2
            print("intron{:<2} len:{:<5} {{{}..{}}}".format(int(n), boundaries[j+1]-boundaries[j], boundaries[j], boundaries[j+1]-1))

args = parse_args()
inFA = args.fasta
refNUM = args.n # -1 for last

FA_aln = AlignIO.read(inFA, 'fasta')

refID = FA_aln[refNUM].id
refREC = str(FA_aln[refNUM].seq)

boundaries, is_exon = get_boundaries()
if args.l:
    print_length()
else:
    #ordered_by_input()
    ordered_by_feature()
