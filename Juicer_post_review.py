#!/usr/bin/env python3

'''
After generating `*.review.assembly` in Juicerbox, generate final fasta and .agp file according to the `*.review.assembly`.
The output fasta is sorted by length.
Overhang Ns are not trimmed.
'''

import sys
from Bio import SeqIO
import textwrap

if len(sys.argv) != 4:
    print('Usage: {} <assembly> <orig_fasta> <out_prefix>', sys.stderr)
    sys.exit(1)

def argsort(seq):
    # https://stackoverflow.com/questions/3382352/equivalent-of-numpy-argsort-in-basic-python
    return sorted(range(len(seq)), key=seq.__getitem__)

in_ass = sys.argv[1]
in_fa = sys.argv[2]
out_file_prefix = sys.argv[3]
out_fasta = out_file_prefix + '.fasta'
out_agp = out_file_prefix + '.agp'
out_fa_prefix = 'HiC_scaffold'
gap_len = 500
gap = 'N' * gap_len

# read *.assembly
with open(in_ass) as ass:
    ass_info = {}
    ass_superscaf = []
    len_superscaf = []
    ref_scaf = ''
    ref_st = 0
    ref_ed = 0
    for line in ass:
        line = line.strip()
        if line.startswith('>'):
            line = line[1:].split(' ')

            last_ref_scaf = ref_scaf

            ref_scaf = line[0].split(':::fragment')[0]
            ref_len = int(line[2])
            if ref_scaf == last_ref_scaf:
                ref_st = ref_ed
            else:
                ref_st = 0
            ref_ed = ref_st + ref_len
            ass_info[line[1]] = [ref_scaf, [ref_st, ref_ed], ref_len]
        else:
            ass_superscaf.append(line.split(' '))
            len_superscaf.append(sum([ass_info[x.lstrip('-')][2] for x in ass_superscaf[-1]]))

# filter short sequences?
sort_superscaf = argsort(len_superscaf)[::-1]

# write .agp, 1-base format
# <obj> <obj_beg> <obj_end> <part_num> <type> <comp_id> <comp_beg> <comp_end> <orientation>
#                                      N      <gap_len> scaffold   yes        na
fa_counter = 1
with open(out_agp, 'w') as f:
    print('##agp-version 2.0', file=f)
    for s in sort_superscaf:
        obj_scaf = '{}_{}'.format(out_fa_prefix, fa_counter)
        obj_st = 0
        obj_ed = 0
        part = 1
        for i in ass_superscaf[s]:
            # insert gap record
            if part > 1:
                obj_ed = obj_st + gap_len
                out = [obj_scaf, obj_st + 1, obj_ed, part, 'N', gap_len, 'scaffold', 'yes', 'Hi_C']
                print('\t'.join(str(x) for x in out), file=f)
                part += 1
                obj_st = obj_ed

            if i.startswith('-'):
                strand = '-'
                i = i[1:]
            else:
                strand = '+'
            comp_rec = ass_info[i]
            comp_scaf = comp_rec[0]
            comp_st = comp_rec[1][0]
            comp_ed = comp_rec[1][1]
            comp_len = comp_rec[2]
            obj_ed = obj_st + comp_len

            out = [obj_scaf, obj_st + 1, obj_ed, part, 'W', comp_scaf, comp_st + 1, comp_ed, strand]
            print('\t'.join(str(x) for x in out), file=f)

            part += 1
            obj_st = obj_ed
        fa_counter += 1

# write fasta
FA = SeqIO.to_dict(SeqIO.parse(in_fa, 'fasta'))
fa_counter = 1
out_seq = []
with open(out_fasta, 'w') as f:
    for s in sort_superscaf:
        for i in ass_superscaf[s]:
            if i.startswith('-'):
                info = ass_info[i[1:]]
                out_seq.append(str(FA[info[0]].seq[info[1][0]:info[1][1]].reverse_complement()))
            else:
                info = ass_info[i]
                out_seq.append(str(FA[info[0]].seq[info[1][0]:info[1][1]]))
        print('>{}_{}'.format(out_fa_prefix, fa_counter), file=f)
        # print(gap.join(out_seq, file=f)
        print(textwrap.fill(gap.join(out_seq), width=80), file=f)
        fa_counter += 1
        out_seq = []
