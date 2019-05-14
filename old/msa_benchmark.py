#!/usr/bin/env python3

import sys, os
from Bio import SeqIO, AlignIO
import numpy as np
import time
import cProfile

# def pid1(s1, s2, type='nuc'):
#     '''
#     convert N into gap, and ignore positions consisting only gaps in both sequences
#     global: consider gap in either one of sequences as mismatch
#     local: ignore all gaps (gapless)
#     '''
#     try:
#         if len(s1) != len(s2):
#             raise ValueError('Unequal length in sequences!')
#         aln_len = glb_len = loc_len = len(s1)
#         s1_len = len(s1.seq.ungap('-'))
#         s2_len = len(s2.seq.ungap('-'))

#         s1.seq = s1.seq.upper()
#         s2.seq = s2.seq.upper()
#         match = 0
#         for pos in range(aln_len):
#             # i, j = str(s1[pos]).upper(), str(s2[pos]).upper()
#             i, j = s1[pos], s2[pos]

#             if type == 'nuc':
#                 # i = i.replace('N', '-')
#                 # j = j.replace('N', '-')
#                 if i == j == '-':
#                     glb_len -= 1
#                     loc_len -= 1
#                 elif i == '-' or j == '-':
#                     loc_len -= 1
#                 elif i == j:
#                     match += 1

#         # glb_id = '{:.4f}'.format(match/glb_len)
#         # loc_id = '{:.4f}'.format(match/loc_len)
#         glb_id = round(match / glb_len, 4)
#         loc_id = round(match / loc_len, 4)
#     except ZeroDivisionError as e:
#         # if s1_len == 0:
#         #     print('Sequence "{}" is empty or contains on gaps.'.format(s1.id), file=sys.stderr)
#         # elif s2_len == 0:
#         #     print('Sequence "{}" is empty or contains on gaps.'.format(s2.id), file=sys.stderr)
#         glb_id = loc_id = float(0)
#     return [s1.id, s2.id, s1_len, s2_len, glb_len, glb_id, loc_len, loc_id]

# def pid2(s1, s2, type='nuc'):
#     '''
#     convert N into gap, and ignore positions consisting only gaps in both sequences
#     global: consider gap in either one of sequences as mismatch
#     local: ignore all gaps (gapless)
#     '''
#     try:
#         if len(s1) != len(s2):
#             raise ValueError('Unequal length in sequences!')
#         aln_len = glb_len = loc_len = len(s1)
#         s1_len = len(s1.ungap('-'))
#         s2_len = len(s2.ungap('-'))

#         s1 = s1.upper()
#         s2 = s2.upper()
#         match = 0
#         for pos in range(aln_len):
#             # i, j = str(s1[pos]).upper(), str(s2[pos]).upper()
#             i, j = s1[pos], s2[pos]

#             if type == 'nuc':
#                 # i = i.replace('N', '-')
#                 # j = j.replace('N', '-')
#                 if i == j == '-':
#                     glb_len -= 1
#                     loc_len -= 1
#                 elif i == '-' or j == '-':
#                     loc_len -= 1
#                 elif i == j:
#                     match += 1

#         glb_id = round(match / glb_len, 4)
#         loc_id = round(match / loc_len, 4)
#     except ZeroDivisionError as e:
#         glb_id = loc_id = float(0)
#     return [s1_len, s2_len, glb_len, glb_id, loc_len, loc_id]

def pid3(s1, s2, type='nuc'):
    '''
    convert N into gap, and ignore positions consisting only gaps in both sequences
    global: consider gap in either one of sequences as mismatch
    local: ignore all gaps (gapless)
    '''
    try:
        if len(s1) != len(s2):
            raise ValueError('Unequal length in sequences!')
        aln_len = glb_len = loc_len = len(s1)
        s1_len = len(s1.replace('-', ''))
        s2_len = len(s2.replace('-', ''))

        s1 = s1.upper()
        s2 = s2.upper()
        match = 0
        if type == 'nuc':
            for pos in range(aln_len):
                # i, j = str(s1[pos]).upper(), str(s2[pos]).upper()
                i, j = s1[pos], s2[pos]

                if i == j == '-':
                    glb_len -= 1
                    loc_len -= 1
                elif i == '-' or j == '-':
                    loc_len -= 1
                elif i == j and i != 'N':
                    match += 1

        glb_id = round(match / glb_len, 4)
        loc_id = round(match / loc_len, 4)
    except ZeroDivisionError as e:
        glb_id = loc_id = float(0)
    return [s1_len, s2_len, glb_len, glb_id, loc_len, loc_id]

# def pid4(s1, s2, type='nuc'):
#     '''
#     convert N into gap, and ignore positions consisting only gaps in both sequences
#     global: consider gap in either one of sequences as mismatch
#     local: ignore all gaps (gapless)
#     '''
#     def comp(i, j):
#         g = False
#         l = False
#         m = False
#         if i == j == '-':
#             g = True
#             l = True
#         elif i == '-' or j == '-':
#             l = True
#         elif i == j and i != 'N':
#             m = True
#         return g, l, m
#     try:
#         if len(s1) != len(s2):
#             raise ValueError('Unequal length in sequences!')
#         aln_len = glb_len = loc_len = len(s1)
#         s1_len = len(s1.replace('-', ''))
#         s2_len = len(s2.replace('-', ''))

#         s1 = s1.upper()
#         s2 = s2.upper()
#         match = 0
#         if type == 'nuc':
#             for pos in range(aln_len):
#                 g, l, m = comp(s1[pos], s2[pos])
#                 if g:
#                     glb_len -= 1
#                 if l:
#                     loc_len -= 1
#                 if m:
#                     match += 1

#         glb_id = round(match / glb_len, 4)
#         loc_id = round(match / loc_len, 4)
#     except ZeroDivisionError as e:
#         glb_id = loc_id = float(0)
#     return [s1_len, s2_len, glb_len, glb_id, loc_len, loc_id]

FA = AlignIO.read('test.large.3seq.aln.fa', 'fasta')

# cProfile.run('pid1(FA[1], FA[2])')
# cProfile.run('pid2(FA[1].seq, FA[2].seq)')
cProfile.run('pid3(str(FA[1].seq), str(FA[2].seq))')
# cProfile.run('pid4(str(FA[1].seq), str(FA[2].seq))')

# seq_num = len(FA)
# total_runs = int(seq_num * (seq_num + 1) / 2)

# st = time.clock()
# outs = []
# # n = 0
# for n1 in range(seq_num):
#     for n2 in range(n1, seq_num):
#         # outs.append(pid1(FA[n1], FA[n2]))
#         outs.append(pid2(FA[n1].seq, FA[n2].seq))
#         # n += 1
# ed = time.clock()
# print('{:.5f}s'.format(ed-st))

###########
'''
nul = open(os.devnull, 'w')
s = 'atgc' * 100000
st = time.clock()
for i in range(len(s)):
    print(s[i], file=nul)
ed = time.clock()
print('{:.5f}s'.format(ed-st))

st = time.clock()
s = list(s)
for i in range(len(s)):
    print(s[i], file=nul)
ed = time.clock()
print('{:.5f}s'.format(ed-st))
'''
