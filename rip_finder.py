#!/usr/bin/env python3

import sys
from collections import Counter, defaultdict
from math import ceil
# from fractions import Fraction

from Bio.SeqIO.FastaIO import SimpleFastaParser
import numpy as np

def generate_consensus(seq_2d):
    '''
    input: seq 2d array
    output consensus 1d array
    '''
    # base_to_flag = {
    #     'A': 1 ,
    #     'T': 2 ,
    #     'C': 4 ,
    #     'G': 8 ,
    #     'W': 3 ,
    #     'M': 5 ,
    #     'Y': 6 ,
    #     'R': 9 ,
    #     'K': 10,
    #     'S': 12,
    #     'H': 7 ,
    #     'D': 11,
    #     'V': 13,
    #     'B': 14,
    #     'N': 15
    # }
    flag_to_base = {
        1 : 'A',
        2 : 'T',
        4 : 'C',
        8 : 'G',
        3 : 'W',
        5 : 'M',
        6 : 'Y',
        9 : 'R',
        10: 'K',
        12: 'S',
        7 : 'H',
        11: 'D',
        13: 'V',
        14: 'B',
        15: 'N'
    }
    consensus = np.zeros(seq_2d.shape[1], dtype='U1')
    max_gap = seq_2d.shape[0] - 2
    # iterate columns
    # for i in range(seq_2d.shape[1]):
    #     count = Counter(seq_2d[:,i])
    #     if count['-'] > max_gap:
    #         consensus[i] = '-'
    #         continue
    #     flag_sum = 0
    #     last_n = 0
    #     for char, n in count.most_common():  # list odered by counts
    #         if char == '-':
    #             continue
    #         if n == last_n:
    #             flag_sum += base_to_flag[char]
    #             continue
    #         elif n < last_n:
    #             break
    #         flag_sum += base_to_flag[char]
    #         last_n = n
    #     consensus[i] = flag_to_base[flag_sum]
    # return consensus

    seq_1d_counts = np.zeros((seq_2d.shape[1], 5), dtype=np.uint16)  # 0-65535
    for n, col in enumerate(seq_2d.T):
        count = Counter(col)
        seq_1d_counts[n] = [
            count['A'],
            count['T'],
            count['C'],
            count['G'],
            count['-']
        ]
        if count['-'] > max_gap:
            consensus[n] = '-'
            continue
        c = seq_1d_counts[n, :4]
        consensus[n] = flag_to_base[sum(np.power(2, np.nonzero(c == c.max())[0]))]
    return seq_1d_counts, consensus

def cal_mutation_rate(seq_2d, consensus):
    valid_bases = 0
    match_bases = 0.
    # for a,b in np.nditer([seq_2d, consensus]):
    #     if a != '-' and b != '-':
    #         valid_bases += 1
    #         match_bases += nuc_convert_chance(str(a), str(b), omit_N=False)
    for n, arr in enumerate(seq_2d.T):
        b = consensus[n]
        for a in arr:
            if a != '-' and b != '-':
                valid_bases += 1
                match_bases += nuc_convert_chance(a, b, omit_N=False)
    mu = 1 - (match_bases / valid_bases)
    min_n = ceil(seq_2d.shape[0] * mu)  # minmal mutation number greater than random mutation
    return mu, min_n

###
# calculate trinucletide RIP
###
def nuc_convert_chance(origin, degeneracy, omit_N=True):
    df = {
        ('A', 'A'): 1,
        ('T', 'T'): 1,
        ('C', 'C'): 1,
        ('G', 'G'): 1,
        ('A', 'M'): 1/2,
        ('A', 'R'): 1/2,
        ('A', 'W'): 1/2,
        ('T', 'K'): 1/2,
        ('T', 'W'): 1/2,
        ('T', 'Y'): 1/2,
        ('C', 'M'): 1/2,
        ('C', 'S'): 1/2,
        ('C', 'Y'): 1/2,
        ('G', 'K'): 1/2,
        ('G', 'R'): 1/2,
        ('G', 'S'): 1/2,
        ('A', 'D'): 1/3,
        ('A', 'H'): 1/3,
        ('A', 'V'): 1/3,
        ('T', 'B'): 1/3,
        ('T', 'D'): 1/3,
        ('T', 'H'): 1/3,
        ('C', 'B'): 1/3,
        ('C', 'H'): 1/3,
        ('C', 'V'): 1/3,
        ('G', 'B'): 1/3,
        ('G', 'D'): 1/3,
        ('G', 'V'): 1/3
    }
    if not omit_N:
        df.update({
            ('A', 'N'): 1/4,
            ('T', 'N'): 1/4,
            ('C', 'N'): 1/4,
            ('G', 'N'): 1/4
        })
    try:
        return df[(origin, degeneracy)]
    except KeyError:  # gap, N or others
        return 0

def trinuc_rip_chance(ref_trinuc, qry_trinuc):
    '''
    inputs: np.chararray()
    assert query is unambiguous

    CMD->CTA = 1 * 1/2 * 1/3
    '''
    if qry_trinuc[1] == 'T':    # ref is 'C'
        ref_char = 'C'
    elif qry_trinuc[1] == 'C':  # ref is 'T'
        ref_char = 'T'
    elif qry_trinuc[1] == 'A':  # ref is 'G'
        ref_char = 'G'
    elif qry_trinuc[1] == 'G':  # ref is 'A'
        ref_char = 'A'
    return nuc_convert_chance(qry_trinuc[0], ref_trinuc[0]) *\
           nuc_convert_chance(ref_char     , ref_trinuc[1]) *\
           nuc_convert_chance(qry_trinuc[2], ref_trinuc[2])

def rc_seq(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

# def rc_seq_chararray(array):
#     complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
#     return np.array([complement.get(base, base) for base in array[::-1]])

# def np_sliding_window(a, window):
#     # https://stackoverflow.com/a/6811241/7859425
#     # https://docs.scipy.org/doc/numpy/reference/generated/numpy.lib.stride_tricks.as_strided.html
#     # create an array mapping, disable writeability for safety
#     shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
#     strides = a.strides + (a.strides[-1],)
#     return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides, writeable=False)

def main():
    ambiguous_base_set = {
        'A': set(['A', 'W', 'M', 'R', 'D', 'H', 'V', 'N']),
        'T': set(['T', 'W', 'K', 'Y', 'B', 'D', 'H', 'N']),
        'C': set(['C', 'S', 'M', 'Y', 'B', 'H', 'V', 'N']),
        'G': set(['G', 'S', 'K', 'R', 'B', 'D', 'V', 'N'])
    }
    trinuc_combinations = [
        'ACA', 'ACT', 'ACC', 'ACG',
        'TCA', 'TCT', 'TCC', 'TCG',
        'CCA', 'CCT', 'CCC', 'CCG',
        'GCA', 'GCT', 'GCC', 'GCG'
    ]

    #########

    aln_fa = sys.argv[1]

    titles = []
    seqs = []
    with open(aln_fa) as f:
        for t, s in SimpleFastaParser(f):
            titles.append(t)
            seqs.append(s.upper())

    # unicode 2d array, dtype='U1'
    # https://stackoverflow.com/a/9493192/7859425
    seq_2d = np.array([list(x) for x in seqs])

    seq_1d_counts, consensus = generate_consensus(seq_2d)
    mu, min_n = cal_mutation_rate(seq_2d, consensus)
    # min_n = 0  # disable filtering by mu

    for i in range(seq_2d.shape[0]):
        seq_ID = titles[i]
        seq_length = np.count_nonzero(seq_2d[i] != '-')

        forward_positive  = defaultdict(lambda: 0)  # C->T, + strand
        backward_positive = defaultdict(lambda: 0)
        forward_negative  = defaultdict(lambda: 0)
        backward_negative = defaultdict(lambda: 0)
        for j in range(1, seq_2d.shape[1] - 2):  # omit first base
            ref_base = consensus[j:j+3]
            qry_base = seq_2d[i,j:j+3]
            # print(i, j, ref_base, qry_base)
            # skip any gaps for now
            if np.any(qry_base == '-'):
                continue
            if ref_base[1] in ambiguous_base_set['C'] and qry_base[1] == 'T' and seq_1d_counts[j+1,1] >= min_n:    # C->T, +
                forward_positive[''.join([qry_base[0], 'C', qry_base[2]])] += \
                    trinuc_rip_chance(ref_base, qry_base)
            elif ref_base[1] in ambiguous_base_set['T'] and qry_base[1] == 'C' and seq_1d_counts[j+1,2] >= min_n:  # T->C, +
                backward_positive[''.join(qry_base)] += \
                    trinuc_rip_chance(ref_base, qry_base)
            elif ref_base[1] in ambiguous_base_set['G'] and qry_base[1] == 'A' and seq_1d_counts[j+1,0] >= min_n:  # G->A, -
                forward_negative[rc_seq(''.join([qry_base[0], 'G', qry_base[2]]))] += \
                    trinuc_rip_chance(ref_base, qry_base)
            elif ref_base[1] in ambiguous_base_set['A'] and qry_base[1] == 'G' and seq_1d_counts[j+1,3] >= min_n:  # A->G, -
                backward_negative[rc_seq(''.join(qry_base))] += \
                    trinuc_rip_chance(ref_base, qry_base)
        summary = []
        for item in trinuc_combinations:
            summary.append(sum([forward_positive[item], backward_positive[item], forward_negative[item], backward_negative[item]])/seq_length)
        # cols: aln_fa, mu, seq_ID, seq_length, ACA, ACT, ACC, ACG, TCA, TCT, TCC, TCG, CCA, CCT, CCC, CCG, GCA, GCT, GCC, GCG
        print('{}\t{}\t{}\t{}\t{}'.format(aln_fa, mu, seq_ID, seq_length, '\t'.join([str(x) for x in summary])))

    # print(forward_positive)
    # print(backward_positive)
    # print(forward_negative)
    # print(backward_negative)
    # summary = []
    # for i in trinuc_combinations:
    #     summary.append(sum([forward_positive[i], backward_positive[i], forward_negative[i], backward_negative[i]]))
    # print('\t'.join([str(x) for x in summary]))

def test():
    aln_fa = sys.argv[1]

    titles = []
    seqs = []
    with open(aln_fa) as f:
        for t, s in SimpleFastaParser(f):
            titles.append(t)
            seqs.append(s.upper())

    # unicode 2d array, dtype='U1'
    # https://stackoverflow.com/a/9493192/7859425
    seq_2d = np.array([list(x) for x in seqs])

    seq_2d_counts, consensus = generate_consensus(seq_2d)
    print(seq_2d_counts, consensus)
    mu, n = cal_mutation_rate(seq_2d, consensus)
    print(mu, n)


if __name__ == '__main__':
    main()
    # test()
