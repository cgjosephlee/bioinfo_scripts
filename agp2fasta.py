#!/usr/bin/env python3
'''
Generate new fasta from old fasta and agp file.
'''

import sys
from Bio import SeqIO
from collections import defaultdict

in_fa = sys.argv[1]
in_agp = sys.argv[2]
gap_char = 'N'
gap_default_len = 100

sequences = SeqIO.to_dict(SeqIO.parse(in_fa, 'fasta'))

# generate an agp template and quit
if in_agp == '-':
    print('## AGP-version 2.0')
    for k, v in sequences.items():
        out = [k, '1', str(len(v)), '1', 'W', k, '1', str(len(v)), '+']
        print('\t'.join(out))
    sys.exit(0)

agp_arrange = defaultdict(list)
with open(in_agp) as f:
    for line in f:
        if line.startswith('#'):
            continue
        line = line.rstrip().split()
        if line[4] == 'W':
            agp_arrange[line[0]].append([line[5],
                                        int(line[6]) - 1,
                                        int(line[7]),
                                        line[8]])
        elif line[4] == 'N':
            # try to avoid possible ID in fasta
            agp_arrange[line[0]].append(['PLACE_GAP', 0, int(line[5]), None])
        elif line[4] == 'U':  # gap of unknown size
            agp_arrange[line[0]].append(['PLACE_GAP', 0, gap_default_len, None])
        else:
            raise NotImplementedError('Not support yet.')

old_ids = sequences.keys()

for chr in agp_arrange.keys():  # assuming order is kept
    seq = ''
    for item in agp_arrange[chr]:
        if item[0] in old_ids:
            if item[3] == '+':
                seq += sequences[item[0]][item[1]:item[2]].seq
            else:
                seq += sequences[item[0]][item[1]:item[2]].reverse_complement().seq
        elif item[0] == 'PLACE_GAP':
            seq += gap_char * item[2]
        else:
            raise KeyError('{} is not found in reference fasta.'.format(item[0]))
    print('>{}\n{}'.format(chr, seq))
    seq = ''
