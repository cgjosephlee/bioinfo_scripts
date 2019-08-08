#!/usr/bin/env python3

'''
After run `run-asm-pipeline-post-review.sh`, generate agp file according to new fasta and `*.review.assembly`.
'''

import sys
import re

if len(sys.argv) != 3:
    print('\n{} <assembly> <fasta> > stdout\n'.format(sys.argv[0]))
    sys.exit(1)

with open(sys.argv[1]) as ass:
    cprops = []
    asm = []
    fragment = {}
    for line in ass:
        line = line.strip()
        if line.startswith('>'):
            line = line[1:].split(' ')
            cprops.append([line[0], int(line[2])])
            # creat a dict for fragmented sequence
            if re.search('fragment_', line[0]):
                frag_id = line[0].split(':::')[0]
                fragment.setdefault(frag_id, []).append(int(line[2]))
        elif line:
            asm.append(list(int(x) for x in line.split(' ')))

with open(sys.argv[2]) as FA:
    id = []
    for line in FA:
        if line.startswith('>'):
            line = line.strip()
            id.append(line[1:])

def comp_coord(n):
    if n > 0:
        orient = str('+')
        n -= 1
    elif n < 0:
        orient = str('-')
        n = abs(n) - 1
    if re.search('fragment_', cprops[n][0]):
        frag_id = cprops[n][0].split(':::')[0]
        frag_n = int(re.match(r'.+fragment_(\d+)\D*', cprops[n][0]).group(1))
        frag_coord = fragment[frag_id]
        if frag_n == 0:
            return [frag_id, 1, frag_coord[0], orient]
        else:
            return [frag_id, sum(frag_coord[0:frag_n-1])+1, sum(frag_coord[0:frag_n]), orient]
    else:
        return [cprops[n][0], 1, cprops[n][1], orient]

# print .agp
# <obj> <obj_beg> <obj_end> <part_num> <type> <comp_id> <comp_beg> <comp_end> <orientation>
#                                      N      <gap_len> scaffold   yes        na
gap_len = 500
if len(asm) == len(id):
    print('##agp-version 2.0')
    for n in range(len(id)):
        out = []
        for portion in range(len(asm[n])):
            # first portion
            if portion == 0:
                first_end = cprops[abs(asm[n][portion]-1)][1]
                out = [id[n], 1, first_end, portion+1, 'W'] + comp_coord(asm[n][portion])
                print('\t'.join(str(x) for x in out))
            elif portion > 0:
                gap_out = [id[n], first_end+1, first_end+gap_len, 2*portion, 'N', gap_len, 'scaffold', 'yes', 'na']
                out = [id[n], first_end+gap_len+1, first_end+gap_len+cprops[abs(asm[n][portion])-1][1], 2*portion+1, 'W'] + comp_coord(asm[n][portion])
                print('\t'.join(str(x) for x in gap_out))
                print('\t'.join(str(x) for x in out))
                first_end += gap_len+cprops[abs(asm[n][portion])-1][1]
            #print('{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}'.format())
else:
    print('sequence numbers did not match fasta file', file=sys.stderr)
    sys.exit(1)

