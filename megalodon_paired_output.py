#!/usr/bin/env python3
'''
Merge +/- pairs of a motif
'''

import sys

fin = sys.argv[1]
minDepth = 4  # of a motif

with open(fin) as f:
    prev_line = None
    for line in f:
        line = line.strip()
        if prev_line:
            l1 = prev_line.split()
            l2 = line.split()
            if l1[0] == l2[0] and l1[2] == l2[1]:  # is a pair
                if int(l1[9]) + int(l2[9]) >= minDepth:  # add up depths of + and -
                    out = l1
                    out[7] = l2[2]
                    out.append(l2[9])
                    out.append(l2[10])
                    out.append(f'{int(l1[9])+int(l2[9])}')
                    out.append(f'{((float(l1[9])*float(l1[10]))+(float(l2[9])*float(l2[10])))/(float(l1[9])+float(l2[9])):.1f}')
                    out[4] = out[-2]
                    print('\t'.join(out))
                prev_line = None
            else:
                prev_line = line
        else:
            prev_line = line
