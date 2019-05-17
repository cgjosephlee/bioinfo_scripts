#!/usr/bin/env python3

# use for samtools or bcftoolfs pipelines

import sys
import gzip

IN = sys.argv[1]

if IN.endswith('gz'):
    f = gzip.open(IN, 'rt')
else:
    f = open(IN)

with f:
    print('CHR\tPOS\tDP\tDP4_1\tDP4_2\tDP4_3\tDP4_4\tMAF')
    for line in f.readlines():
        line = line.strip()
        if line.startswith('#'):
            continue
        line = line.split('\t')
        info = line[9].split(':')  # FORMAT column

        # skip INDELs
        # if info[0] == 'INDEL':
        #     continue

        # skip homo
        if info[0] != '0/1':
            continue

        # filter low genotype quality (phred score?)
        if int(info[1]) < 30:
            continue

        # chr = line[0]
        # pos = line[1]

        maf = float(info[6].rstrip('%'))
        maf = 100 - maf if maf > 50 else maf
        print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2f}'.format(line[0], line[1],
                                                          info[2],  # SDP
                                                          info[10], info[11], info[12], info[13],
                                                          maf))
