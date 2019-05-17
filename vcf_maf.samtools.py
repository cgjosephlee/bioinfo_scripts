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
        info = line[7].split(';')
        # skip INDELs
        if info[0] == 'INDEL':
            continue

        chr = line[0]
        pos = line[1]
        features = {}
        for item in info:
            if item.startswith('DP'):
                item = item.replace('=', ',')
                item = item.split(',')
                features[item[0]] = [int(x) for x in item[1:]]
        maf = (features['DP4'][0]+features['DP4'][1])/sum(features['DP4'])
        maf = 1 - maf if maf > 0.5 else maf
        print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4f}'.format(chr, pos,
                                                          features['DP'][0],
                                                          features['DP4'][0], features['DP4'][1], features['DP4'][2], features['DP4'][3],
                                                          maf))
