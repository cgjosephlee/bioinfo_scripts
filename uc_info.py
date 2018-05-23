#!/usr/bin/env python3

import sys

seqs = 0
clus = 0
sings = 0
clu_size = []

with open(sys.argv[1]) as fin:
    for line in fin.readlines():
        if line.startswith('C'):
            line = line.strip()
            line = line.split('\t')
            seqs += int(line[2])
            clus += 1
            clu_size.append(int(line[2]))
            if line[2] == '1':
                sings += 1

print('''
      Seqs  {}
  Clusters  {}
  Max size  {}
  Avg size  {:.1f}
  Min size  {}
Singletons  {}, {:.1%} of seqs, {:.1%} of clusters\n'''.format(seqs,
                                                               clus,
                                                               max(clu_size),
                                                               sum(clu_size)/clus,
                                                               min(clu_size),
                                                               sings,
                                                               sings/seqs,
                                                               sings/clus))
