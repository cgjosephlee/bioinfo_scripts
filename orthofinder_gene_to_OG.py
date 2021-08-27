#!/usr/bin/env python3
'''
Map gene ID to OG.

Run in 'OrthoFinder/Results_Jun30/Orthogroups'
Generate 'Genes_to_Orthogroups.tsv'
'''

import re

file_og = 'Orthogroups.tsv'
file_og_ua = 'Orthogroups_UnassignedGenes.tsv'
fout = open('Genes_to_Orthogroups.tsv', 'w')

with open(file_og) as f:
    line = next(f)
    line = line.strip().split()
    names = {}
    for n, name in enumerate(line[1:]):
        names[n] = name
    for line in f:
        line = line.rstrip('\n')
        line = re.split('\t', line)
        assert len(line) == len(names)+1
        OG = line[0]
        for n, genes in enumerate(line[1:]):
            if genes:
                genes = re.split(r',\s', genes)
                for g in genes:
                    print(f'{names[n]}\t{g}\t{OG}', file=fout)
with open(file_og_ua) as f:
    line = next(f)
    for line in f:
        line = line.rstrip('\n')
        line = re.split('\t', line)
        assert len(line) == len(names)+1
        OG = line[0]
        for n, genes in enumerate(line[1:]):
            if genes:
                genes = re.split(r',\s', genes)
                for g in genes:
                    print(f'{names[n]}\t{g}\t{OG}', file=fout)
fout.close()
