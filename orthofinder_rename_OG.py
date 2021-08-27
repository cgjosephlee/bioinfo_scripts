#!/usr/bin/env python3
'''
Rename gene ID to species ID in OG fastas.

Run in 'OrthoFinder/Results_Jun30'
Optional 'WorkingDirectory/SpeciesIDs_rename.txt'
Generate 'fastas.xxxxx'
'''

import os
from glob import glob

dir_in = 'Single_Copy_Orthologue_Sequences'
dir_out = f'fastas.{os.getpid()}'
os.mkdir(dir_out)

SeqID = 'WorkingDirectory/SequenceIDs.txt'

SpeID = 'WorkingDirectory/SpeciesIDs_rename.txt'
if not os.path.exists(SpeID):
    SpeID = 'WorkingDirectory/SpeciesIDs.txt'

SpeID_dict = {}
with open(SpeID) as f:
    for line in f:
        line = line.strip().split()
        SpeID_dict[line[0][:-1]] = line[1]

SeqID_dict = {}
with open(SeqID) as f:
    for line in f:
        line = line.strip().split()
        SeqID_dict[line[1]] = SpeID_dict[line[0].split('_')[0]]

fa_list = glob(dir_in+'/*fa')
for i in fa_list:
    fin = open(i)
    fout = open(os.path.join(dir_out, os.path.basename(i)), 'w')
    for line in fin:
        line = line.strip()
        if line.startswith('>'):
            print('>'+SeqID_dict[line[1:]], file=fout)
        else:
            print(line, file=fout)
    fin.close()
    fout.close()
