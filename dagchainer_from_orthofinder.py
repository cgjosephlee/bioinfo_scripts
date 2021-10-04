#!/usr/bin/env python3
'''
Run in "OrthoFinder/Results_Jun25/WorkingDirectory"
Require "SpeciesBED.txt"
    0: /path/to/bed
Optional "BlastList.txt"
Require followings in $PATH:
    filter_repetitive_matches.pl
    run_DAG_chainer.pl
    parallel
'''

import sys
import os
import gzip
from glob import glob
import subprocess as sp

useRank = True  # unrank is not finished
proc = 1  # run dagchainer with parallel
# careful to use more processes, IO error may happen

##### Step 1

File_SequenceID = 'SequenceIDs.txt'
# if os.path.isfile('SpeciesIDs_rename.txt'):
#     File_SequenceID = 'SpeciesIDs_rename.txt'
print(f'Read "{File_SequenceID}".', file=sys.stderr)
OF_ID_dict = {}
with open(File_SequenceID) as f:
    for line in f:
        line = line.strip().split()
        OF_ID_dict[line[0][:-1]] = line[1]

File_SpeciesBED = 'SpeciesBED.txt'
print(f'Read "{File_SpeciesBED}".', file=sys.stderr)
bed_list = []
with open(File_SpeciesBED) as f:
    for line in f:
        line = line.strip().split()
        if os.path.isfile(line[1]):
            bed_list.append(line[1])
        else:
            print(f'File not found: {line[1]}', file=sys.stderr)

# assume IDs are unique across all files
bed_dict = {}
n = 1
for b in bed_list:
    with open(b) as f:
        for line in f:
            line = line.strip().split()
            if line[3] == '+':
                bed_dict[line[4]] = (line[0], int(line[1])+1, int(line[2]), n)
            else:
                bed_dict[line[4]] = (line[0], int(line[2]), int(line[1])+1, n)
            n += 1
    # rank is continued

# read "BlastList.txt" list if exist
if os.path.isfile('BlastList.txt'):
    print('Read "BlastList.txt".', file=sys.stderr)
    blast_list = []
    with open('BlastList.txt') as f:
        for line in f:
            blast_list.append(line.strip())
else:
    blast_list = glob('Blast*.txt.gz')

for b in blast_list:
    qry, sbj = b[5:-7].split('_')
    if qry == sbj:
        continue
    print(b, file=sys.stderr)
    fin = gzip.open(b, 'rt')
    fout = open(f'DAGchainer{qry}_{sbj}.input', 'w')
    for line in fin:
        line = line.strip().split()
        id0 = OF_ID_dict[line[0]]
        id1 = OF_ID_dict[line[1]]
        id0_info = bed_dict[id0]
        id1_info = bed_dict[id1]
        evalue = line[10]
        if useRank:
            print(f'{id0_info[0]}\t{id0}\t{id0_info[3]}\t{id0_info[3]}\t{id1_info[0]}\t{id1}\t{id1_info[3]}\t{id1_info[3]}\t{evalue}', file=fout)
        else:
            print(f'{id0_info[0]}\t{id0}\t{id0_info[1]}\t{id0_info[2]}\t{id1_info[0]}\t{id1}\t{id1_info[1]}\t{id1_info[2]}\t{evalue}', file=fout)
    fin.close()
    fout.close()

##### Step 2

if useRank:
    # dagchainer ranked
    cmd = f"""ls DAGchainer*input | parallel -j {proc} 'filter_repetitive_matches.pl 5 < {{}} > {{}}Flt'"""
    print(cmd, file=sys.stderr)
    sp.run(cmd, shell=True, check=True)
    cmd = f"""ls DAGchainer*inputFlt | parallel -j {proc} 'run_DAG_chainer.pl -i {{}} -Z 12 -D 10 -g 1 -A 5' > dagchainer.log"""
    print(cmd, file=sys.stderr)
    sp.run(cmd, shell=True, check=True)
else:
    # TODO
    raise NotImplementedError

# output: daginput_filtered.aligncoords

##### Step 3

dag_list = glob('DAGchainer*inputFlt.aligncoords')

for d in dag_list:
    print(d, file=sys.stderr)
    fin = open(d)
    fout = open(f'{d[:-12]}.synblocks', 'w')
    if useRank:
        try:
            next(fin)
        except StopIteration:
            print('Empty file.', file=sys.stderr)
            fin.close()
            fout.close()
            continue
        firstline = next(fin)
        lastline = ''
        for line in fin:
            if line.startswith('##'):
                firstline = firstline.split()
                lastline = lastline.split()
                # subject
                GeneStB = bed_dict[firstline[1]]
                GeneEdB = bed_dict[lastline[1]]
                # query
                GeneStA = bed_dict[firstline[5]]
                GeneEdA = bed_dict[lastline[5]]
                if int(firstline[6]) < int(lastline[6]):  # ordinal query
                    print(f'{GeneStB[0]}\t{min(GeneStB[1:3])}\t{max(GeneEdB[1:3])}\t'
                          f'{GeneStA[0]}\t{min(GeneStA[1:3])}\t{max(GeneEdA[1:3])}',
                          file=fout)
                else:  # reversed query
                    print(f'{GeneStB[0]}\t{min(GeneStB[1:3])}\t{max(GeneEdB[1:3])}\t'
                          f'{GeneStA[0]}\t{max(GeneStA[1:3])}\t{min(GeneEdA[1:3])}',
                          file=fout)
                firstline = next(fin)
            else:
                lastline = line
        # last record
        firstline = firstline.split()
        lastline = lastline.split()
        # subject
        GeneStB = bed_dict[firstline[1]]
        GeneEdB = bed_dict[lastline[1]]
        # query
        GeneStA = bed_dict[firstline[5]]
        GeneEdA = bed_dict[lastline[5]]
        if int(firstline[6]) < int(lastline[6]):  # ordinal query
            print(f'{GeneStB[0]}\t{min(GeneStB[1:3])}\t{max(GeneEdB[1:3])}\t'
                  f'{GeneStA[0]}\t{min(GeneStA[1:3])}\t{max(GeneEdA[1:3])}',
                  file=fout)
        else:  # reversed query
            print(f'{GeneStB[0]}\t{min(GeneStB[1:3])}\t{max(GeneEdB[1:3])}\t'
                  f'{GeneStA[0]}\t{max(GeneStA[1:3])}\t{min(GeneEdA[1:3])}',
                  file=fout)
    else:
        # TODO
        raise NotImplementedError
    fin.close()
    fout.close()
