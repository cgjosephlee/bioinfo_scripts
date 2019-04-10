#!/usr/bin/env python3
'''
usage: check_fq_dups.py [fastq]
'''

import sys
import pysam
from collections import Counter
import gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def fastq_name_pysam(infile):
    '''
    stoped before the malformed record and not return error
    '''
    names = []
    # faster than list comprehension
    with pysam.FastxFile(infile) as f:
        for entry in f:
            # print(entry)
            names.append(entry.name)
    return names

def fastq_name_biopy(infile):
    '''
    return error if record is malformed
    return full title line
    '''
    names = []
    # handle gzipped file
    if infile.endswith('.gz'):
        f = gzip.open(infile, 'rt')
    else:
        f = open(infile)
    try:
        for entry in FastqGeneralIterator(f):
            # title, seq, qual
            # print(entry[0])
            names.append(entry[0].split()[0])
        f.close()
    except ValueError as e:
        print('Found malformed record in {}!'.format(IN))
        print(e)
        f.close()
        sys.exit(1)
    return names

def deduplicate_list(inlist):
    # faster than for loop
    dups = [[item, count] for item, count in Counter(inlist).items() if count > 1]
    return dups

IN = sys.argv[1]

# names = fastq_name_pysam(IN)
names = fastq_name_biopy(IN)
dups = deduplicate_list(names)

print('Found {} sequences in {}.'.format(len(names), IN))
if len(dups) != 0:
    print('Found dupclicated records!')
    for i in dups:
        print('{}\t{}'.format(i[1], i[0]))
        sys.exit(1)
else:
    print('No duplicate was found!')

