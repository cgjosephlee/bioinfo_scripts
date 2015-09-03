#!/usr/bin/env python

"""
Length filter for fasta files

author: Hisn-han Lee
date: Sep, 2015
"""

import sys
from Bio import SeqIO

usage = """
Usage: length_filter.py [in.fasta] [length cutoff]
Output to stdout
"""

def Usage():
    print usage
    sys.exit()
    
if len(sys.argv) == 1:
    Usage()

IN1 = open(sys.argv[1], "r")

for seq in SeqIO.parse(IN1, "fasta"):
    if len(seq) > int(sys.argv[2]):
        print ">" + seq.id
        print seq.seq

IN1.close()
