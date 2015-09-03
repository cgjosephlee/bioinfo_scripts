#!/usr/bin/env python

"""
Header name filter for fasta files

author: Hisn-han Lee
date: Sep, 2015
"""

import sys
from Bio import SeqIO

usage = """
Usage: seq_name_filter.py [in.fasta] [list.txt]
One name per line in the list
Output to stdout
"""

def Usage():
    print usage
    sys.exit()
    
if len(sys.argv) == 1:
    Usage()

IN1 = open(sys.argv[1], "r")
IN2 = open(sys.argv[2], "r")

names = [name.rstrip("\n") for name in IN2]
IN2.close()
# print names

for seq in SeqIO.parse(IN1, "fasta"):
    if seq.id in names:
    	print ">" + seq.id
    	print seq.seq

IN1.close()
