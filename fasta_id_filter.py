#!/usr/bin/env python

"""
Header name filter for fasta files

author: Hsin-han Lee
date: Sep, 2015
"""

import sys
from Bio import SeqIO

usage = """
Usage: seq_name_filter.py [function] <in.fasta> <list.txt>
[function]
    keep (k)    : keep the sequences in the list
    discard (d) : discard the sequences in the list
One name per line in the list
Output to stdout
"""


def Usage():
    print usage
    sys.exit()


def read_list():
    LIST = open(sys.argv[3], "r")
    names = [name.rstrip("\n") for name in LIST]
    LIST.close()
    return names

if len(sys.argv) != 4:
    Usage()

elif sys.argv[1] in ("keep", "k"):
    names = read_list()
    FILE = open(sys.argv[2], "r")
    for seq in SeqIO.parse(FILE, "fasta"):
        if seq.id in names:
            print ">" + seq.id
            print seq.seq
    FILE.close()

elif sys.argv[1] in ("discard", "d"):
    names = read_list()
    FILE = open(sys.argv[2], "r")
    for seq in SeqIO.parse(FILE, "fasta"):
        if seq.id not in names:
            print ">" + seq.id
            print seq.seq
    FILE.close()

else:
    Usage()
