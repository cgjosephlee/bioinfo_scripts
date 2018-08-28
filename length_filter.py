#!/usr/bin/env python3

"""
author: Hsinhan Lee
date:   Sep, 2015
update: Aug, 2018
"""

import sys
from Bio import SeqIO

usage = """
Discard the sequences shorter than the cutoff
Usage: length_filter.py <in.fasta> <length cutoff>
Output to stdout
"""

if len(sys.argv) == 1:
    print(usage)
    sys.exit(1)

with open(sys.argv[1], "r") as IN:
    for seq in SeqIO.parse(IN, "fasta"):
        if len(seq) > int(sys.argv[2]):
            print('>{}\n{}'.format(seq.id, seq.seq))
