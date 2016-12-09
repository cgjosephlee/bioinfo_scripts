#!/usr/bin/env python

import sys
from Bio import SeqIO

IN = open(sys.argv[1],"r")
FA = SeqIO.to_dict(SeqIO.parse(IN, "fasta"))
IN.close()

for ID in FA.keys():
    NAME = ID.replace("|", "_") + ".fasta"
    OUT = open(NAME,"w")
    SeqIO.write(FA[ID], OUT, "fasta")
    OUT.close()
