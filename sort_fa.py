#!/usr/bin/env python

import sys
import re

with open(sys.argv[1]) as f:
    header = re.compile("^>")
    seq = {}
    for line in f:
        line = line.strip()
        if header.match(line):
            ID = line[1:]
        else:
            if ID not in seq:
                seq[ID] = line
            else:
                seq[ID] += line

for title in sorted(seq):
    print ">%s\n%s" % (title,seq[title])
