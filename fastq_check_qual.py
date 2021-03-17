#!/usr/bin/env python3
'''
HiSeq X Ten and NovaSeq: 2, 12, 23 and 37 ('#', ',', ':', 'F')

from wiki (https://en.wikipedia.org/wiki/FASTQ_format#Encoding):
  SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
  ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................
  ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
  .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ.....................
  LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL....................................................
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
  !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
  |                         |    |        |                              |                     |
 33                        59   64       73                            104                   126
  0........................26...31.......40
                           -5....0........9.............................40
                                 0........9.............................40
                                    3.....9..............................41
  0.2......................26...31........41
  0..................20........30........40........50..........................................93
 S - Sanger        Phred+33,  raw reads typically (0, 40)
 X - Solexa        Solexa+64, raw reads typically (-5, 40)
 I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
 J - Illumina 1.5+ Phred+64,  raw reads typically (3, 41)
 L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)
 P - PacBio        Phred+33,  HiFi reads typically (0, 93)
'''

import sys
import gzip

max = 10000

infile = sys.argv[1]

if infile.endswith('.gz'):
    handle = gzip.open(infile, 'rt')
else:
    handle = open(infile)

quals = set()
n = 0

next(handle)
next(handle)
next(handle)
for line in handle:
    quals.update(list(line.strip()))  # quality line
    next(handle)
    next(handle)
    next(handle)

    n += 1
    if n == max:
        break

handle.close()

q = ''.join(sorted([x for x in quals]))
print(len(q))
print(q)
