#!/usr/bin/env python3

import sys
import pysam

def get_strand(rec):
    if not rec.is_reverse:
        return '0'
    else:
        return '16' # anything other than 0 is reverse

def get_contact(R1,R2):
    if not R1.query_name == R2.query_name:
        sys.exit()
    if R1.reference_name < R2.reference_name:
        A,B = R2,R1
    else:
        A,B = R1,R2
    # <readname> <str1> <chr1> <pos1> <frag1> <str2> <chr2> <pos2> <frag2> <mapq1> <mapq2>
    out = [A.query_name,
           get_strand(A),
           A.reference_name,
           A.reference_start+1,
           '0', # <frag1>
           get_strand(B),
           B.reference_name,
           B.reference_start+1,
           '1', # <frag2>
           A.mapping_quality,
           B.mapping_quality
           ]
    return out

if len(sys.argv) == 1:
    print('make_hic_contacts.py [s/bam sorted by queryname] (to stdout)')
    sys.exit(0)
elif sys.argv[1].endswith('sam'):
    readtype = 'r'
elif sys.argv[1].endswith('bam'):
    readtype = 'rb'
else:
    print('invalid file type')
    sys.exit(1)

with pysam.AlignmentFile(sys.argv[1], readtype) as infile:
    if not infile.header['HD']['SO'] == 'queryname':
        print('need to sort bam by queryname')
        sys.exit(1)
    readA, readB = '', ''
    contacts = []
    for rec in infile:
        if rec.is_read1 and not readA:
            readA = rec
        elif rec.is_read1:
            contacts.append(get_contact(readA, readB))
            readA = rec
        elif rec.is_read2 and rec.query_name == readA.query_name:
            readB = rec
    contacts.append(get_contact(readA, readB))

contacts = sorted(contacts, key=lambda x: (x[2], x[6]))

for item in contacts:
    print('\t'.join(str(x) for x in item))

