#!/usr/bin/env python3
'''
When you change the chr names and don't want to re-run mappings.
'''

import sys
import pysam

LIST_IN = sys.argv[2]

ID_mappings = {}
with open(LIST_IN) as fin:
    for line in fin:
        line = line.strip().split()
        ID_mappings[line[0]] = line[1]

# auto-detect sam/bam
BAMIN = pysam.AlignmentFile(sys.argv[1])

# sam header
header = BAMIN.header.to_dict()
for i in header['SQ']:
    i['SN'] = ID_mappings[i['SN']]

# stdout, may pipe to samtools view -b
BAMOUT = pysam.AlignmentFile('-', 'w', header=header)

# sam records
for rec in BAMIN:
    # RNAME, RNEXT are referenced to ID in header, so they are automatically modified by pysam
    # SA tag
    if rec.has_tag('SA'):
        t = rec.get_tag('SA').split(';')
        nt = ''
        for j in t:
            if j:
                k = j.partition(',')
                j = ID_mappings[k[0]] + ',' + k[2]
            nt = nt + ';' + j
        rec.set_tag('SA', nt, value_type='Z')
    # bwa XA tag
    if rec.has_tag('XA'):
        t = rec.get_tag('XA').split(';')
        nt = ''
        for j in t:
            if j:
                k = j.partition(',')
                j = ID_mappings[k[0]] + ',' + k[2]
            nt = nt + ';' + j
        rec.set_tag('XA', nt, value_type='Z')
    # any other tags? CC?
    BAMOUT.write(rec)

BAMIN.close()
BAMOUT.close()
