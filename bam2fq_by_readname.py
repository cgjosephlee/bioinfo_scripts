#!/usr/bin/env python3
'''
Aim to extract hardclipped reads in bam file.

Q:
pysam not iterating whole fastq?
https://github.com/pysam-developers/pysam/issues/738
Re-gzip fastq may solve. Or try using zgrep.
'''

import sys
import pysam

bam_in = sys.argv[1]
fq_in = sys.argv[2]

print('Reading bam...', file=sys.stderr)
read_ids = []
with pysam.AlignmentFile(bam_in) as f:
    for rec in f:
        read_ids.append(rec.query_name)

read_ids = set(read_ids)
print('{} unique read ids in bam.'.format(len(read_ids)), file=sys.stderr)
# for k in read_ids:
#     print(k)
# sys.exit(1)

print('Reading fastq...', file=sys.stderr)
with pysam.FastxFile(fq_in) as f:
    for rec in f:
        if rec.name in read_ids:
            # print(rec.name)
            print(str(rec))
            read_ids.discard(rec.name)

if len(read_ids) != 0:
    print('{} reads are not in the fastq file.'.format(len(read_ids)), file=sys.stderr)
    # for k in read_ids:
    #     print(k)
else:
    print('Done!', file=sys.stderr)
