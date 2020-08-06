#!/usr/bin/env python3
'''
Estimate genome QV.
'''

import argparse
import sys
import tempfile
import subprocess as sp
import gzip
from cyvcf2 import VCF
from math import log

parser = argparse.ArgumentParser(description='Calculate genome QV.')
parser.add_argument('VCF', help='VCF from freebayes or bcftools')
parser.add_argument('CovBed', help='BED from mosdepth')
parser.add_argument('--minDP', type=int, default=10,
                    help='default: 10')
parser.add_argument('--minQ', type=int, default=30,
                    help='default: 30')
parser.add_argument('--region',
                    help='subregion file in bedtools supporting format')
args = parser.parse_args()

inVCF = args.VCF
inBED = args.CovBed

minDP = args.minDP
minQ  = args.minQ

if args.region:
    inREG = args.region
    handle_VCF = tempfile.NamedTemporaryFile(suffix='.vcf')
    handle_BED = tempfile.NamedTemporaryFile(suffix='.bed')
    # print(f'Create {handle_VCF.name} and {handle_BED.name}.', file=sys.stderr)
    cmd = ['bedtools', 'intersect', '-header', '-a', inVCF, '-b', inREG]
    sp.run(cmd, stdout=handle_VCF, check=True)
    handle_VCF.seek(0)
    handle_VCF = VCF(handle_VCF.name)
    cmd = ['bedtools', 'intersect', '-a', inBED, '-b', inREG]
    sp.run(cmd, stdout=handle_BED, check=True)
    handle_BED.seek(0)
else:
    handle_VCF = VCF(inVCF)
    if inBED.endswith('.gz'):
        handle_BED = gzip.open(inBED, 'rt')
    else:
        handle_BED = open(inBED)

# print('Reading vcf...', file=sys.stderr)
# assume single sample
diff_base = 0
num_SNP = 0
num_INDEL = 0
for var in handle_VCF:
    if var.QUAL < minQ:
        continue
    # if var.gt_depths[0] < minDP:  # freebayes
    #     continue
    if var.INFO['DP'] < minDP:  # bcftools
        continue
    if not var.gt_types[0] == 3:  # HOM_ALT (1/1)
        continue

    if var.is_snp:
        diff_base += 1
        num_SNP += 1
    elif var.is_indel:
        # diff_base += abs(len(var.REF) - len(var.ALT))
        diff_base += 1
        num_INDEL += 1
    else:
        continue

# print('Reading bed...', file=sys.stderr)
# handle same-value-compressed bed
total_base = 0
for line in handle_BED:
    line = line.strip().split()
    if int(line[3]) >= minDP:
        total_base += (int(line[2]) - int(line[1]))

handle_VCF.close()
handle_BED.close()

QV = -10 * log(diff_base / total_base, 10)

print('SNP:{}, INDEL:{}, Diff:{}, Total:{}, QV:{:.2f}'.format(num_SNP, num_INDEL, diff_base, total_base, QV))
