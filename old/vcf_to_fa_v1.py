#!/usr/bin/env python3
'''
v1
Generate aligned fasta for phylogenetic analysis from merged vcf.

MAF?
LD?
'''

import sys
import argparse
import math
import random
from cyvcf2 import VCF
import numpy as np
from tqdm import tqdm

parser = argparse.ArgumentParser(description='Generate aligned fasta for phylogenetic analysis from merged vcf.')
parser.add_argument('vcf',
                    help='vcf file')
parser.add_argument('--missing-rate', type=float, default=0.1,
                    help='maximal missing rate of a site')
parser.add_argument('--min-sample', type=int,
                    help='minimal samples of a site (ignore missing rate if this is given)')
parser.add_argument('--haplotype', choices=[0, 1, 2], default=1,
                    help='representative haplotype if phased')
parser.add_argument('--missing-to-ref', action='store_true',
                    help='assign REF to missing call "./."')
parser.add_argument('--iupac', action='store_true',
                    help='use IUPAC code for HET, otherwise choose at random')
parser.add_argument('--seed',
                    help='random seed')
args = parser.parse_args()

in_vcf = args.vcf
handle = VCF(in_vcf)
'''
gts012=False
var.gt_types = {
    0: 'HOM_REF',
    1: 'HET',
    2: 'UNKNOWN',
    3: 'HOM_ALT'
}
'''

IGNORE_PHASING = False
HAPLOTYPE = None
if args.haplotype == 0:
    IGNORE_PHASING = True
    raise NotImplementedError('just don\'t phase')
elif args.haplotype == 1:
    HAPLOTYPE = 0  # 0: H1, 2: H2, index of string 'A|G'
elif args.haplotype == 2:
    HAPLOTYPE = 2
MISSING_TO_REF = args.missing_to_ref
HET_TO_IUPAC = args.iupac

SAMPLES = handle.samples
NUM_SAMPLE = len(SAMPLES)
MAX_MISSING_RATE = args.missing_rate
if args.min_sample:
    MIN_SAMPLE = args.min_sample
    MAX_MISSING_RATE = None
else:
    MIN_SAMPLE = int(NUM_SAMPLE - math.ceil(NUM_SAMPLE * MAX_MISSING_RATE))

out = '''\
Num. samples: {}
Max. missing rate: {}
Min. samples of a site: {}
HAPLOTYPE: {}
MISSING_TO_REF: {}
HET_TO_IUPAC: {}'''.format(
    NUM_SAMPLE,
    MAX_MISSING_RATE,
    MIN_SAMPLE,
    args.haplotype,
    MISSING_TO_REF,
    HET_TO_IUPAC
)
print(out, file=sys.stderr)

# diploid only
convert_het = {
    'A/C': 'M',
    'C/A': 'M',
    'A/G': 'R',
    'G/A': 'R',
    'A/T': 'W',
    'T/A': 'W',
    'C/G': 'S',
    'G/C': 'S',
    'C/T': 'Y',
    'T/C': 'Y',
    'G/T': 'K',
    'T/G': 'K'
}

random.seed(args.seed)
random_pick_index = [0, 2]  # REF, ALT

print('Start reading...', file=sys.stderr)

SITES = 0  # all sites, including invariant sites
VAR_SITES = 0  # excluding indels
VAR_SITES_PASS = 0
ALL_REF = []
ALL_PHASE = []
ALL_TYPES = []
ALL_BASES = []
for var in handle:
    # DEBUG
    # if SITES == 1e5:
    #     break
    if VAR_SITES_PASS == 5e3:
        break

    SITES += 1
    if SITES % 5e4 == 0:
        print('{} processed.'.format(SITES), file=sys.stderr)
        # print(VAR_SITES, VAR_SITES_PASS, file=sys.stderr)

    # skip REF only site
    if not var.ALT:
        continue
    # skip REF+missing only
    if set(var.gt_types) == set([0, 2]):
        continue
    # skip indel currently
    if var.is_indel:
        continue
    VAR_SITES += 1
    # skip if too many missing value
    if var.num_called < MIN_SAMPLE:
        continue
    VAR_SITES_PASS += 1

    # main logic
    # print(var.CHROM, var.end, var.REF, var.ALT, var.gt_bases, var.gt_types, file=sys.stderr)
    # print(set(var.gt_types), file=sys.stderr)
    ALL_REF.append(var.REF)
    ALL_PHASE.append(var.gt_phases.copy())
    # not really sure what's happening, but this do solve
    # https://stackoverflow.com/q/55019837/7859425
    ALL_TYPES.append(var.gt_types.copy())
    ALL_BASES.append(var.gt_bases.copy())

out = '''\
All sites: {}
Variant sites: {} (excluding indels)
Pass variant sites: {}
Pass rate: {:.4f}'''.format(
    SITES,
    VAR_SITES,
    VAR_SITES_PASS,
    VAR_SITES_PASS / VAR_SITES
)
print(out, file=sys.stderr)

print('Preparing output...', file=sys.stderr)
seq_array = np.empty((VAR_SITES_PASS, NUM_SAMPLE), dtype='U1')
for i in tqdm(range(VAR_SITES_PASS)):
    P = ALL_PHASE[i]
    T = ALL_TYPES[i]
    B = ALL_BASES[i]
    for j in range(NUM_SAMPLE):
        # phased HET
        if P[j]:
            seq_array[i][j] = B[j][HAPLOTYPE]
        # missing
        elif T[j] == 2:
            if MISSING_TO_REF:
                seq_array[i][j] = ALL_REF[i]
            else:
                seq_array[i][j] = 'N'
        # unphased HET
        elif T[j] == 1:
            if HET_TO_IUPAC:
                seq_array[i][j] = convert_het[B[j]]
            else:
                seq_array[i][j] = B[j][random.choice(random_pick_index)]
        # HOMO REF/ALT
        elif T[j] == 0 or T[j] == 3:
            seq_array[i][j] = B[j][0]

print('Now printing...', file=sys.stderr)
# for Seq in seq_array:
#     print(''.join(Seq))
# for Seq in seq_array.T:
#     print(''.join(Seq))
for Name, Seq in zip(SAMPLES, seq_array.T):
    print('>{}\n{}'.format(
        Name,
        ''.join(Seq)
    ))

print('Finish.', file=sys.stderr)
