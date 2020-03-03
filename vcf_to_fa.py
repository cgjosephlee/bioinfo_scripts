#!/usr/bin/env python3
'''
v2
Generate aligned fasta for phylogenetic analysis from merged vcf.

MAF?
LD?
'''

import sys
import argparse
import math
import random
import gzip
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
if in_vcf.endswith('.gz'):
    handle = gzip.open(in_vcf, 'rt')
else:
    handle = open(in_vcf)

IGNORE_PHASING = False
HAPLOTYPE = None
if args.haplotype == 0:
    IGNORE_PHASING = True
    raise NotImplementedError('Just don\'t phase.')
elif args.haplotype == 1:
    HAPLOTYPE = 0  # 0: H1, 2: H2, index of string 'A|G'
elif args.haplotype == 2:
    HAPLOTYPE = 2
MISSING_TO_REF = args.missing_to_ref
HET_TO_IUPAC = args.iupac

# read header
for line in handle:
    if line.startswith('#CHROM'):
        SAMPLES = line.strip().split('\t')[9:]
        break
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
amb_base = {
    frozenset(('A', 'C')): 'M',
    frozenset(('A', 'G')): 'R',
    frozenset(('A', 'T')): 'W',
    frozenset(('C', 'G')): 'S',
    frozenset(('C', 'T')): 'Y',
    frozenset(('T', 'G')): 'K'
}

random.seed(args.seed)
random_pick_index = [0, 2]  # REF, ALT

print('Start reading...', file=sys.stderr)

valid_base = set(('A', 'T', 'C', 'G'))
base = {
    '0': '',  # REF
    '1': '',
    '2': '',
    '3': ''
}
all_str = []

SITES = 0  # all sites, including invariant sites
VAR_SITES = 0  # excluding indels
VAR_SITES_PASS = 0
for line in handle:
    # DEBUG
    # if SITES == 5e5:
    #     break
    # if VAR_SITES_PASS == 5e5:
    #     break

    SITES += 1
    if SITES % 5e5 == 0:
        print('{} processed.'.format(SITES), file=sys.stderr)
        # print(VAR_SITES, VAR_SITES_PASS, file=sys.stderr)

    line = line.rstrip().split('\t')

    # parse allele, and skip indel, null ALT
    # REF
    if not len(line[3]) == 1:
        continue
    base['0'] = line[3]
    # DEBUG
    # base['1'] = ''
    # base['2'] = ''
    # base['3'] = ''
    # ALT
    if line[4] == '.':
        continue
    elif line[4] in valid_base:
        base['1'] = line[4]
    else:
        for n, b in enumerate(line[4].split(','), 1):
            if not len(b) == 1:
                break
            base[str(n)] = b
        if not len(b) == 1:  # if break in loop
            continue
    VAR_SITES += 1

    # main logic
    tmp_str = ''
    valid = NUM_SAMPLE
    for g in line[9:]:
        # missing
        if g[0] == '.':
            valid -= 1
            if valid < MIN_SAMPLE:
                break
            if MISSING_TO_REF:
                tmp_str += base['0']
            else:
                tmp_str += 'N'
        # HOM
        elif g[0] == g[2]:
            tmp_str += base[g[0]]
        # phased HET
        elif g[1] == '|':
            tmp_str += base[g[HAPLOTYPE]]
        # unphased HET
        elif g[0] != g[2]:
            if HET_TO_IUPAC:
                tmp_str += amb_base[frozenset((base[g[x]] for x in (0, 2)))]
            else:
                tmp_str += base[g[random.choice(random_pick_index)]]
        else:
            print(line, file=sys.stderr)
            raise ValueError('Malformed line.')
    else:  # inner loop did not break
        # DEBUG
        if len(tmp_str) != NUM_SAMPLE:
            print(line, file=sys.stderr)
            print(base, file=sys.stderr)
            print(tmp_str, file=sys.stderr)
            raise ValueError('Sample number error.')
        all_str.append(tmp_str)
        VAR_SITES_PASS += 1
    # continue  # inner loop did break
handle.close()

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

# seq_array = np.empty((VAR_SITES_PASS, NUM_SAMPLE), dtype='U1')

print('Now printing...', file=sys.stderr)
for i in tqdm(range(NUM_SAMPLE)):
    tmp_str = ''
    for j in range(VAR_SITES_PASS):
        tmp_str += all_str[j][i]
    print('>{}\n{}'.format(SAMPLES[i], tmp_str))

print('Finish.', file=sys.stderr)
