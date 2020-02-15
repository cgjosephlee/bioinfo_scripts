#!/usr/bin/env python3

import sys
from cyvcf2 import VCF

class PhasedBlock:
    def __init__(self):
        self.chr = None
        self.id = None
        self.uid = None
        self.start = None
        self.end = None
        self.len = None
        self.sites = []
        self.n_sites = 0

    def add(self, var):  # Variant class
        if not self.chr:
            self.chr = var.CHROM
        if not self.id:
            self.id = var.format('PS')[0][0]  # int
            self.uid = '{}:{}'.format(self.chr, self.id)
        if not self.start:
            self.start = var.POS
        self.end = var.POS
        self.len = self.end - self.start + 1
        self.sites.append(var.POS)
        self.n_sites += 1

in_VCF = sys.argv[1]
VCF_handle = VCF(in_VCF)

# intialize
chr_stats = {}
try:
    for i, j in zip(VCF_handle.seqnames, VCF_handle.seqlens):
        chr_stats[i] = {}
        chr_stats[i]['len'] = j
        chr_stats[i]['n_variants'] = 0
        chr_stats[i]['n_snps'] = 0
        chr_stats[i]['n_phased_snps'] = 0
        chr_stats[i]['n_indel'] = 0
        chr_stats[i]['n_phased_indel'] = 0
        chr_stats[i]['blocks'] = {}
except AttributeError:
    print('No chr info in vcf header, read from fasta index.', file=sys.stderr)
    in_fai = sys.argv[2]
    with open(in_fai) as f:
        for line in f:
            line = line.split()
            i = line[0]
            j = int(line[1])
            chr_stats[i] = {}
            chr_stats[i]['len'] = j
            chr_stats[i]['n_variants'] = 0
            chr_stats[i]['n_snps'] = 0
            chr_stats[i]['n_phased_snps'] = 0
            chr_stats[i]['n_indel'] = 0
            chr_stats[i]['n_phased_indel'] = 0
            chr_stats[i]['blocks'] = {}

# assuming only one sample
for var in VCF_handle:
    chr = var.CHROM
    chr_stats[chr]['n_variants'] += 1
    if var.is_snp:
        chr_stats[chr]['n_snps'] += 1
    if var.is_indel:
        chr_stats[chr]['n_indel'] += 1
    if var.gt_phases[0]:  # is_phased
        if var.is_snp:
            chr_stats[chr]['n_phased_snps'] += 1
        if var.is_indel:
            chr_stats[chr]['n_phased_indel'] += 1
        block_id = var.format('PS')[0][0]  # int
        chr_stats[chr]['blocks'].setdefault(block_id, PhasedBlock()).add(var)
VCF_handle.close()

for k, v in chr_stats.items():
    b = v['blocks'].values()
    n_blocks = len(b)
    if n_blocks > 0:
        block_sizes = [x.len for x in b]
        max_block_size = max(block_sizes)
        min_block_size = min(block_sizes)
        sum_block_size = sum(block_sizes)
    else:
        block_sizes = []
        max_block_size = 0
        min_block_size = 0
        sum_block_size = 0
    # chr chr_len n_blocks max_block min_block block_total block_perc n_variants n_SNP n_phased_SNP n_INDEL n_phased_INDEL
    print(
        k,
        v['len'],
        n_blocks,
        max_block_size,
        min_block_size,
        sum(block_sizes),
        sum_block_size / v['len'],  # may be interleaved
        v['n_variants'],
        v['n_snps'],
        v['n_phased_snps'],
        v['n_indel'],
        v['n_phased_indel'],
        sep='\t'
    )
