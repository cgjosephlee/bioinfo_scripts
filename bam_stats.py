#!/usr/bin/env python3

import argparse
import pysam
from collections import defaultdict

class explain_flag():
    def __init__(self, n):
        assert isinstance(n, int)
        self.n                = n
        self.bit              = bin(n)[2:].zfill(12)  # 12 digits
        self.is_paired        = False  # 1
        self.is_proper_pair   = False  # 2
        self.is_unmapped      = False  # 4
        self.mate_is_unmapped = False  # 8
        self.is_reverse       = False  # 16
        self.mate_is_reverse  = False  # 32
        self.is_read1         = False  # 64
        self.is_read2         = False  # 128
        self.is_secondary     = False  # 256
        self.is_qcfail        = False  # 512
        self.is_duplicate     = False  # 1024
        self.is_supplementary = False  # 2048

        if self.bit[-1] == '1':
            self.is_paired        = True  # 1
        if self.bit[-2] == '1':
            self.is_proper_pair   = True  # 2
        if self.bit[-3] == '1':
            self.is_unmapped      = True  # 4
        if self.bit[-4] == '1':
            self.mate_is_unmapped = True  # 8
        if self.bit[-5] == '1':
            self.is_reverse       = True  # 16
        if self.bit[-6] == '1':
            self.mate_is_reverse  = True  # 32
        if self.bit[-7] == '1':
            self.is_read1         = True  # 64
        if self.bit[-8] == '1':
            self.is_read2         = True  # 128
        if self.bit[-9] == '1':
            self.is_secondary     = True  # 256
        if self.bit[-10] == '1':
            self.is_qcfail        = True  # 512
        if self.bit[-11] == '1':
            self.is_duplicate     = True  # 1024
        if self.bit[-12] == '1':
            self.is_supplementary = True  # 2048

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='bam statistics')
    parser.add_argument('bam',
                        help='bam file')
    parser.add_argument('-s', action='store_true',
                        help='print single line format')
    args = parser.parse_args()

    in_bam = args.bam

    flag_counts = defaultdict(int)
    with pysam.AlignmentFile(in_bam, 'rb') as fin:
        for rec in fin:
            flag_counts[rec.flag] += 1

    total = 0
    secondary = 0
    supplementary = 0
    mapped = 0
    proper_paired = 0
    fw_strand = 0
    rv_strand = 0
    duplicate = 0
    failQC = 0
    for k, v in flag_counts.items():
        f = explain_flag(k)
        if not (f.is_secondary or f.is_supplementary):
            total += v
            if not f.is_unmapped:
                mapped += v
                if f.is_proper_pair:
                    proper_paired += v
                if f.is_reverse:
                    rv_strand += v
                else:
                    fw_strand += v
            if f.is_qcfail:
                failQC += v
            if f.is_duplicate:
                duplicate += v
        else:
            if f.is_secondary:
                secondary += v
            if f.is_supplementary:
                supplementary += v

    p_mapped        = round(mapped / total, 6)
    p_proper_paired = round(proper_paired / total, 6)
    p_fw_strand     = round(fw_strand / total, 6)
    p_rv_strand     = round(rv_strand / total, 6)
    p_duplicate     = round(duplicate / total, 6)
    p_failQC        = round(failQC / total, 6)

    if args.s:
        # total_reads secondary supplementary mapped_reads mapped_ratio proper_paired proper_paired_ratio fw_strand fw_ratio rv_strand rv_ratio duplicated duplicated_ratio failQC failQC_ratio
        print(
            total,
            secondary,
            supplementary,
            mapped, p_mapped,
            proper_paired, p_proper_paired,
            fw_strand, p_fw_strand,
            rv_strand, p_rv_strand,
            duplicate, p_duplicate,
            failQC, p_failQC,
            sep='\t'
        )
    else:
        out = '''\
Total reads:\t{}
Secondary:\t{}
Supplementary:\t{}
Mapped reads:\t{}\t({:.4%})
Proper pairs:\t{}\t({:.4%})
Forward mapped:\t{}\t({:.4%})
Reverse mapped:\t{}\t({:.4%})
Duplicates:\t{}\t({:.4%})
FailQC:\t{}\t({:.4%})
'''.format(
            total,
            secondary,
            supplementary,
            mapped, p_mapped,
            proper_paired, p_proper_paired,
            fw_strand, p_fw_strand,
            rv_strand, p_rv_strand,
            duplicate, p_duplicate,
            failQC, p_failQC
        )

        print(out)
        print()
        print('# total_reads secondary supplementary mapped_reads mapped_ratio proper_paired proper_paired_ratio fw_strand fw_ratio rv_strand rv_ratio duplicated duplicated_ratio failQC failQC_ratio')
        print(
            total,
            secondary,
            supplementary,
            mapped, p_mapped,
            proper_paired, p_proper_paired,
            fw_strand, p_fw_strand,
            rv_strand, p_rv_strand,
            duplicate, p_duplicate,
            failQC, p_failQC,
            sep='\t'
        )
