#!/usr/bin/env python3

# followd by falcon_phase, separate phased diploid fastas to different files
# execute in falcon run directory

import sys
import os
from glob import glob
from Bio import SeqIO

falcon_fa = '2-asm-falcon/p_ctg.fa'
phase_prefix = 'test2'
phased_fa = '5-phase/' + phase_prefix + '.diploid_phased.fasta'

phased_ctg = [os.path.basename(x) for x in glob('3-unzip/0-phasing/*')]

phased_FA = SeqIO.to_dict(SeqIO.parse(phased_fa,'fasta'))
phase0_FA = open('5-phase/'+phase_prefix+'.phase0.fasta', 'w')
phase1_FA = open('5-phase/'+phase_prefix+'.phase1.fasta', 'w')
hp_phase0_FA = open('5-phase/'+phase_prefix+'.haploid_phase0.fasta', 'w')
hp_phase1_FA = open('5-phase/'+phase_prefix+'.haploid_phase1.fasta', 'w')

for rec in SeqIO.parse(falcon_fa,'fasta'):
    p_id = rec.id.split(' ')[0]
    if p_id not in phased_ctg:
        print('>{}\n{}'.format(p_id, rec.seq), file = hp_phase0_FA)
        print('>{}\n{}'.format(p_id, rec.seq), file = hp_phase1_FA)
    else:
        p_id_0 = p_id + '_0'
        p_id_1 = p_id + '_1'
        print('>{}\n{}'.format(p_id_0, phased_FA[p_id_0].seq), file = phase0_FA)
        print('>{}\n{}'.format(p_id_1, phased_FA[p_id_1].seq), file = phase1_FA)
        print('>{}\n{}'.format(p_id_0, phased_FA[p_id_0].seq), file = hp_phase0_FA)
        print('>{}\n{}'.format(p_id_1, phased_FA[p_id_1].seq), file = hp_phase1_FA)

phase0_FA.close()
phase1_FA.close()
hp_phase0_FA.close()
hp_phase1_FA.close()

