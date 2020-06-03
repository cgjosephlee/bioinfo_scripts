#!/bin/bash
set -e
# https://github.com/dfguan/purge_dups

pri_asm=$1
reads=$2
thread=16
prefix=out

# python pd_config.py asm.fasta ‘pwd’ <pb folder> <10x folder left blank> asm 
split_fa $pri_asm > $pri_asm.split
minimap2 -t $thread -x asm5 -DP $pri_asm.split $pri_asm.split > $pri_asm.split.self.paf 
minimap2 -t $thread -x map-ont $pri_asm $reads > read.$jobid.paf
pbcstat read.*.paf 
hist_plot.py PB.stat PB.stat.png
calcuts PB.stat > cutoffs 2> calcults.log 
purge_dups -2 -T cutoffs -c PB.base.cov $pri_asm.split.self.paf > dups.bed 2> purge_dups.log 
get_seqs -p $prefix dups.bed $pri_asm.fasta
