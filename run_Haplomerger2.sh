#!/bin/bash
# /mnt/nas1/ijt/perlscripts_jit/runHaplomerger2.sh
# must be excuted in a sub-folder under HaploMerger folder
set -e

input=$1
thread=16
prefix=out

HM_dir=/home/hhl/local/HaploMerger2_20180603
chainNet_dir=$HM_dir/chainNet_jksrc20100603_ubuntu64bit
lastz_dir=$HM_dir/lastz_1.02.00_unbuntu64bit

export PATH=$PATH:$chainNet_dir:$lastz_dir
ulimit -n 655350

if [[ ! -s ${prefix}.fa.gz ]]; then
    # fasta header may not contain space or special character
    $HM_dir/bin/faDnaPolishing.pl --legalizing --maskShortPortion=1 --noLeadingN --removeShortSeq=1 < $input 2> cleaning.log | awk '/^>/ {print $1} !/^>/ {print $0}' > input.clean.fa
    $HM_dir/winMasker/windowmasker -mk_counts -infmt fasta -sformat obinary -in input.clean.fa -out input.clean.count
    $HM_dir/winMasker/windowmasker -outfmt fasta -ustat input.clean.count -in input.clean.fa -out input.clean.masked.fa
    gzip input.clean.masked.fa
    mv input.clean.masked.fa.gz $prefix.fa.gz
fi

cp $HM_dir/project_template/all_lastz.ctl $HM_dir/project_template/scoreMatrix.q $HM_dir/project_template/hm.batchB* .
sed -i -r "s/threads=[0-9]+/threads=${thread}/" hm.batchB*

bash hm.batchB1.initiation_and_all_lastz            $prefix
bash hm.batchB2.chainNet_and_netToMaf               $prefix
bash hm.batchB3.haplomerger                         $prefix
bash hm.batchB4.refine_unpaired_sequences           $prefix
bash hm.batchB5.merge_paired_and_unpaired_sequences $prefix

# cleaning
rm -rf ${prefix}.${prefix}x.result/raw.axt
