#!/bin/bash

# Hsin-Han Lee, 2015/8
# calculate the SNPs from vcf format file

if [ "$#" -eq 0 ]; then
    echo "Usage: $0 [vcf file(.gz)]"
    exit
else
    # check if input is gzipped
    if [[ ${1} == *gz ]]; then
        >&2 echo -e "input is gzipped\n" # echo to STDERR
        FILE=$(zcat $1 | grep -v '#' | grep -iv 'indel')
    else
        >&2 echo -e "input is not gzipped\n" # echo to STDERR
        FILE=$(cat $1 | grep -v '#' | grep -iv 'indel')
    fi
fi

total=$(wc -l <<< "$FILE") # <<< is here string
cal=$(awk -v total="$total" 'substr($10,1,1) == substr($10,3,1) {homo++ } substr($10,1,1) != substr($10,3,1) {hetero++} END {print homo "\t" hetero "\t" hetero/total }' <<< "$FILE")
# main func, parse the 10th column which looks like "0/1:49,0,152", 1/1 means homozygous, 0/1 means heterozygous
homo=$(cut -f 1 <<< "$cal")
hetero=$(cut -f 2 <<< "$cal")
ratio=$(cut -f 3 <<< "$cal")

echo -e "total SNPs:" $total
echo -e "homozygous alleles:" $homo
echo -e "heterozygous alleles:" $hetero
echo "heterozygosity rate:" $ratio

#awk 'substr($10,1,1) = substr($10,3,1) {print $0}' $1
#awk 'substr($10,1,1) != substr($10,3,1) {print $0}' $1
