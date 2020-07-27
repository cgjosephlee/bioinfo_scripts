#!/bin/bash
set -e

guppy_ver=4.0.11
guppy_model=dna_r9.4.1_450bps_hac.cfg

DIR=$1
IN=$1/reads
OUT=$1/basecall_guppy_$guppy_ver
minLen=1000
# FQ=$1/$DIR.fq.gz
FQ_flt=$1/$DIR.min${minLen}.fq.gz

singularity exec --nv --bind /data:/data --bind /mnt/nas1:/mnt/nas1 /mnt/nas1/hhl/vm/guppy_${guppy_ver}.simg \
    guppy_basecaller -c $guppy_model -i $IN -s $OUT -x 'cuda:all' -q 0 -r
# find $OUT -name '*fastq' | xargs cat | pigz -p 6 > $FQ
find $OUT -name '*fastq' | \
    xargs cat | \
    /mnt/nas1/hhl/bin/My_scripts/fastn_length_filter.py --fq -c $minLen - | \
    pigz --best -p 6 \
    > $FQ_flt
