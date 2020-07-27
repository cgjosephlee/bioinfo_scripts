#!/bin/bash
set -e

guppy_ver=4.0.11
guppy_model=dna_r9.4.1_450bps_hac.cfg

DIR=$1
IN=$1/reads
OUT=$1/basecall_guppy_$guppy_ver
minLen=1000
# FQ=$1/$DIR.fq.gz
# FQ_flt=$1/$DIR.min${minLen}.fq.gz

KIT=$2

# singularity exec --nv --bind /data:/data --bind /mnt/nas1:/mnt/nas1 /mnt/nas1/hhl/vm/guppy_${guppy_ver}.simg \
#     guppy_basecaller -c $guppy_model -i $IN -s $OUT -x 'cuda:all' -q 0 -r --barcode_kits "${KIT}" --trim_barcodes
for i in `find $OUT -name 'barcode*' -type d`; do
    b=`basename $i`
    echo $b
    cat $i/*fastq | /mnt/nas1/hhl/bin/My_scripts/fastn_length_filter.py --fq -c $minLen - | pigz --best -p 6 > $OUT/${DIR}_${b}.min${minLen}.fq.gz
done

./run_guppy_barcode_summary.py $OUT/sequencing_summary.txt $minLen
