#!/bin/bash
set -e

guppy_ver=5.0.13
model=r9.4
barcode_kit=
minLen=1000
minQ=10

while getopts "m:b:l:q:" opt; do
    case $opt in
        m) model=$OPTARG ;;
        b) barcode_kit=$OPTARG ;;
        l) minLen=$OPTARG ;;
        q) minQ=$OPTARG ;;
    esac
done

shift $((OPTIND-1))

if [[ "$model" == "r9.4" ]]; then
    guppy_model=dna_r9.4.1_450bps_sup.cfg
elif [[ "$model" == "r9.5" ]]; then
    guppy_model=dna_r9.5_450bps.cfg
elif [[ "$model" == "r10" ]]; then
    guppy_model=dna_r10_450bps_hac.cfg
elif [[ "$model" == "r10.3" ]]; then
    guppy_model=dna_r10.3_450bps_sup.cfg
elif [[ "$model" == "RNA" ]]; then
    guppy_model=rna_r9.4.1_70bps_hac.cfg
else
    guppy_model=$model
fi

barcode_mode=false
if [[ ! -z "$barcode_kit" ]]; then
    barcode_mode=true
    barcode_opts="--barcode_kits ${barcode_kit} --trim_barcodes --chunks_per_runner 150"
    # reduce chunks to use less GPU mem (for GTX1080 Ti 12G)
fi

DIR=$1
IN=$DIR/reads
OUT=$DIR/basecall_guppy_$guppy_ver
FQ_flt=$DIR/$DIR.min${minLen}.Q${minQ}.fq.gz
threads=6  # pigz

# singularity exec --nv --bind /data:/data --bind /mnt/nas1:/mnt/nas1 /mnt/nas1/hhl/vm/guppy_${guppy_ver}.simg \
#     guppy_basecaller -c $guppy_model -i $IN -s $OUT -x 'cuda:all' -q 0 -r --disable_qscore_filtering "$barcode_opts"

# built-in
guppy_basecaller -c $guppy_model -i $IN -s $OUT -x 'cuda:all' -q 0 -r --disable_qscore_filtering $barcode_opts

if [[ "$barcode_mode" == false ]]; then
    find $OUT -name '*fastq' | \
        xargs cat | \
        nanoq -l $minLen -q $minQ | \
        pigz --best -p $threads \
        > $FQ_flt
else
    for i in `find $OUT -name 'barcode*' -type d`; do
        b=`basename $i`
        echo $b
        FQ_flt=$DIR/${DIR}_${b}.min${minLen}.Q${minQ}.fq.gz
        cat $i/*fastq | \
            nanoq -l $minLen -q $minQ | \
            pigz --best -p $threads \
            > $FQ_flt
    done
fi
