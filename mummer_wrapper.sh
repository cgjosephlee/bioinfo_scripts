#!/bin/bash
set -e

IDY=90
LEN=1000
SORTBY='d'
# d: default, do nothing
# l: sort by length
# n: sort by name
# N: sort by name in version order

while getopts "i:l:b:" opt; do
    case $opt in
        i) IDY=$OPTARG ;;
        l) LEN=$OPTARG ;;
        b) SORTBY=$OPTARG ;;
    esac
done

shift $((OPTIND-1))
DELTA=$1
DELTA_flt=${DELTA%.delta}.i${IDY}.l${LEN}.delta

delta-filter -i $IDY -l $LEN $DELTA > $DELTA_flt

if [[ "$SORTBY" == "d" ]]; then
    # default
    mummerplot --large --png -p $DELTA_flt $DELTA_flt
elif [[ "$SORTBY" == "l" ]]; then
    # sort by length
    grep '>' $DELTA_flt | cut -d ' ' -f 1,3 | sed 's/>//g' | sort -k2,2nr -k1,1 | uniq > ${DELTA_flt}.xorder
    grep '>' $DELTA_flt | cut -d ' ' -f 2,4 | sort -k2,2nr -k1,1 | uniq > ${DELTA_flt}.yorder
    mummerplot --large --png -R ${DELTA_flt}.xorder -Q ${DELTA_flt}.yorder -p $DELTA_flt $DELTA_flt
elif [[ "$SORTBY" == "n" ]]; then
    grep '>' $DELTA_flt | cut -d ' ' -f 1,3 | sed 's/>//g' | sort -k1,1 | uniq > ${DELTA_flt}.xorder
    grep '>' $DELTA_flt | cut -d ' ' -f 2,4 | sort -k1,1 | uniq > ${DELTA_flt}.yorder
    mummerplot --large --png -R ${DELTA_flt}.xorder -Q ${DELTA_flt}.yorder -p $DELTA_flt $DELTA_flt
elif [[ "$SORTBY" == "N" ]]; then
    grep '>' $DELTA_flt | cut -d ' ' -f 1,3 | sed 's/>//g' | sort -k1,1V | uniq > ${DELTA_flt}.xorder
    grep '>' $DELTA_flt | cut -d ' ' -f 2,4 | sort -k1,1V | uniq > ${DELTA_flt}.yorder
    mummerplot --large --png -R ${DELTA_flt}.xorder -Q ${DELTA_flt}.yorder -p $DELTA_flt $DELTA_flt
fi

# clean up
rm -f ${DELTA_flt}.fplot ${DELTA_flt}.rplot ${DELTA_flt}.gp ${DELTA_flt}.xorder ${DELTA_flt}.yorder
