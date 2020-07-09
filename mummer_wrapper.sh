#!/bin/bash
set -e

DELTA=$1
IDY=${2:=90}
LEN=${3:=1000}
DELTA_flt=${DELTA%.delta}.i${IDY}.l${LEN}.delta

delta-filter -i $IDY -l $LEN $DELTA > $DELTA_flt

# default
# mummerplot --large --png -p $DELTA_flt $DELTA_flt

# sort by length
grep '>' $DELTA_flt | cut -d ' ' -f 1,3 | sed 's/>//g' | sort -k2,2nr -k1,1 | uniq > ${DELTA_flt}.xorder
grep '>' $DELTA_flt | cut -d ' ' -f 2,4 | sort -k2,2nr -k1,1 | uniq > ${DELTA_flt}.yorder
mummerplot --large --png -R ${DELTA_flt}.xorder -Q ${DELTA_flt}.yorder -p $DELTA_flt $DELTA_flt

# clean up
rm -f ${DELTA_flt}.fplot ${DELTA_flt}.rplot ${DELTA_flt}.gp ${DELTA_flt}.xorder ${DELTA_flt}.yorder
