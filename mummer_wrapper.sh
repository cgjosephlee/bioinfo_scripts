#!/bin/bash
set -e

DELTA=$1
IDY=${2:=90}
LEN=${3:=1000}
DELTA_flt=${DELTA%.delta}.i${IDY}.l${LEN}.delta

delta-filter -i $IDY -l $LEN $DELTA > $DELTA_flt
mummerplot --large --png -p $DELTA_flt $DELTA_flt
