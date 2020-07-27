#!/bin/bash
set -x

DIR=$1

mv $DIR/reads $DIR/reads_single
mkdir $DIR/reads
single_to_multi_fast5 -i $DIR/reads_single -s $DIR/reads -f $DIR -n 8000 --recursive
