#!/bin/bash

# Hsin-Han Lee

if [ "$#" -eq 0 ]; then
    echo "Usage: $0 [fasta file] [length cutoff]"
    echo "Tips: Set length cutoff to 0 will only make multipe lines fasta to one line"
    echo "      Default output to stdout"
    exit

elif [[ $1 != @(*.fa|*.fasta) ]] ; then
    echo "ERROR: invalid file"
    exit

elif [ -z "$2" ] || ! [[ "$2" =~ ^-?[0-9]+$ ]]; then
    echo "ERROR: invalid length cutoff"
    exit

else
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $1 | awk -v LEN="$2" '!/^>/ { next } { getline seq } length(seq) >= LEN { print $0 "\n" seq }'

fi
