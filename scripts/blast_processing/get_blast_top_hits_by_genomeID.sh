#!/bin/bash

# CLIs:
# $1: blast file to filter; output of blast2gen.py which adds metadata. Assumes evalues are in column 4, and genome IDs (assembly accessions) are in column 9.
# $2: output file name. To this file, writes only the top hit (i.e. lowest evalue) per genome ID.

# example run:
# bash get_blast_top_hits_by_genomeID.sh "${HOME}/data/blast_outputs/fha1_genome_info_paOnly.tsv" "${HOME}/data/blast_outputs/fha1_genome_info_paOnly_top.tsv"


# The first sort is to ensure that the file is sorted by genome ID and then by evalue (because the input file may not be in order),
# the second is to get the top hit for each genome ID.
# sort -t$'\t' -k9,9 -k4,4g "$1" | sort -u -t$'\t' -k9,9 > "$2"

# wait, I need to get rid of the header
tail -n +2 "$1" | sort -t$'\t' -k9,9 -k4,4g | sort -u -t$'\t' -k9,9 > "$2"