#!/bin/bash

# Inputs:
# $1: BLAST output file (-outfmt "6 sgi sseqid sseq evalue stitle")
# $2: output filename

# Output:
# file containing the list of unique locus tags from column 2 of $1

# Example run:
# bash get_locus_tags.sh ${HOME}/data/PA3565_orthologs_top_evalueThreshold_1e-50.txt ${HOME}/data/locus_tags_synteny_unfiltered.txt

cut -f 2 $1 | cut -d "|" -f 2 | sort | uniq  > $2