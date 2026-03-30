#!/bin/bash

# For use with run_local_blast.sh.
# Input: BLAST file with -outfmt "6 sallgi sallseqid sseq evalue salltitles",
# in which the sallseqid column (col 2) is formatted as "genomeID-proteinID".
# This is the case if using a BLAST database custom-built from a genome database using the provided script (TODO)
# (May have additional columns, e.g. organism column, to the right.)

# Output: overwrite input so that col 1 contains genome IDs, and col 2 contains protein IDs

blast=$1
temp=$(mktemp)

awk -F'\t' 'BEGIN { OFS="\t" } {split($2, arr, "-"); print arr[1], arr[2], $0}' "$blast" | cut -f 1,2,5- > "$temp" 
# separate 1st and 2nd column, append the entire file, and use cut to remove the old cols 1 and 2 (which are temporarily cols 3 and 4)

mv "$temp" "$blast"