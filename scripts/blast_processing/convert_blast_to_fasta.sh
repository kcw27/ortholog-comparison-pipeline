#!/bin/bash

# Example run:
# $blastdir/convert_blast_to_fasta.sh -b $outfile -r -g -c "genome" -o $rundir/sample_fasta.fasta

print_help() {
    echo "Usage: $0 -b <blast> -r -e -c <'genome'|'locus'> -o <outfile> -h"
    echo "        Output: FASTA file (with headers >ID|title|evalue)"
    echo ""
    echo "  -b    Required. Path to BLAST file in which the first five columns are sgi sseqid sseq evalue stitle, and optionally the sixth column is locus tag."
    echo "  -r    Optional. Include this flag to remove the first row, assumed to be the header (default: false)."
    echo "  -e    Optional. Extract if the sequence ID needs to be extracted from within | symbols; do not use otherwise."
    echo "  -g    Optional. Remove gap characters ("-") (default: false)."
    echo "  -c    Optional. Concatenate IDs for sequence name. Default behavior: protein ID used as sequence name."
    echo "                    -c 'genome': Use if you want to concatenate 'genomeID-sequenceID'."
    echo "                    -c 'locus': Use if you want to concatenate 'genomeID-sequenceID-locusTag'. Locus tag is assumed to be in column 6."
    echo "                        If you need to rearrange the input file to satisfy this requirement, use awk."
    echo "                        Example where locus tag is in column 18: awk -F'\t' -v OFS='\t' '{ tmp=$6; $6=$18; $18=tmp; print }' $infile > $outfile"
    echo "  -o    Optional. output filename (by default, ${blast%.*}.fasta)."
    echo "  -h    Show help message and exit"
}

# set default values for optional arguments
header=false
extract=false
gapRemoval=false
concat="no"

# Parse options with getopts (not GNU getopts, as that causes issues on Mac systems)
while getopts "b:regc:o:h" opt; do
  case $opt in
    b) blast="$OPTARG" ;;
    r) header="true" ;;
    e) extract="true" ;;
    g) gapRemoval="true" ;;
    c) concat="$OPTARG" ;;
    o) outfile="$OPTARG" ;;
    h) print_help; exit 0 ;;
    *) print_help; exit 1 ;; #echo "Invalid option"; exit 1 ;;
  esac
done

if [[ -z "$blast" ]]; then
  echo "Error: BLAST file not provided. Use -b or --blast to specify input." >&2
  exit 1
fi

# Assign default outfile only if not set
output="${outfile:-${blast%.*}.fasta}"
mkdir -p $(dirname $output) # make the outdir if it doesn't exist already

> "$output"

echo "Writing results to ${output}..."

IFS=$'\n'

# If there's a header in the file, skip it.
if [[ "$header" = true ]]; then
  tail_input=$(tail -n +2 "$blast")
else
  tail_input=$(cat "$blast")
fi

# Use awk to process the input with proper empty field handling
awk -F'\t' -v extract="$extract" -v gapRemoval="$gapRemoval" -v concat="$concat" '
{
    genome = $1
    id = $2
    sequence = $3
    evalue = $4
    title = $5
    locus = $6

    if (extract == "true") {
        sub(/^[^|]*\|/, "", id)  # Remove everything before first |
        sub(/\|.*$/, "", id)      # Remove everything after first |
    }
    
    # If requested, remove gap characters from sequences (col 3)
    if (gapRemoval == "true") {
        gsub("-", "", sequence)
    }

    # Generate the appropriate fasta format based on concat mode
    if (concat == "genome") {
        printf ">%s-%s %s|%s\n%s\n", genome, id, title, evalue, sequence
    } else if (concat == "locus") {
        printf ">%s-%s-%s %s|%s\n%s\n", genome, id, locus, title, evalue, sequence
    } else {
        printf ">%s %s|%s\n%s\n", id, title, evalue, sequence
    }
}
' <<< "$tail_input" >> "$output"



echo "Finished converting BLAST to FASTA: ${output}"