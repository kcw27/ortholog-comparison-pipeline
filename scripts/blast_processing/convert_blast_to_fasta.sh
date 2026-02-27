#!/bin/bash

# Example run:
# $blastdir/convert_blast_to_fasta.sh -b $outfile -r -c "genome" -o $rundir/sample_fasta.fasta

print_help() {
    echo "Usage: $0 -b <blast> -r -e -c <'genome'|'locus'> -o <outfile> -h"
    echo "        Output: FASTA file (with headers >ID|title|evalue)"
    echo ""
    echo "  -b    Required. Path to BLAST file in which the first five columns are sgi sseqid sseq evalue stitle, and optionally the sixth column is locus tag."
    echo "  -r    Optional. Include this flag to remove the first row, assumed to be the header (default: false)."
    echo "  -e    Optional. Extract if the sequence ID needs to be extracted from within | symbols; do not use otherwise."
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
concat="no"

# Parse options with getopts (not GNU getopts, as that causes issues on Mac systems)
while getopts "b:rec:o:h" opt; do
  case $opt in
    b) blast="$OPTARG" ;;
    r) header="true" ;;
    e) extract="true" ;;
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

## define the write_fasta() function differently depending on concatenation mode
#if [[ "$concat" = "genome" ]]; then
#  echo "Concatenating genome to subject ID"
#  write_fasta() {
#    printf ">%s-%s %s|%s\n%s\n" "$1" "$2" "$3" "$4" "$5" >> "$output"
#  }
#elif [[ "$concat" = "locus" ]]; then
#  echo "Concatenating genome and locus tag to subject ID"
#  write_fasta() {
#    printf ">%s-%s-%s %s|%s\n%s\n" "$1" "$2" "$6" "$3" "$4" "$5" >> "$output"
#  }
#else
#  echo "Not concatenating genome to subject ID"
#  write_fasta() {
#    printf ">%s %s|%s\n%s\n" "$2" "$3" "$4" "$5" >> "$output"
#  }
#fi

#while read -r next; do
#  # Extract columnsin one go
#  #read id title evalue sequence < <(awk -F '\t' '{print $2, $5, $4, $3}' <<< "$next")
#  #IFS=$'\t' read -r col1 col2 col3 col4 col5 <<< "$next"
#  IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col_rest <<< "$next"
#  # there may be more than 5 columns, so if there are, they get dumped into col_rest. (If there aren't, col_rest will just be empty.)
#  
#  genome="$col1"
#  id="$col2"
#  title="$col5"
#  evalue="$col4"
#  sequence="$col3"
#  locus="$col6" # no harm if there's no 6th column in the file
#
#  if [[ "$extract" = true ]]; then
#    id="${id#*|}"
#    id="${id%%|*}"
#  fi
#
#  
#  write_fasta "$genome" "$id" "$title" "$evalue" "$sequence" "$locus"
#
#done <<< "$tail_input"

# Use awk to process the input with proper empty field handling
awk -F'\t' -v extract="$extract" -v concat="$concat" '
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

# old example runs:
# on synteny-filtered output (-hd is required because it has a header):
# bash $blastscripts/convert_blast_to_fasta.sh -b $datadir/PA3565_67_topHitsPerGenome.blast -hd

# on blast output that hasn't been through metadata processing (i.e. blast2gen.py or intersect_blast_and_synteny_with_metadata.R) and therefore has no headers (-e)
# bash $blastscripts/get_blast_top_hits.sh $datadir/PA3565_with_orgs_long.blast $datadir/PA3565_top_byOrganism.blast

# on blast output that has been through a filtering script that left it without a header:
# /home/kcw2/ortholog-comparison-pipeline/test_data/PA3565_top_byGenomeID.fasta

# on blast2gen.py output which hasn't gone through further filtering (-hd and -e)
# bash $blastscripts/convert_blast_to_fasta.sh -b $datadir/PA3565_with_orgs_long_annotatedByBlast2gen.blast -e -hd

# with specified outname:
# bash "/home/kcw2/ortholog-comparison-pipeline/scripts/blast_processing/convert_blast_to_fasta.sh" -b "/home/kcw2/data/testing/PA3565_nr_small_foo.txt" -e -g -o "/home/kcw2/data/testing/PA3565_nr_small_foo_bar.fasta"

## Parse flags
#while [[ $# -gt 0 ]]; do
#  case "$1" in
#    -b|--blast)
#      blast="$2"
#      shift 2
#      ;;
#    -hd|--header)
#      header=true
#      shift
#      ;;
#    -e|--extract)
#      extract=true
#      shift
#      ;;
#    -g|--concat_genome)
#      concat_genome=true
#      shift
#      ;;
#    -o|--outfile)
#      outfile="$2"
#      shift 2
#      ;;
#    -h|--help)
#      echo "I will add a help string later"
#      exit 1
#      ;;
#    *)
#      echo "Unknown option: $1" >&2
#      exit 1
#      ;;
#  esac
#done