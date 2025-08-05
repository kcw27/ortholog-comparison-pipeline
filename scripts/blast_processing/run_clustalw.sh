#!/bin/bash

# CLIs:
# $1: multi-FASTA to align with ClustalW
# $2: alignment type (either "PROTEIN" or "DNA")
# $3: output format (e.g. "FASTA" or "PHYLIP")
# $4: output filename

# Output: preprocesses the input multi-FASTA ($1) to keep only non-empty sequences with unique names,
# then performs a ClustalW alignment using $2 and $3 as parameters, and saves the output to the specified output filename $4.

# Example run:
# bash run_clustalw.sh "${HOME}/data/POG002101.fasta" "PROTEIN" "FASTA" "${HOME}/data/POG002101_aligned.fasta"
# bash run_clustalw.sh "${HOME}/data/PA3565_orthologs_65_66_67_top_evalueThreshold_1e-50.fasta" "PROTEIN" "FASTA" "${HOME}/data/PA3565_orthologs_65_66_67_aligned.fasta"

tempdir="$(dirname $4)/temp_clustalw"
mkdir -p $tempdir
echo "Saving temp files to ${tempdir}. Log will also be written there."
log="${tempdir}/run_clustalw_benchmarking_log.txt"
> $log
# echo "Benchmarking: $(wc -l $1 | cut -d ' ' -f 1) lines in input file" # cut so it doesn't display the filename too
#echo "Benchmarking: number of lines in input file is $(wc -l $1)" # actually, leave the filename on
#echo "Please refer to the input file to determine the original number of sequences."

# get sequences with unique names in the fasta
#sed -e '/^>/s/$/@/' -e 's/^>/#/' "$1" | tr -d '\n' | tr "#" "\n" | tr "@" "\t" > "${tempdir}/temp.tsv"
# This line creates a tab-separated table with two columns: sequence name and sequence
# @ is placed at the end of header lines, and the leading > is replaced with # to mark the start of header lines
# (because ">" is risky when dealing with text processing)
# Newlines are removed, then restored wherever "#" occurs. "@" are replaced with tabs.

# trying a version that doesn't replace > with #, since # appears in some FASTA headers
awk '/^>/ {if (seq) print hdr "\t" seq; hdr=$0; seq=""} !/^>/ {seq=seq $0} END {print hdr "\t" seq}' "$1" > "${tempdir}/temp.tsv"
echo "Benchmarking: number of sequences in original file was $(wc -l ${tempdir}/temp.tsv | cut -d ' ' -f 1)" >> $log 
#echo "Benchmarking: number of sequences in original file was $(expr $(wc -l ${tempdir}/temp.tsv | cut -d ' ' -f 1) / 2)" >> $log # incorrect; don't need to divide by 2

# sort -u -t $'\t' -k 2,2 "${tempdir}/temp.fasta" > "${tempdir}/temp_sorted.fasta"
# BUG: this originally only kept unique sequences (col 2), but we actually wanted to filter out duplicate sequence names (col 1) instead

sort -u -t $'\t' -k 1,1 "${tempdir}/temp.tsv" > "${tempdir}/temp_unique.tsv"
# Sort by sequence name (column 1) so that the -u flag will keep only records with unique sequence names
echo "Benchmarking: number of sequences with unique names is $(wc -l ${tempdir}/temp_unique.tsv)" >> $log 
#echo "Benchmarking: number of sequences with unique names is $(expr $(wc -l ${tempdir}/temp_unique.tsv | cut -d ' ' -f 1) / 2)" >> $log # incorrect; don't need to divide by 2

awk -F'\t' '{print ">"$1"\n"$2}' "${tempdir}/temp_unique.tsv" > "${tempdir}/temp_unique.fasta"
# Converts the table back to FASTA: a line for header, then a line for sequence

# get non-empty sequences in the fasta; also removes malformed entries
awk 'BEGIN {RS=">"; ORS=""} length($2) > 0 {print ">"$0}' "${tempdir}/temp_unique.fasta" > "${tempdir}/temp_unique_nonempty.fasta"
# Records must start with > and contain a sequence
echo "Benchmarking: the number of sequences in the FASTA used for the alignment, ${tempdir}/temp_unique_nonempty.fasta, is $(expr $(wc -l ${tempdir}/temp_unique_nonempty.fasta | cut -d ' ' -f 1) / 2)" >> $log
# Divided the linecount by 2 to get the number of sequences because the headers and sequences are on separate lines in the FASTA

clustalw -INFILE="${tempdir}/temp_unique_nonempty.fasta" -TYPE="$2" -OUTPUT="$3" -OUTFILE="$4"

# remove temp files
# I decided to pipe "yes" to rm instead of using -f
# Depending on the system, rm might be aliased to rm -i, so this ensures that the files are removed
# yes | rm "${HOME}/temp.fasta" "${HOME}/temp_sorted.fasta" "${HOME}/temp_unique.fasta" "${HOME}/temp_unique_sorted.fasta"
