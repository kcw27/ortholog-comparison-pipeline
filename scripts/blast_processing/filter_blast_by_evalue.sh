#!/bin/bash

# Inputs:
# -t (--threshold); optional with default value 1e-30
# -f (--files): list of BLAST output files (evalue in column 4, e.g. -outfmt "6 sallgi sallseqid sseq evalue salltitles") to filter

# Outputs: filters the input files to only include evalues (column 4) which are smaller than the provided threshold.
# Saves them to output files named like so: "${file%.*}_evalueThreshold_${THRESHOLD}.txt"

# Example run:
# bash filter_blast_by_evalue.sh -t "1e-50" -f ${HOME}/data/PA3565_orthologs_top.txt ${HOME}/data/PA3565_orthologs_65_67_top.txt ${HOME}/data/PA3565_orthologs_65_66_67_top.txt

# Default threshold:
THRESHOLD="1e-30"

# Parse options
while [[ $# -gt 0 ]]; do
    case "$1" in
        -t|--threshold)
            THRESHOLD="$2"
            shift 2
            ;;
        -f|--files)
            FILES=("${@:2}")
            break
            ;;
        *)
            echo "Invalid option: $1"
            exit 1
            ;;
    esac
done

# Check if file arguments were provided
if [[ ${#FILES[@]} -eq 0 ]]; then
    echo "Usage: $0 [-t <threshold>] -f <file1> <file2> ..."
    exit 1
fi

# Process each file
for file in "${FILES[@]}"; do
    output="${file%.*}_evalueThreshold_${THRESHOLD}.txt"
    awk -v threshold="$THRESHOLD" '$4 < threshold' "$file" > "$output" # -v is the safest way to pass a variable to awk
    echo "Processed: $file -> $output"
done
