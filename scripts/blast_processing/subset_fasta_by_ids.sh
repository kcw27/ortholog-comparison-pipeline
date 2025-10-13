#!/bin/bash

# This script was AI-generated. I haven't tested it yet.
# Potential use case: in advance of alignment_and_tree_wrapper.sh, subset dataset to sequences of interest
# to produce a more legible phylogenetic tree

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <fasta_file> <id_list_file> <output_file>"
    exit 1
fi

# Input files
FASTA_FILE="$1"
ID_LIST_FILE="$2"
OUTPUT_FILE="$3"

# Validate that the input files exist
if [ ! -f "$FASTA_FILE" ]; then
    echo "Error: FASTA file '$FASTA_FILE' not found."
    exit 1
fi

if [ ! -f "$ID_LIST_FILE" ]; then
    echo "Error: ID list file '$ID_LIST_FILE' not found."
    exit 1
fi

# Create a temporary file to store IDs as a regex pattern
ID_PATTERN_FILE=$(mktemp)
awk '{print "^>" $0 "$"}' "$ID_LIST_FILE" > "$ID_PATTERN_FILE"

# Extract sequences from the FASTA file
awk -v id_file="$ID_PATTERN_FILE" '
BEGIN {
    # Load IDs into an array for quick lookup
    while ((getline id < id_file) > 0) {
        ids[id] = 1
    }
    close(id_file)
}
{
    if ($0 ~ /^>/) {
        # Check if the current header matches any ID
        match_found = ($0 in ids)
    }
    if (match_found) {
        print $0
    }
}
' "$FASTA_FILE" > "$OUTPUT_FILE"

# Clean up temporary file
rm -f "$ID_PATTERN_FILE"

echo "Subset FASTA file created: $OUTPUT_FILE"
