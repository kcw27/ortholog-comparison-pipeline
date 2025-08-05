#!/bin/bash

### CLIs:
# $1: alignment FASTA to convert to PHYLIP. Assumes that the only thing in the headers are sequence names, which should be the case for ClustalW outputs anyway.
# $2: sequence type, i.e. alphabet for seqmagick {dna,dna-ambiguous,rna,rna-ambiguous,protein}
# $3: outdir

### Outputs: to the outdir ($2), saves the following:
## name_map.tsv, a table with two columns. The first column features short sequence names that won't be truncated by PHYLIP format,
# and the second column features long sequence names (i.e. the original sequence names from the FASTA $1).
## renamed.fasta, a copy of $1 but with short sequence names substituted in for the long sequence names.
## PHYLIP version of $1, named aligned.phy.

# Example run:
# bash fasta_to_phylip.sh "${HOME}/data/POG002101_aligned.fasta" "protein" "${HOME}/data/POG002101_seqs/"
# bash fasta_to_phylip.sh "${HOME}/data/PA3565_orthologs_65_66_67_aligned.fasta" "protein" "${HOME}/data/PA3565_66_67_seqs/"

# Install seqmagick if it hasn't been already
# pip install seqmagick

# Set up variables
infile="$1"
seqtype="$2"
outdir="${3%/}" # remove trailing slash if present
mapfile="${outdir}/name_map.tsv"
renamedfile="${outdir}/renamed.fasta" # sequence names renamed from long to short
#outfile="${outdir}/$(basename $infile | cut -d . -f 1).phy" # names the output PHYLIP after the input FASTA
outfile="${outdir}/aligned.phy" # decided to standardize the filename for ease of use

# Make the outdir if needed
mkdir -p $outdir

## Make the name map table
grep ">" "$infile" | sed 's/>//' | awk '{printf "seq%04d\t%s\n", NR, $0}' > "$mapfile"

## Save a copy of the input fasta $1 with long sequence names switched to short sequence names
# Create an associative array to map long names to short names
declare -A name_map
while IFS=$'\t' read -r short long; do
    name_map["$long"]="$short"
done < "$mapfile"

# Rewrite the FASTA with substituted names
awk -v OFS="" -v mapfile="$mapfile" '
BEGIN {
    # Load name_map into memory
    while ((getline < mapfile) > 0) {
        short = $1
        long = $2
        map[long] = short
    }
}

{
    if ($0 ~ /^>/) { # find FASTA sequence names, which follow ">"
        split($0, parts, " ")
        gsub(/^>/, "", parts[1])
        if (parts[1] in map) {
            print ">" map[parts[1]]
        } else {
            print $0  # fallback if mapping not found
        }
    } else {
        print $0
    }
}' "$infile" > "$renamedfile"

## Convert renamed.fasta to strict PHYLIP with seqmagick convert
seqmagick convert --alphabet "$seqtype" "$renamedfile" "$outfile"

echo "Name map, FASTA with renamed sequences, and output PHYLIP file have been written to ${outdir}!"