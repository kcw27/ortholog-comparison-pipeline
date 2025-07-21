#!/bin/bash

# Command line inputs:
# $1: name of file containing a mix of RefSeq (GCF) and GenBank (GCA) assembly accession numbers,
# which is the output of convert_organisms_to_accessions.sh.
# $2: outdir to which the two output files will be written

# Outputs: one file containing RefSeq accessions and one file containing genbank accessions.
# Files will be placed in $2, named after the basename from $1 with _refseq or _genbank added to their names.

# example run:
# bash sort_accessions.sh ~/data/PA3565_accessions_small.txt ~/data

name_root=$(basename "${1}" | cut -d "." -f 1) 
# basename to strip the filepath
# assume the name of the file is the part that cones before the first period

cat $1 | grep "GCF_" > "${2%/}/${name_root}_refseq.txt" # if there's a "/" at the end of $2, strip it so there's no redundant slash; also add the suffix to the filename
cat $1 | grep "GCA_" > "${2%/}/${name_root}_genbank.txt"

echo "Wrote outputs to ${2%/}/${name_root}_refseq.txt and ${2%/}/${name_root}_genbank.txt."