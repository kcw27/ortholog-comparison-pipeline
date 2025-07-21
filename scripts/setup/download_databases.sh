#!/bin/bash

# https://github.com/kblin/ncbi-genome-download

#cd /home/share
#ncbi-genome-download --section genbank bacteria --parallel 12 

# Modifying this code to only download complete genomes for organisms from the BLAST
# blast_PA3565_nr.sh -> convert_organisms_to_accessions.sh -> sort_accessions.sh
# The two outputs of sort_accessions.sh are the inputs to the lines below.

# CLIs:
# $1: file containing a list of GenBank accessions
# $2: file containing a list of RefSeq accessions
# $3: directory in which to download the databases

# I should later set this up so that the user can use flags rather than relying on the order of CLIs
# Also set up an optional argument for the section to search, with bacteria being the default
gb=$1
ref=$2
dbdir=$3

mkdir -p $dbdir
cd $dbdir # so that the databases are placed in the correct location
#cd "${HOME}/data" # since these aren't complete database downloads, I won't put them in the shared directory

# the commands below assume the user wants to download bacterial genomes
# yet another thing that could be refactored into flags
ncbi-genome-download --section genbank --assembly-accessions "${gb}" bacteria --parallel 12 --progress-bar &
#ncbi-genome-download --section genbank --assembly-accessions "${gb}" bacteria --parallel 12 --progress-bar --dry-run # for testing
echo "Downloading genbank with job id ${!}"

ncbi-genome-download --section refseq --assembly-accessions "${ref}" bacteria --parallel 12 --progress-bar &
#ncbi-genome-download --section refseq --assembly-accessions "${ref}" bacteria --parallel 12 --progress-bar --dry-run # for testing
echo "Downloading refseq with job id ${!}"

echo "Moving genomes to ${dbdir}..."
mv "${dbdir}/genbank/bacteria/*" "${dbdir}"
mv "${dbdir}/refseq/bacteria/*" "${dbdir}"

echo "Removing now-empty directories..."
rm -rf "${dbdir}/genbank"
rm -rf "${dbdir}/refseq"