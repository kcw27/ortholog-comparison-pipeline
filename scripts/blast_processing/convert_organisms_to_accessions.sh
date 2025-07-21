#!/bin/bash

source activate blast_env # activate the conda environment in which BLAST is installed
# The line above may fail if you're not in the base conda environment

# Command line inputs:
# $1: name of input file, which is a command line BLAST output with -outfmt "6 sgi sseqid sseq evalue stitle"
# $2: name of output file, which is a list of assembly accessions (GCAs and GCFs) corresponding to the organism names in $1
# Note that because genomes are not reliably labeled with strains, this will pull assembly accessions for each genus + species combo
# in $1, without regard to strain.

# Example run:
# bash convert_organisms_to_accessions.sh ~/data/PA3565_orthologs.txt ~/data/PA3565_accessions.txt

# First, get the list of unique organisms.
# This is done by extracting any text found between square brackets [],
# then taking the first two words in the organism name (i.e. genus and species) to ignore strain,
# then sort | uniq to remove duplicates,
# then writing the list to a temporary text file.
echo "Will retrieve assembly numbers for each organism name in the BLAST input."

#grep -oP '(?<=\[).+(?=\])' $1 | cut -d " " -f 1,2 | sort | uniq > ~/temp.txt

# refactor to deal with multiple titles per record (salltitles is planned to be part of the new BLAST outfmt expected by the pipeline)
cat $1 | cut -f 5 | tr ';' '\n' | grep -oP '(?<=\[).+(?=\])' | sort | uniq > ~/temp.txt
echo "Unique organism list obtained."

# create a second temp file to which assembly accessions will be written (will contain redundant records) 
> ~/temp2.txt # starts as empty

echo "Querying NCBI..."
# iterate through each genus + species in the temp file, querying NCBI to get the corresponding AssemblyAccession

IFS=$'\n' # temp.txt is newline-separated.
# Need to set up the loop this way, otherwise it only reads the first line
for next in $(cat ~/temp.txt); do
  echo "Processing: ${next}"
  #esearch -db assembly -query "${next}" -retmax 10 | \
  esearch -db assembly -query "${next} AND complete genome[filter]" | \
  efetch -format docsum | \
  xtract -pattern DocumentSummary -element AssemblyAccession >> ~/temp2.txt
  
  sleep 2 # avoid overloading NCBI
done

echo "NCBI queries completed. Processing outputs..."

# retmax is not working as intended, but in theory, it caps the number of unique identifiers returned to 10.
# If you would like to limit the number of hits per query, just insert "head -n $upper_limit" right before it writes to ~/temp2.txt.

# Important: use "${line}" rather than $line so that it interprets the entire string as one argument rather than two separate arguments

# AssemblyAccession comes in two columns: the assembly number and the organism name.
# I will take only the assembly number, but if you'd like to know the organism name, you can remove the "cut" step.
# Write the list of unique identifiers to the output file.
cat ~/temp2.txt | cut -f 1 | sort | uniq > $2

# remove temp files
rm ~/temp.txt 
#rm ~/temp2.txt # Let's keep temp2.txt for now.

echo "Done retrieving assembly numbers for each organism name."