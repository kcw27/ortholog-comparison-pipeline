#!/bin/bash

# To find sequences that appear in both the blast dataset and synteny hit dataset without relying on consistent protein IDs,
# perform a pairwise BLAST of the latter against the former, and filter the output to hits with >=99 pident.
# We have to use the synteny hit dataset as the query and the blast hit dataset as the subject, because we're taking a subset of the former
# (i.e. appears in the synteny context) based on which of its sequences resembles sequences that appear in the latter (i.e. actually looks like the query sequence).
# The question is "does this sequence from the synteny hits dataset look anything like any sequence from the dataset I got from blasting my protein of interest against a large database?"

# Output:
# To the specified output filename, save a version of the full output from the pairwise BLAST (i.e. intersection)
# Also, saves a fasta which has been filtered to exclude sequences with <$pident.

# Example run:
# TODO
# for testing:
# bash "/home/kcw2/ortholog-comparison-pipeline/scripts/synteny_search/intersect_blast_and_synteny_pairwiseBlastApproach.sh" -b /home/kcw2/data/blast_outputs/pseudomonas_aeruginosa_PA3565_orgs.blast -s /home/kcw2/data/results_65_67/synteny_summary.tsv -o /home/kcw2/data/example_run/foo.blast



print_help() {
    echo "Usage: $0 -b <blast> -s <synteny> -p <pident> -c<'genome'|> -o <outname> -h"
    echo "        Outputs: to $outname, writes the output of the pairwise BLAST of $synteny as query against $blast as subject."
    echo "                 Rationale for this direction of BLAST is explained in comments at the top of this script."
    echo "                 Also writes a fasta with $pident cutoff- look for a file with a name similar to outname."
    echo "  -b    Required. Path to BLAST file in which the first five columns are sgi sseqid sseq evalue stitle."
    echo "  -s    Required. Path to summary file from synteny search via synteny_wrapper.sh (which calls find_synteny_hits.sh)."
    echo "  -p    Optional. Pident threshold for filtering (default = 99)."
    echo "  -c    Optional. Concatenation mode for convert_blast_to_fasta.sh (default = 'no')."
    echo "  -o    Required. Name of output file."
    echo "  -h    Show help message and exit"
}

# set default values for optional arguments
pident="99"
concat="no"

# Parse options
while getopts "b:s:p:c:o:h" opt; do
  case $opt in
    b) blast="$OPTARG" ;;
    s) synteny="$OPTARG" ;;
    p) pident="$OPTARG" ;;
    c) concat="$OPTARG" ;;
    o) outname="$OPTARG" ;;
    h) print_help; exit 0 ;;
    *) print_help; exit 1 ;; #echo "Invalid option"; exit 1 ;;
  esac
done

currdir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # get location of current script
scriptsdir="${currdir}/../blast_processing" # this is where I put the scripts that this script calls

mkdir -p $(dirname $outname) # make the outdir if it doesn't exist already

### Convert to FASTA
echo "Converting input files to FASTA after preprocessing..."

echo "Processing $blast"
tempBlast=$(mktemp)
tempBlastFasta=$(mktemp)
echo "Saving tempBlastFasta to ${tempBlastFasta}..."

# remove gap characters "-" because command line BLAST doesn't recognize them.
# only write a line if the sequence hasn't been encountered before; this ensures the FASTA only contains unique seqs
awk 'BEGIN {OFS="\t"} {
  gsub("-", "", $3);
  if (!seen[$3]++) print
}' "$blast" > "$tempBlast"
# very important to set the delimiter as tab; awk defaults to space delimiter
bash "${scriptsdir}/convert_blast_to_fasta.sh" -b "$tempBlast" -o "$tempBlastFasta" &
pid1=$!


echo "Processing $synteny"
tempSynteny=$(mktemp)
tempSyntenyFasta=$(mktemp)
echo "Saving tempSyntenyFasta to ${tempSyntenyFasta}..."

# need to reorder columns so it's compatible with convert_blast_to_fasta.sh
# or rather, I think the awk command just selects the columns specified
awk -F'\t' 'BEGIN {OFS="\t"} {print $1, $7, $8, "na", $5, $6}' "$synteny" > "$tempSynteny"
bash "${scriptsdir}/convert_blast_to_fasta.sh" -b "$tempSynteny" -o "$tempSyntenyFasta" -c "$concat" -r & # -r to remove header
pid2=$!


wait $pid1 $pid2
yes | rm $tempBlast
yes | rm $tempSynteny

### Make blast db from $blast file
# Though it is possible to blast two fastas against each other, making a db allows for use of the -num_threads argument.
# source activate blast_env # for the final pipeline, have the user activate their own blast environment prior to running this script

# assume it's a protein sequence database
temp_dir=$(mktemp -d)
echo "Temporary directory for blast db created at: $temp_dir"
mkdir -p "${temp_dir}/temp_db"
cd $temp_dir # need to be one directory above the blast db being made
makeblastdb -in "$tempBlastFasta" -input_type fasta -dbtype prot -title temp_db -out "${temp_dir}/temp_db" &
pid3=$!
wait $pid3
yes | rm -f $tempBlastFasta

### Perform the pairwise BLAST
# conda deactivate # currently, local blast is in the base conda environment away from blast_env
# edit: local blast IS in blast_env???
echo "Saving pairwise blast outputs..."

tempBlastOutputs=$(mktemp) # going to be using awk to add a header later, so for now I'm saving it to a temp file
procs_to_use=$(( $(nproc) / 4 ))
max=10
blastp -query "$tempSyntenyFasta" -db "${temp_dir}/temp_db" -max_target_seqs "$max" -num_threads "$procs_to_use" -outfmt "6 qseqid qseq sallseqid sseq pident evalue bitscore" -out "$tempBlastOutputs" &
pid4=$!

wait $pid4

### Delete the blast db that was created
#cleanup-blastdb-volumes.py -db temp_db -dbtype prot
rm -rf "$temp_dir"

### Filter the output file from the pairwise blast to only include lines with >=$pident
tempFilteredIDs=$(mktemp)
awk -v p="$pident" '$5 >= p' "$tempBlastOutputs" | cut -f 1 | sort | uniq > "$tempFilteredIDs" # IDs were in column 1 of the blast that was just done, pident was in column 5

### Filter the synteny fasta to only include the lines from the filtered IDs list
# This filtered FASTA is the final output of this script. You can still left join the synteny metadata to the seq IDs in that filtered FASTA; it's a subset.
fastaname="${outname%.*}_pident${pident}.fasta"
seqtk subseq "$tempSyntenyFasta" "$tempFilteredIDs" > "$fastaname"
echo "Saved filtered FASTA to $fastaname"

# Add colnames to the file saved to $outname
awk 'BEGIN { OFS="\t"; print "qseqid", "qseq", "sallseqid", "sseq", "pident", "evalue", "bitscore" } { print }' "$tempBlastOutputs" > "$outname"

yes | rm -f $tempSyntenyFasta 