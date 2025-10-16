#!/bin/bash

##### TODO:
### 0: Important!!!: update tempBlastFasta and tempSyntenyFasta to use mktemp; I don't want to overwrite these temp files
# 1: reverse the direction of the blast: the blast DB should be built from $blast, and the query should be syntenyFasta
# 2: --num-threads should be calculated based on num cores
# 3: make a text file of list of IDs from the output of the pairwise BLAST for which pident ($5) >= 80
### e.g., awk '$5 >= 80' /home/kcw2/data/blast_outputs/pseudomonas_aeruginosa_PA3565_67_intersectBlastSynteny_syntenyQuery_blastOutputSubject.blast | cut -f 1 | sort | uniq > /home/kcw2/data/blast_outputs/ids_from_pairwiseBlast_pidentThreshold80.txt 
# 4: filter the syntenyFasta using seqtk subseq and that text file (#3 in TODO list)
### So, don't delete syntenyFasta until the very end.
# That filtered FASTA is the final output of this script. You can still left join the synteny metadata to the seq IDs in that filtered FASTA; it's a subset.


# To find sequences that appear in both the blast dataset and synteny hit dataset without relying on consistent protein IDs,
# perform a pairwise BLAST of the former against the latter, and filter the output to hits with >=99 pident.
# We have to use the blast dataset as the query and the synteny hit dataset as the subject, because the blast dataset features
# aligned sequences that are as a result shorter. 100 pident is possible if a shorter sequence is aligned against a longer sequence
# (and would indicate a perfect match between the blast dataset sequence and a sequence in the synteny hit dataset),
# but 100 pident is not achievable if a longer sequence (i.e. from the synteny hit dataset) is aligned against a shorter sequence.

# Input:
# $1: BLAST hits file (-outfmt "6 sallgi sallseqid sseq evalue salltitles")
# $2: synteny hits file from synteny pipeline
# $3: output filename

# Output:
# To the specified output filename, save a version of the output from the pairwise BLAST (i.e. intersection) which has been filtered to exclude <99 pident.
# Also save the full pairwise BLAST output to ""${outname%.*}"_full.blast"

# Example run:
# TODO
# for testing:
# bash "/home/kcw2/ortholog-comparison-pipeline/scripts/synteny_search/intersect_blast_and_synteny_pairwiseBlastApproach.sh" /home/kcw2/data/blast_outputs/pseudomonas_aeruginosa_PA3565_orgs.blast /home/kcw2/data/results_65_67/synteny_summary.tsv foo

blast=$1
synteny=$2
outname=$3

currdir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # get location of current script
scriptsdir="${currdir}/../blast_processing" # this is where I put the scripts that this script calls

### Convert to FASTA
echo "Converting input files to FASTA after preprocessing..."

echo "Processing $blast"
tempBlast=$(mktemp)

# eventually replace this with $(mktemp) and delete it afterward, but for now I want to check the FASTAs
tempBlastFasta="/home/kcw2/data/temp/blastFasta.fasta"

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

# eventually replace this with $(mktemp) and delete it afterward, but for now I want to check the FASTAs
tempSyntenyFasta="/home/kcw2/data/temp/syntenyFasta.fasta"

# need to reorder columns so it's compatible with convert_blast_to_fasta.sh
awk -F'\t' 'BEGIN {OFS="\t"} {print $1, $7, $8, "na", $5}' "$synteny" > "$tempSynteny"
bash "${scriptsdir}/convert_blast_to_fasta.sh" -b "$tempSynteny" -o "$tempSyntenyFasta" -g &
pid2=$!


wait $pid1 $pid2
yes | rm $tempBlast
yes | rm $tempSynteny

### Make blast db from synteny file
# Though it is possible to blast two fastas against each other, making a db allows for use of the -num_threads argument.
source activate blast_env # TODO: need to change this for the final pipeline

# assume it's a protein sequence database
mkdir -p ~/temp_db # sorry, I'll clean my code up later
cd ~

makeblastdb -in "$tempSyntenyFasta" -input_type fasta -dbtype prot -title temp_db -out /home/kcw2/temp_db/temp_db &
pid3=$!
wait $pid3
yes | rm -f $tempSyntenyFasta

### Run local blast
conda deactivate # currently, local blast is in the base conda environment away from blast_env


#cd ~/temp_db # so that local blast can find the database

# NOTE: for now, I'm not using $outname. I'll incorporate it once I decide on a value of -max_target_seqs.
# I'll also calculate the number of processors to use for -num_threads instead of hard-coding it.

max=10
blastp -query "$tempBlastFasta" -db /home/kcw2/temp_db/temp_db -max_target_seqs "$max" -num_threads 9 -outfmt "6 qseqid qseq sallseqid sseq pident evalue bitscore" -out "$outname" &
pid4=$!

wait $pid4
yes | rm -f $tempBlastFasta
#cd .. # out of temp_db directory

### Delete the blast db that was created
#cleanup-blastdb-volumes.py -db temp_db -dbtype prot
rm /home/kcw2/temp_db/temp_db.*

### TODO: filter the single output file to >=99 pident or an otherwise appropriate threshold.