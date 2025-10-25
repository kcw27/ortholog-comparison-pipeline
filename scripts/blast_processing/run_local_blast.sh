#!/bin/bash

# CLIs:
# 1: path to BLAST db
# 2: database name (blastp -db)
# 3: query FASTA
# 4: output name

# Output: local blast of $3 against $2, saved to $4

# Example run:
# bash /home/kcw2/ortholog-comparison-pipeline/scripts/blast_processing/run_local_blast.sh "/home/kcw2/temp_db/" temp_db /home/kcw2/data/PAO1_PA3565.fasta /home/kcw2/data/blast_test/pa3565_test.blast &
# bash /home/kcw2/ortholog-comparison-pipeline/scripts/blast_processing/run_local_blast.sh /home/share/pa/pa_db pa /home/kcw2/data/PAO1_PA3565.fasta $HOME/PA3565_orthologs_test.blast & # for other people to test
# bash run_local_blast.sh /home/share/nr nr /home/kcw2/data/PAO1_PA3565.fasta /home/kcw2/data/PA3565_orthologs.txt &

cd $1 # so that it can find the blast database
procs_to_use=$(( $(nproc) / 4 ))
mkdir -p $(dirname $4) # make the outdir if it doesn't exist already
scriptsdir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # get location of current script; this is where the other scripts are found too

tempBlast=$(mktemp)
echo "Running BLAST... Temp file at $tempBlast"
blastp -db "$2" -query "$3" -max_target_seqs 500000 -show_gis -num_threads $procs_to_use -outfmt "6 sallgi sallseqid sseq evalue salltitles" -out $tempBlast

# add organisms column
tempBlastOrgs=$(mktemp)
echo "Adding organisms... Temp file at $tempBlastOrgs"
bash "${scriptsdir}/add_organism_column.sh" "$tempBlast" "$tempBlastOrgs"

# save BLAST output as long format
echo "Expanding BLAST output file to long format..."
bash "${scriptsdir}/expand_blast_output.sh" -b "$tempBlastOrgs" -o "$4"