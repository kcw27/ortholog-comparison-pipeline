#!/bin/bash

# CLIs:
# 1: path to BLAST db
# 2: database name (blastp -db)
# 3: query FASTA
# 4: output name

# Output: local blast of $3 against $2, saved to $4

# Example run:
# bash /home/kcw2/ortholog-comparison-pipeline/scripts/blast_processing/run_local_blast.sh /home/share/pa/pa_db /home/kcw2/data/PAO1_PA3565.fasta $HOME/PA3565_orthologs_test.txt & # for other people to test
# bash run_local_blast.sh /home/share/nr /home/kcw2/data/PAO1_PA3565.fasta /home/kcw2/data/PA3565_orthologs.txt &

cd $1 # so that it can find the blast database
blastp -db "$2" -query "$3" -max_target_seqs 500000 -show_gis -num_threads 12 -outfmt "6 sallgi sallseqid sseq evalue salltitles" -out $4