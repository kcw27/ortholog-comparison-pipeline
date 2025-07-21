#!/bin/bash

# Input: BLAST output file (-outfmt "6 sgi sseqid sseq evalue stitle")

# Output: FASTA file (with headers >ID|title|evalue)

# Example run:
# bash convert_blast_to_fasta.sh ${HOME}/data/PA3565_orthologs_65_66_67_top_evalueThreshold_1e-50.txt

blast=$1
output="${blast%.*}.fasta"
> $output
echo "Writing results to ${output}..."

IFS=$'\n' # input file is newline-separated.
# Need to set up the loop this way, otherwise it only reads the first line

# Serial implementation:
for next in $(cat $blast); do
  id=$(echo $next | cut -f 2 | cut -d "|" -f 2) # extract the id from within the pipe symbols 
  title=$(echo $next | cut -f 5)
  evalue=$(echo $next | cut -f 4)
  echo ">${id} ${title}|${evalue}" >> $output
  echo $next | cut -f 3 >> $output # this is the sequence corresponding to the header
done

# Parallelization: I don't think this code actually writes the correct header-sequence pairs to the file, so use this with caution
#export output  # Ensure output file is accessible across subshells
    
#cat "$blast" | xargs -I {} -P "$(nproc)" bash -c '
    #id=$(echo "{}" | cut -f 2 | cut -d "|" -f 2)
    #title=$(echo "{}" | cut -f 5)
    #evalue=$(echo "{}" | cut -f 4)
    #seq=$(echo "{}" | cut -f 3)
    # need to use a single echo command to write the output to the file, otherwise the header might be separated from the sequence
    #echo -e ">${id}|${title}|${evalue}\n${seq}" >> '"$output"' # assumes there aren't escape sequences in the text, which is risky



echo "Finished converting BLAST to FASTA!"