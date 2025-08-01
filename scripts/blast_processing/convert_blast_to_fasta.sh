#!/bin/bash

# Inputs:
# Mandatory: BLAST output file (-outfmt "6 sgi sseqid sseq evalue stitle")
# Optional: header (true if there's a header in the BLAST file, which the program will skip; false otherwise)
# Optional: extract (true if the sequence ID needs to be extracted from within | symbols, false otherwise)

# Output: FASTA file (with headers >ID|title|evalue)

# Example runs:
# on synteny-filtered output (-hd is required because it has a header):
# bash $blastscripts/convert_blast_to_fasta.sh -b $datadir/PA3565_67_topHitsPerGenome.blast -hd

# on blast output that hasn't been through metadata processing (i.e. blast2gen.py or intersect_blast_and_synteny_with_metadata.R) and therefore has no headers (-e)
# bash $blastscripts/get_blast_top_hits.sh $datadir/PA3565_with_orgs_long.blast $datadir/PA3565_top_byOrganism.blast

# on blast output that has been through a filtering script that left it without a header:
# /home/kcw2/ortholog-comparison-pipeline/test_data/PA3565_top_byGenomeID.fasta

# on blast2gen.py output which hasn't gone through further filtering (-hd and -e)
# bash $blastscripts/convert_blast_to_fasta.sh -b $datadir/PA3565_with_orgs_long_annotatedByBlast2gen.blast -e -hd

# set default values for optional arguments
header=false
extract=false

# Parse flags
while [[ $# -gt 0 ]]; do
  case "$1" in
    -b|--blast)
      blast="$2"
      shift 2
      ;;
    -hd|--header)
      header=true
      shift
      ;;
    -e|--extract)
      extract=true
      shift
      ;;
    -h|--help)
      echo "I will add a help string later"
      exit 1
      ;;
    *)
      echo "Unknown option: $1" >&2
      exit 1
      ;;
  esac
done

if [[ -z "$blast" ]]; then
  echo "Error: BLAST file not provided. Use -b or --blast to specify input." >&2
  exit 1
fi

output="${blast%.*}.fasta"
> "$output"

echo "Writing results to ${output}..."

IFS=$'\n'

# If there's a header in the file, skip it.
if [[ "$header" = true ]]; then
  tail_input=$(tail -n +2 "$blast")
else
  tail_input=$(cat "$blast")
fi

while read -r next; do
  # Extract columns 2, 5, 4, 3 in one go
  #read id title evalue sequence < <(awk -F '\t' '{print $2, $5, $4, $3}' <<< "$next")
  #IFS=$'\t' read -r col1 col2 col3 col4 col5 <<< "$next"
  IFS=$'\t' read -r col1 col2 col3 col4 col5 col_rest <<< "$next"
  # there may be more than 5 columns, so if there are, they get dumped into col_rest. (If there aren't, col_rest will just be empty.)
  
  id="$col2"
  title="$col5"
  evalue="$col4"
  sequence="$col3"

  if [[ "$extract" = true ]]; then
    id=$(cut -d "|" -f 2 <<< "$id")
  fi

#  echo "&&&"
#  printf ">%s %s|%s\n%s\n" "$id" "$title" "$evalue" "$sequence"
  printf ">%s %s|%s\n%s\n" "$id" "$title" "$evalue" "$sequence" >> "$output"
  
#  echo "***"
#  echo $id $title $evalue
#  echo $sequence
#  echo "///"
#  echo ">${id} ${title}|${evalue}"
#  echo "$sequence"
  #echo ">${id} ${title}|${evalue}" >> $output
  #echo "$sequence" >> $output # this is the sequence corresponding to the header
done <<< "$tail_input"


echo "Finished converting BLAST to FASTA: ${output}"


#### Example run:
#### bash convert_blast_to_fasta.sh ${HOME}/data/PA3565_orthologs_65_66_67_top_evalueThreshold_1e-50.txt
###
###blast=$1
###output="${blast%.*}.fasta"
###> $output
###echo "Writing results to ${output}..."
###
###IFS=$'\n' # input file is newline-separated.
#### Need to set up the loop this way, otherwise it only reads the first line
###
#### Serial implementation, old version:
###for next in $(cat $blast); do
###  id=$(echo $next | cut -f 2 | cut -d "|" -f 2) # extract the id from within the pipe symbols 
###  title=$(echo $next | cut -f 5)
###  evalue=$(echo $next | cut -f 4)
###  echo ">${id} ${title}|${evalue}" >> $output
###  echo $next | cut -f 3 >> $output # this is the sequence corresponding to the header
###done
###
#### Parallelization: I don't think this code actually writes the correct header-sequence pairs to the file, so use this with caution
####export output  # Ensure output file is accessible across subshells
###    
####cat "$blast" | xargs -I {} -P "$(nproc)" bash -c '
###    #id=$(echo "{}" | cut -f 2 | cut -d "|" -f 2)
###    #title=$(echo "{}" | cut -f 5)
###    #evalue=$(echo "{}" | cut -f 4)
###    #seq=$(echo "{}" | cut -f 3)
###    # need to use a single echo command to write the output to the file, otherwise the header might be separated from the sequence
###    #echo -e ">${id}|${title}|${evalue}\n${seq}" >> '"$output"' # assumes there aren't escape sequences in the text, which is risky
###
###
###
###echo "Finished converting BLAST to FASTA!"