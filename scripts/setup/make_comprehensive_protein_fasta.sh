#!/bin/bash

# inputs:
# $1: a directory like the output of kblin's ncbi-genome-download, in which the subdirectories are named after genomes, and each subdirectory contains one .gbff.gz file
# $2: name of fasta to which to write all proteins in the genomes. Names in fasta are genome_id-protein_id followed by

# Prior to running, please activate the conda environment for https://github.com/AstrobioMike/bit

genomedir=${1%/}
master_fasta=$2
> $master_fasta # wipes file clean if it already exists

totalgenomes=$(ls $genomedir | wc -l | cut -f 1 -d " ")

sum=0
genomecounter=1
for dir in "${genomedir}/"*; do
  genome=$(basename "$dir")
  echo "Genome $genomecounter out of $totalgenomes: $genome"

  gb=$(echo "$dir"/*.gbff.gz)
  gb=$(echo "${gb%% *}") # if for some reason there are multiple genbank files in a dir, only process the first one

  # # unzip the file, process it, rezip
  # # it'll complain if a file has already been unzipped, but ultimately it doesn't matter
  # gunzip $gb
  # gb_unzipped=$(echo "$dir"/*.gbff)

  # gunzip takes an unreasonably long time; zcat is faster
  gb_unzipped="${dir}/temp.gbff" # delete this later
  zcat $gb > $gb_unzipped
  
  bit-genbank-to-AA-seqs -i $gb_unzipped -o "$dir/temp_aa.fasta"
  organism=$(cat $gb_unzipped | grep "\bORGANISM\b" | head -n 1 | awk '{$1=$1;print}' | cut -d " " -f 2-)
  
  # update the fasta titles to follow the desired format   
  awk '
    /^>/ {
        old = substr($0, 2)                 # drop leading ">"
    
        n = split(old, a, "_")              # split on underscores
        last = a[n]                         # last field
    
        # Case 1: last field contains letters: use only last field as protein ID (e.g. EGH06075.1)
        if (last ~ /[A-Za-z]/) {
            protein = last
        }
        # Case 2: last field has no letters: use last two fields as protein ID (e.g. WP_003122138.1)
        else {
            protein = a[n-1] "_" a[n]
        }
    
        print ">" genome "-" protein " " old " [" organism "]"
        next
    }
    { print }
    ' genome="$genome" organism="$organism" "$dir/temp_aa.fasta" >> "$master_fasta"

  
  contrib=$(wc -l "$dir/temp_aa.fasta" | cut -f 1 -d " ") # line count from this gb file
  sum=$(echo $(( $sum + $contrib )) )

  rm -f "$dir/temp_aa.fasta"
  # gzip "$gb_unzipped"
  rm -f "${dir}/temp.gbff" # hard-coding this instead of referring to $gb_unzipped in case I revert the approach and forget to uncomment this
  
  genomecounter=$(( $genomecounter + 1 ))
done

echo "Line count that $master_fasta should have: $sum"
wc -l $master_fasta
# when I ran this on a very large dataset, the two numbers differed, but as far as I could tell, there was nothing wrong with the FASTA???