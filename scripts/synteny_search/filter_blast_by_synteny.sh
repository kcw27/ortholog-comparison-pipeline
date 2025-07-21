#!/bin/bash

# CAUTION: parallel execution causes this to run the risk of interleaved outputs.
# In my own experience, it hasn't happened for this script. If your output has the wrong number of columns,
# uncomment the non-parallel version and comment out the parallel version.

# CLIs:
# $1: BLAST output file to filter. Produced by command line BLAST, with -outfmt "6 sgi sseqid sseq evalue stitle"
# IMPORTANT: to increase the likelihood that the parallel calls write to the output file correctly, the protein IDs in $1 must be unique,
# i.e. grep only writes one line at a time.
# This should generally be the case for BLAST output files anyway.
# Based on filter_metadata_by_locus_tags.sh, I think interleaving in the output tends to happen when grep returns multiple lines.
# $2: the same $input_file as was used for the synteny analysis. Will get the synteny outdirs from column 2, and write to these outdirs.

# Output: in each synteny analysis outdir specified in $2, produces a synteny-filtered BLAST output file.
# That is, the BLAST file saved to the outdir only contains records whose locus tags appear in synteny_summary.tsv

# Example run:
# bash filter_blast_by_synteny.sh /home/kcw2/data/PA3565_orthologs.txt /home/kcw2/data/synteny_input.tsv

blast=$1
input_file=$2

IFS=$'\n' # $input_file is newline-separated.
# Need to set up the loop this way, otherwise it only reads the first line
for next in $(cat "$input_file" | cut -f 2); do # column 2 of $input_file contains synteny outdirs
  echo "Intersecting BLAST hits with synteny hits from ${next}..."
  synteny_summary="${next}/synteny_summary.tsv"
  
  name=$(basename "$blast") # to name the filtered output after the input file
  output_file="${next}/${name%.*}_synteny_filtered.tsv" # remove the extension from the basename, assuming extension starts at the first period
  # For information on which synteny structure was used to produce the filtered file, refer to the directory where the file was saved, as the name itself isn't informative.
  > $output_file # create the output file
  
  # non-parallel version if you want to be safe
  #for pid in $(cat "$synteny_summary" | cut -f 7); do # column 7 of $synteny_summary contains protein IDs
    #grep $pid $blast >> $output_file
  #done
  
  # Parallelize with xargs, going through unique protein IDs in $blast (which is column 7
  #cat "$synteny_summary" | cut -f 7 | xargs -P $(nproc) -I {} grep {} "$blast" >> "$output_file"
  cat "$synteny_summary" | cut -f 7 | sort | uniq | xargs -P $(nproc) -I {} grep {} "$blast" >> "$output_file" # use sort and uniq to remove duplicate protein IDs
  
  

  echo "Finished writing to ${output_file}"
done

echo "Done!"