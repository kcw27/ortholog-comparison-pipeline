#!/bin/bash

# Example run:
# bash find_synteny_hits.sh -g "${HOME}/data/fha1_dbs" -i "${HOME}/data/synteny_input_fha1.tsv" -h "${HOME}/data/hmms_of_interest_fha1.txt" &
# for testing:
# bash find_synteny_hits.sh -g "${HOME}/data/genbank_toy/bacteria/" -i "${HOME}/data/synteny_input_example.tsv" -h "${HOME}/data/hmms_of_interest.txt"

# Old example run:
# a run with actual data, submitted as a job with &
# bash find_synteny_hits.sh "${HOME}/data/genome_db" "${HOME}/data/synteny_input.tsv" "${HOME}/data/hmms_of_interest.txt" &
# for testing:
# bash find_synteny_hits.sh "${HOME}/data/genbank_toy/bacteria/" "${HOME}/data/synteny_input_example.tsv" "${HOME}/data/hmms_of_interest.txt"


#outfile="${1%/}/synteny_summary.tsv"

# Create output file once- no headers, to make it easier to parse later
# If there were header names, they would be:
# Genome ID, Contig ID, Organism, Isolation source list, Title list, Locus tag, Protein ID, Protein sequence
# The last two columns pertain to an individual protein from synteny_matched.tsv that matches 
# NEW: as the very last column (after protein sequence), added a column for the sequencing technology

#> "$outfile"

#!/bin/bash

# CLIs:
# $1 - genome_db: absolute path to directory whose subdirectories contain .gbff files
# $2 - input_file: absolute path to file whose second column has Pynteny outdirs
# $3 - file containing names of HMMs of interest
# Must use names from either the first column or second column of the hmm metadata file that you passed to the pynteny search call

### genome_db structure:
# genome_db="~/data/genome_db/" # its direct children are subdirectories named after accessions

# - $genome_db
# -- GCF_001646865.1
# --- GCF_001646865.1_rk21_genomic.gbff.gz 
# -- GCF_001646945.1
# --- GCF_001646945.1_AVMART05_1.0_genomic.gbff.gz 
# -- GCF_001647025.1
# --- GCF_001647025.1_EnAtr_2.0_genomic.gbff.gz 

### Output directory structure (for a given row in $input_file; each row in $input_file should have a unique outdir)
# query="<(PF00126|PF03466) 1  >(PF08240|PF13602)"
# outdir="~/data/results_65_67"

# - $outdir
# - synteny_summary.tsv (will be written by find_synteny_hits.sh)
# -- genomes
# --- GCF_001646865.1
# ---- synteny_matched.tsv
# ---- Other output files and directories
# --- GCF_001646945.1
# ---- synteny_matched.tsv
# ---- Other output files and directories
# --- GCF_001647025.1
# ---- synteny_matched.tsv
# ---- Other output files and directories

# Initialize variables
genome_db=""
input_file=""
hmms=""

# Parse command-line flags
while getopts "g:i:h:d:m:" opt; do
  case $opt in
    g) export genome_db="${OPTARG%/}/" ;;  # ensures the path to this dir ends with exactly one slash
    i) export input_file="$OPTARG" ;;
    h) export hmms="$OPTARG" ;; # hmms of interest; only write metadata pertaining to these hmms
    \?) echo "Usage: $0 -g genome_db -i input_file -h hmms" >&2
        exit 1 ;;
  esac
done


#export genome_db=${1%/}
#export input_file=${2%/}
#export hmms=$3

fname="synteny_matched.tsv" # Pynteny output file

# Function to extract metadata
write_metadata() {
  output_genome=${1%/}
  outdir=${2%/}  # $outdir was passed explicitly into the function
  genome_id=$(basename "$output_genome")
  database_genome="${genome_db}/${genome_id}"
  
  #echo "Processing: $output_genome"
  #echo "$output_genome $genome_id $database_genome"
  #cd $database_genome
  #ls
  
  # Only need to get metadata if synteny_matched.tsv exists.
  # synteny_matched.tsv does not exist for genomes that didn't contain one or more of the HMMs in the pynteny search query from run_pynteny.sh.
  if [[ -e "${output_genome}/synteny_matched.tsv" ]]; then
    # From synteny_matched.tsv, get the lists of contig_ids and locus_tags corresponding to any HMMs listed in $hmms
    contig_ids=$(grep -E "$(cat $hmms | grep -v "^$" | paste -sd '|')" "${output_genome}/synteny_matched.tsv" | cut -f 1)
    locus_tags=$(grep -E "$(cat $hmms | grep -v "^$" | paste -sd '|')" "${output_genome}/synteny_matched.tsv" | cut -f 2)
    
    # Only do the following if there were matches found, i.e. locus_tags (and therefore contig_ids as well) is not empty
    if [[ -n "${locus_tags}" ]]; then
      # In $database_genome, search one level down for first matching .gbff* file
      genbank_file=$(find "$database_genome" -mindepth 1 -maxdepth 1 -name "*.gbff*" -print -quit)
    
      # First, get the information that is the same across all loci in locus_tags.
      organism=$(zcat $genbank_file | grep "\bORGANISM\b" | head -n 1 | awk '{$1=$1;print}' | cut -d " " -f 2-)
      
      # old code for isolation_source and titles that only retrieved the first line of each hit
      #isolation_source=$(zcat $genbank_file | grep "isolation_source" | sort | uniq | tr -d \\n | tr -s " ") # may be a list, but if it exists, there's usually just one unique value
      #titles=$(zcat $genbank_file | grep "\bTITLE\b" | sort | uniq | tr -d \\n | tr -s " ") # a list, though the elements are probably just "Direct Submission" and an actual title
      
      # new code that gets all lines of each isolation_source and titles hit
      #isolation_source=$(zcat $genbank_file | sed -n '/isolation_source="/,/"/p' | tr -d '\n' | tr -s " ")
      # remove all newlines, then squeeze spaces
      #titles=$(zcat $genbank_file | sed -n '/\bTITLE\b/,/\bJOURNAL\b/p' | grep -v "JOURNAL" | tr -d '\n' | tr -s " ")
      # finds matches between lines featuring TITLE and JOURNAL
      # (with word boundaries), then uses inverse grep to remove the lines featuring JOURNAL, then removes newlines and squeezes spaces
      
      # further-updated code that is supposed to filter out extraneous lines and remove duplicates:
      # after finding matches and collapsing into a single line, reintroduce newlines (only between separate entries, not within entries)
      # then grep for the desired header (e.g. "^/isolation_source" to only keep lines that start with /isolation_source;
      # for some reason, sed include lines that came directly after the intended isolation source match), and use sort and uniq.
      isolation_source=$(zcat $genbank_file | sed -n '/isolation_source="/,/"/p' | tr -d '\n' | tr -s " " \
        | sed 's|" /|"\n/|g' | grep "^/isolation_source" | sort | uniq | tr '\n' ' ')
      # remove all newlines, then squeeze spaces, then reintroduce newlines only between separate entries to use sort and uniq,
      # then remove newlines again by replacing with spaces
      titles=$(zcat $genbank_file | sed -n '/\bTITLE\b/,/\bJOURNAL\b/p' | grep -v "JOURNAL" | tr -d '\n' | tr -s " " \
        | sed 's|TITLE|\nTITLE|g' | grep "^TITLE" | sort | uniq | tr '\n' ' ')
      # finds matches between lines featuring TITLE and JOURNAL
      # (with word boundaries), then uses inverse grep to remove the lines featuring JOURNAL, then removes newlines and squeezes spaces
      # then reintroduce newlines only between separate entries to use sort and uniq,
      # then remove newlines again by replacing with spaces
      
      # get the sequencing technology
      seq_tech=$(zcat $genbank_file | grep -oP '(?<=Sequencing Technology  :: ).*' | head -n 1)
      
      # Now iterate through contig_ids and locus_tags simultaneously, pulling the corresponding metadata and then writing it to the output file
      paste <(echo "$contig_ids" | tr ' ' '\n') <(echo "$locus_tags" | tr ' ' '\n') | while read contig locus; do
       # echo  Contig ID: $contig, Locus Tag: $locus" # write $contig and $locus to the file later
      
        # the many awk commands are to strip unnecessary text from the result
        protein_id=$(zcat "$genbank_file" | awk -v tag="$locus" '$0 ~ tag, /  gene  /' | \
          awk '/protein_id="/ {flag=1; sub(/.*protein_id="/, "")} flag && /"/ {flag=0; sub(/".*/, ""); print} flag' | \
          tr -d '[:space:]')
        sequence=$(zcat "$genbank_file" | awk -v tag="$locus" '$0 ~ tag, /  gene  /' | \
          awk '/translation="/ {flag=1; sub(/.*translation="/, "")} flag && /"/ {flag=0; sub(/".*/, ""); print} flag' | \
          tr -d '[:space:]') # in the .gbff, the sequence was split over multiple lines; this combines it into a single line

        # write metadata for this protein to the summary file, tab-separated
        #echo "$genome_id	$contig	$organism	$isolation_source	$title	$locus	$protein_id	$sequence"
        #echo "Should write to ${outdir}/synteny_summary.tsv"
        echo "$genome_id	$contig	$organism	$isolation_source	$titles	$locus	$protein_id	$sequence	$seq_tech" >> "${outdir}/synteny_summary.tsv"
      done
    #else
      #echo "No hits found for ${output_genome}"
    fi
    
    #echo $genbank_file
  #else
    #echo "$output_genome has no synteny_matched.tsv"
  fi
  
  #if [[ -n "$gbff_file" ]]; then
    #zcat "$gbff_file" | grep "protein_id" >> "${outdir}/synteny_summary.tsv"
    #echo "$output_genome ${outdir}/synteny_summary.tsv"
  #fi
}

# Iterate over column 2 of the input file and parallelize metadata extraction
cut -f2 "$input_file" | while read -r outdir; do
  outdir=${outdir%/}
  > "${outdir}/synteny_summary.tsv"  # Create a blank summary file for this synteny structure
  
  export -f write_metadata
  
  # Go through directories two levels down from the Pynteny output directory to access outputs for each genome,
  # because the directory one level down is just the "genomes" directory created to organize them.
  # But also pass $outdir itself as an argument so that the corresponding synteny_summary.tsv file can be found.
  find "$outdir" -mindepth 2 -maxdepth 2 -type d | parallel -j $(nproc)  --eta --progress --joblog write_metadata.log 'write_metadata {}' "$outdir"
  
done


