#!/bin/bash

# Example run:
# bash find_synteny_hits.sh -g "${HOME}/data/fha1_dbs" -i "${HOME}/data/synteny_input_fha1.tsv" -L "${HOME}/data/hmms_of_interest_fha1.txt" &
# for testing:
# bash find_synteny_hits.sh -g "${HOME}/data/genbank_toy/bacteria/" -i "${HOME}/data/synteny_input_example.tsv" -L "${HOME}/data/hmms_of_interest.txt"


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
hmms_list=""

# Parse command-line flags
while getopts "g:i:L:" opt; do
  case $opt in
    g) export genome_db="${OPTARG%/}/" ;;  # ensures the path to this dir ends with exactly one slash
    i) export input_file="$OPTARG" ;;
    L) export hmms_list="$OPTARG" ;; # hmms of interest; only write metadata pertaining to these hmms
    \?) echo "Usage: $0 -g genome_db -i input_file -h hmms" >&2
        exit 1 ;;
  esac
done


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
    #contig_ids=$(grep -E "$(cat $hmms | grep -v "^$" | paste -sd '|')" "${output_genome}/synteny_matched.tsv" | cut -f 1)
    #locus_tags=$(grep -E "$(cat $hmms | grep -v "^$" | paste -sd '|')" "${output_genome}/synteny_matched.tsv" | cut -f 2)
    
    # now require all items of $hmms_list to be observed on the same line in order to report a match
    matches=$(awk '
      BEGIN {
        while ((getline < "'"$hmms_list"'") > 0) {
          if ($0 != "") patterns[++n] = $0
        }
      }
      {
        for (i = 1; i <= n; i++) {
          if (index($0, patterns[i]) == 0) next
        }
        print
      }
    ' "${output_genome}/synteny_matched.tsv")
    
    contig_ids=$(echo "$matches" | cut -f 1)
    locus_tags=$(echo "$matches" | cut -f 2)

    
    # Only do the following if there were matches found, i.e. locus_tags (and therefore contig_ids as well) is not empty
    if [[ -n "${locus_tags}" ]]; then
      # In $database_genome, search one level down for first matching .gbff* file
      genbank_file=$(find "$database_genome" -mindepth 1 -maxdepth 1 -name "*.gbff*" -print -quit)
    
      # First, get the information that is the same across all loci in locus_tags.
      organism=$(zcat $genbank_file | grep "\bORGANISM\b" | head -n 1 | awk '{$1=$1;print}' | cut -d " " -f 2-)
      
      # filter out extraneous lines and remove duplicates:
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
      
      isolation_site=$(zcat "$genbank_file" |
        awk '
          /Isolation Site/ {
            if (done) next   # ignore additional matches
            inblock = 1
      
            line = $0
            sub(/^[[:space:]]+/, "", line)
            sub(/.*::[[:space:]]*/, "", line)
            print line
            next
          }
      
          inblock {
            if ($0 ~ /^ {12}[^[:space:]]/) {
              inblock = 0
              done = 1
              next
            }
      
            line = $0
            sub(/^[[:space:]]+/, "", line)
            sub(/.*::[[:space:]]*/, "", line)
            print line
          }
        ' |
        sed ':a;N;$!ba;s/\n/ /g' |
        tr '\t' ' ' |      # <-- new: convert tabs to spaces
        tr -s ' '          # squeeze spaces
      )

      # get the sequencing technology
      seq_tech=$(zcat $genbank_file | awk '$0 ~/Sequencing Technology/ { sub(/^[^:]*[: ]+/, ""); print }' | head -n 1 )
      
      # get assembly method and coverage too
      asm_method=$(zcat $genbank_file | awk '$0 ~/Assembly Method/ { sub(/^[^:]*[: ]+/, ""); print }' | head -n 1 )
      genome_coverage=$(zcat $genbank_file | awk '$0 ~/Genome Coverage/ { sub(/^[^:]*[: ]+/, ""); print }' | head -n 1)
      
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
        echo "$genome_id	$contig	$organism	$isolation_source	$titles	$locus	$protein_id	$sequence	$seq_tech	$asm_method	$genome_coverage	$isolation_site" >> "${outdir}/synteny_summary.tsv"
      done
    fi
  fi
}

# Iterate over column 2 of the input file and parallelize metadata extraction
cut -f2 "$input_file" | while read -r outdir; do
  outdir=${outdir%/}
  echo "genome_id	contig	organism	isolation_source	titles	locus_tag	protein_id	sequence	seq_tech	assembly_method	genome_coverage	isolation_site" > "${outdir}/synteny_summary.tsv"  # Create a blank summary file for this synteny structure
  # renaming the "locus" column to "locus_tag" for consistency with blast2gen.py, which has been updated to include a locus tag column.
  export -f write_metadata
  
  # Go through directories two levels down from the Pynteny output directory to access outputs for each genome,
  # because the directory one level down is just the "genomes" directory created to organize them.
  # But also pass $outdir itself as an argument so that the corresponding synteny_summary.tsv file can be found.
  find "$outdir" -mindepth 2 -maxdepth 2 -type d | parallel -j $(nproc)  --eta --progress --joblog write_metadata.log 'write_metadata {}' "$outdir"
  
done


