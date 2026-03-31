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
    echo "Usage: $0 -b <blast> -s <synteny> -m <multiple> -p <pident> -q <qcovs> -c<'genome'|'locus'> -tempdir -o <outdir> -h"
    echo "        Outputs: to -o, writes the output of the pairwise BLAST of $synteny as query against $blast as subject."
    echo "                 Rationale for this direction of BLAST is explained in comments at the top of this script."
    echo "                 Also writes a fasta with -p and -q cutoff- look for a file with a name similar to outname."
    echo "  -b    Required. Path to BLAST file in which the first five columns are sgi sseqid sseq evalue stitle."
    echo "  -s    Required. Either a path to a single summary file from synteny search via synteny_wrapper.sh (which calls find_synteny_hits.sh),"
    echo "                  or a path to a synteny input file in which the listed outdirs are assumed to contain synteny_summary.tsv files."
    echo "  -m    Optional. If true, treats the -s input as a file referring to (potentially) multiple outdirs and thus multiple summary files (default false)."
    echo "  -p    Optional. Pident threshold for filtering, discarding any hits with pident lower than this (default = 99)."
    echo "  -q    Optional. qcovs threshold for filtering, discarding any hits with pident lower than this (default = 90)."
    echo "  -c    Optional. Concatenation mode for convert_blast_to_fasta.sh (default = 'no')."
    echo "  -k    Optional. Keep FASTA of all sequences in $synteny (default = 'false')."
    echo "  -t    Optional. Directory for temporary files (by default, use mktemp)"
    echo "  -o    Optional. Name of output directory (default $(pwd))."
    echo "  -h    Show help message and exit"
}

### Setup
# set default values for optional arguments
multiple="false"
export pident="99"
export qcovs="90"
export concat="no"
export keep="false"
export tempdir="" # if empty, will later use mktemp
export outdir=$(pwd)

# Parse options
while getopts "b:s:mp:q:c:kt:o:h" opt; do
  case $opt in
    b) export blast="$OPTARG" ;;
    s) synteny="$OPTARG" ;;
    m) multiple="true" ;;
    p) pident="$OPTARG" ;;
    q) qcovs="$OPTARG" ;;
    c) concat="$OPTARG" ;;
    k) keep="true" ;;
    t) tempdir="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    h) print_help; exit 0 ;;
    *) print_help; exit 1 ;; #echo "Invalid option"; exit 1 ;;
  esac
done

currdir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # get location of current script
export scriptsdir="${currdir}/../blast_processing" # this is where I put the scripts that this script calls
mkdir -p $outdir # make the outdir if it doesn't exist already

# if no tempdir specified, just use mktemp
if [[ "$tempdir" == "" ]]; then
  tempdir=$(mktemp -d)
fi
mkdir -p $tempdir

### Functions
make_blast_db_from_blastfile() { 
  # Input: global variables $blast and $tempdir
  # Output: BLAST db at ${tempdir}/temp_db/, titled temp_db
  echo "Will create BLAST db at ${tempdir}/temp_db/, named temp_db"

  echo "Processing $blast"
  tempBlast="${tempdir}/tempBlast.blast"
  tempBlastFasta="${tempdir}/tempBlast.fasta"
  echo "Saving tempBlastFasta to ${tempBlastFasta}..."

  # remove gap characters "-" because command line BLAST doesn't recognize them.
  # only write a line if the sequence hasn't been encountered before; this ensures the FASTA only contains unique seqs
  # remove gaps BEFORE checking for uniqueness- sequences with gaps may be identical without gaps
  awk 'BEGIN {OFS="\t"} {
    gsub("-", "", $3);
    if (!seen[$3]++) print
  }' "$blast" > "$tempBlast"
  # very important to set the delimiter as tab; awk defaults to space delimiter
  bash "${scriptsdir}/convert_blast_to_fasta.sh" -b "$tempBlast" -o "$tempBlastFasta"

  # assume it's a protein sequence database
  currWorkingDir=$(pwd) # save this; especially if running this script in Nextflow, need to be careful about current directory location
  cd $tempdir # need to be one directory above the blast db being made
  makeblastdb -in "$tempBlastFasta" -input_type fasta -dbtype prot -title temp_db -out "${tempdir}/temp_db/temp_db"
  cd $currWorkingDir

  # clean up temp files
  # rm -f $tempBlastFasta
  # rm -f $tempBlast
}


blast_synteny_against_db() {
  # Input: path to a single synteny search summary, optionally a subdir (useful if you have multiple synteny summary files)
  # Also, global variables $tempdir, $pident, $qcovs, $concat, $keep
  # Outputs: synteny FASTA filtered according to $pident and $qcovs (the primary output),
  # and additionally, pairwise BLAST file and (if $keep is true) the full FASTA for the synteny summary file
  syntenyFile=$1
  local subdir="${2:-""}"
  local blast_outname="${outdir}/${subdir}pairwiseSyntenyBlast_pident${pident}_qcovs${qcovs}.blast"
  local fasta_outname="${outdir}/${subdir}filteredSynteny_pident${pident}_qcovs${qcovs}.fasta"
  echo "blast outname $blast_outname"
  echo "fasta outname $fasta_outname"

  # Create a unique temp directory for this job so if running multiple jobs at once, temp files don't clobber each other
  local jobtmp
  jobtmp=$(mktemp -d "${tempdir}/jobXXXXXX")
  local tempSynteny="${jobtmp}/tempSynteny.tsv"
  local tempSyntenyFasta="${jobtmp}/tempSynteny.fasta"
  local tempBlastOutputs="${jobtmp}/tempBlastOutputs.blast"
  local tempFilteredIDs="${jobtmp}/tempFilteredIDs.txt"

  echo "Processing $syntenyFile"
  echo "Saving tempSyntenyFasta to ${tempSyntenyFasta}..."

  # need to reorder columns so it's compatible with convert_blast_to_fasta.sh
  awk -F'\t' 'BEGIN {OFS="\t"} {print $1, $7, $8, "na", $5, $6}' "$syntenyFile" > "$tempSynteny"
  bash "${scriptsdir}/convert_blast_to_fasta.sh" -b "$tempSynteny" -o "$tempSyntenyFasta" -c "$concat" -r # -r to remove header
  
  ### Perform the pairwise BLAST
  echo "Saving pairwise blast outputs..."
  procs_to_use=$(( $(nproc) / 4 ))
  max=10
  blastp -query "$tempSyntenyFasta" -db "${tempdir}/temp_db/temp_db" -max_target_seqs "$max" -num_threads "$procs_to_use" -outfmt "6 qseqid qseq sallseqid sseq pident qcovs evalue bitscore" -out "$tempBlastOutputs"

  ### Filter the output file from the pairwise blast to only include lines with >=$pident
  awk -v p="$pident" -v q="$qcovs" '$5 >= p && $6 >= q' "$tempBlastOutputs" | cut -f 1 | sort | uniq > "$tempFilteredIDs" # IDs were in column 1 of the blast that was just done, pident was in column 5

  ### Filter the synteny fasta to only include the lines from the filtered IDs list
  # This filtered FASTA is the final output of this script. You can still left join the synteny metadata to the seq IDs in that filtered FASTA; it's a subset.
  seqtk subseq "$tempSyntenyFasta" "$tempFilteredIDs" > "$fasta_outname"
  echo "Saved filtered FASTA to $fasta_outname"

  # Add colnames to the file saved to $outname
  awk 'BEGIN { OFS="\t"; print "qseqid", "qseq", "sallseqid", "sseq", "pident", "qcovs", "evalue", "bitscore" } { print }' "$tempBlastOutputs" > "$blast_outname"

  # if requested, keep a copy of the full FASTA too
  if [[ $keep == "true" ]]; then
    local syntenyFastaName="${outdir}/${subdir}synteny_full.fasta"
    cp $tempSyntenyFasta $syntenyFastaName
    echo "Saved synteny summary as FASTA to $syntenyFastaName"
  else
    rm -f $tempSyntenyFasta 
  fi
}

### Call functions
# Build blast DB
export blastdb_address="${tempdir}/temp_db/temp_db"
echo "BLAST db should be created at $blastdb_address"
make_blast_db_from_blastfile # call the function

# Blast synteny summary/summaries against the blast DB
# Check $multiple and decide: one run against $synteny, or parallelized runs over dirs stored in $synteny col2?
if [[ "$multiple" != "true" ]]; then
  # just need to run this once
  echo "Single synteny file: $synteny"
  blast_synteny_against_db "$synteny"
else
  echo "Multiple synteny files stored in column 2 of $synteny"
  pids=()

  while read -r dir; do
      currSynteny="${dir}/synteny_summary.tsv"
      echo "Current synteny file: $currSynteny"
      subdir=$(basename "$dir")
      echo "subdir name is $subdir"
      mkdir -p "${outdir}/${subdir}"

      blast_synteny_against_db "$currSynteny" "${subdir}/" &
      pids+=($!)
  done < <(cut -f2 "$synteny")

  wait "${pids[@]}"
fi

# move blast DB to outdir if keeping
if [[ "$keep" == "true" ]]; then
  mkdir -p ${outdir}/temp_db/
  mv -f ${tempdir}/temp_db/* ${outdir}/temp_db/
  echo "Saved temp BLAST db to ${outdir}/temp_db/"
fi

# delete remaining temp files
rm -rf $tempdir