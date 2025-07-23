#!/bin/bash

source activate pynteny_env # activate the conda environment for Pynteny

### Inputs:
# $1: a genome database directory in which to perform the synteny searches.

# $2: a TSV with two columns:
# Column 1 contains a synteny structure query for pynteny search
# Column 2 contains an outdir (use an absolute path!) to which the results from that pynteny search will be written

# $3: a newline-separated text file listing HMMs that represent genes of interest.
# e.g. if you're looking for metadata related to PA3565 orthologs, $3 contains PF00126 and PF03466,
# and any gene with either HMM will be written to the output file of find_synteny_hits.sh, synteny_summary.tsv 

### Example runs: submit as a job because it will take a while.
# bash synteny_wrapper.sh -g "${HOME}/data/fha1_dbs" -i "${HOME}/data/synteny_input_fha1.tsv" -h "${HOME}/data/hmms_of_interest_fha1.txt" -d "${HOME}/data/hmms_fha1_search/hmms_fha1" -m "${HOME}/data/hmms_fha1_search/hmms_fha1.tsv" &
# for testing:
# bash synteny_wrapper.sh -g "${HOME}/data/genbank_toy/bacteria/" -i "${HOME}/data/synteny_input_fha1.tsv" -h "${HOME}/data/hmms_of_interest_fha1.txt" -d "${HOME}/data/hmms_fha1_search/hmms_fha1" -m "${HOME}/data/hmms_fha1_search/hmms_fha1.tsv"

### Old example runs: 
# bash synteny_wrapper.sh "${HOME}/data/genome_db" "${HOME}/data/synteny_input.tsv" "${HOME}/data/hmms_of_interest.txt" &
# for testing:
# bash synteny_wrapper.sh "${HOME}/data/genbank_toy/bacteria/" "${HOME}/data/synteny_input_example.tsv" "${HOME}/data/hmms_of_interest.txt" &

### Input directory structure:
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
hmms_dir=""
hmms_metadata=""

# Parse command-line flags
while getopts "g:i:h:d:m:" opt; do
  case $opt in
    g) genome_db="${OPTARG%/}/" ;;  # ensures the path to this dir ends with exactly one slash
    i) input_file="$OPTARG" ;;
    h) hmms="$OPTARG" ;; # hmms of interest; only write metadata pertaining to these hmms
    d) hmms_dir="$OPTARG" ;; # directory of HMMs for Pynteny to use in its synteny search
    m) hmms_metadata="$OPTARG" ;; # metadata associated with aforementioned HMM directory
    \?) echo "Usage: $0 -g genome_db -i input_file -h hmms -d hmms_dir -m hmms_metadata" >&2
        exit 1 ;;
  esac
done

wrapperdir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # get location of current script
scriptsdir="${wrapperdir}/synteny_search" # this is where I put the scripts
log="${wrapperdir}/synteny_wrapper_log.txt"
echo "CLIs: ${genome_db}, ${input_file}, ${hmms}, ${hmms_dir}, ${hmms_metadata}" > $log

start_time="$(date -u +%s)"

# bash "${scriptsdir}/run_pynteny.sh" $genome_db $input_file $hmms_dir $hmms_metadata &
bash "${scriptsdir}/run_pynteny.sh" -g $genome_db -i $input_file -d $hmms_dir -m $hmms_metadata &
pid1=$!  # Capture process ID (PID) of run_pynteny.sh
echo "Run Pynteny to find synteny hits... Job ID ${pid1}, time $(date -u +%s)" >> $log
wait $pid1  # Stall execution until run_pynteny.sh finishes


#bash "${scriptsdir}/find_synteny_hits.sh" $genome_db $input_file $hmms &
bash "${scriptsdir}/find_synteny_hits.sh" -g $genome_db -i $input_file -h $hmms &
pid2=$!  # Capture PID of find_synteny_hits.sh
echo "Done running Pynteny! Summarizing Pynteny outputs... Job ID ${pid2}, time $(date -u +%s)" >> $log
wait $pid2  # Stall execution until find_synteny_hits.sh finishes

end_time="$(date -u +%s)"




elapsed="$(($end_time-$start_time))"
echo "Finished! Elapsed: ${elapsed} seconds." >> $log
echo "Start time: ${start_time}" >> $log
echo "End time: ${end_time}" >> $log