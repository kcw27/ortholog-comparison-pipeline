#!/bin/bash

# CLIs:
# $1: $multifasta, a FASTA of orthologs to align
# $2: $seqtype, either "protein" or "dna"
# $3: $outdir in which to save outputs. Subdirectories of $outdir will be alignment/ and trees/.

# Main outputs:
# alignment in FASTA and PHYLIP formats (in alignment/ subdirectory)
# _bipartition tree (in trees/ subdirectory)

# Example run:
# bash alignment_and_tree_wrapper.sh "${HOME}/data/PA3565_orthologs_65_66_67_top_evalueThreshold_1e-50.fasta" "protein" "${HOME}/data/PA3565_align_and_tree/PA3565_66_67" &

start_time="$(date -u +%s)"

### Setup
# Assign variable names
multifasta=$1
seqtype=$2
outdir="${3%/}"
mkdir -p $outdir

# Get root name of $multifasta
rootname=$(basename $multifasta | cut -d . -f 1)
# echo $rootname

# Get synonyms of $seqtype for function calls
# For ClustalW, just convert it to uppercase
seqtype_clustalw=$(echo "$seqtype" | tr '[:lower:]' '[:upper:]')

# For RAxML, pick an appropriate substitution model based on seqtype
if [[ "$seqtype" == "protein" ]]; then
    submodel="PROTGAMMABLOSUM62"
elif [[ "$seqtype" == "dna" ]]; then
    submodel="GTRGAMMA"
else
    echo "Invalid seqtype; must be either protein or dna"
    exit 1
fi
# echo $submodel

# Make alignment and tree subdirs of $outdir
aligndir="${outdir}/alignment"
treedir="${outdir}/trees"
mkdir -p $aligndir
mkdir -p $treedir

### Run scripts
scriptsdir="$HOME" # this is where I put the scripts
log="${outdir}/alignment_and_tree_wrapper_log.txt"
echo "CLIs for alignment_and_tree_wrapper.sh: ${multifasta}, ${seqtype}, ${outdir}" > $log

# Run run_clustalw.sh (if alignment doesn't already exist)
aligned_fasta="${aligndir}/aligned.fasta" # output of run_clustalw.sh

if [[ -f "$aligned_fasta" ]]; then
  echo "Alignment already exists at ${aligned_fasta}. Skipping ClustalW alignment." >> $log
else
  bash "${scriptsdir}/run_clustalw.sh" "$multifasta" "$seqtype_clustalw" "FASTA" "$aligned_fasta" &
  pid1=$! # capture process ID of the job
  echo "Run run_clustalw.sh with CLIs ${multifasta} ${seqtype_clustalw} FASTA ${aligned_fasta}; Job ID ${pid1}, time $(date -u +%s)" >> $log
  wait $pid1  # Stall execution until command finishes
fi

# Run fasta_to_phylip.sh
bash "${scriptsdir}/fasta_to_phylip.sh" "$aligned_fasta" "$seqtype" "$aligndir" &
pid2=$! # capture process ID of the job
echo "Run fasta_to_phylip.sh with CLIs ${aligned_fasta} ${seqtype} ${aligndir}; Job ID ${pid2}, time $(date -u +%s)" >> $log
wait $pid2  # Stall execution until command finishes
aligned_phylip="${aligndir}/aligned.phy" # PHYLIP output is aligned.phy in $aligndir

# Run run_raxml.sh
bash "${scriptsdir}/run_raxml.sh" "${rootname}.txt" "$aligned_phylip" "$submodel" "$treedir" &
pid3=$! # capture process ID of the job
echo "Run run_raxml.sh with CLIs ${rootname} ${aligned_phylip} ${submodel} ${treedir}; Job ID ${pid3}, time $(date -u +%s)" >> $log
wait $pid3  # Stall execution until command finishes

end_time="$(date -u +%s)"
elapsed="$(($end_time-$start_time))"
echo "Finished! Elapsed: ${elapsed} seconds." >> $log
echo "Start time: ${start_time}" >> $log
echo "End time: ${end_time}" >> $log