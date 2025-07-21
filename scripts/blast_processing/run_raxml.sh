#!/bin/bash

# CLIs:
# $1: $outname (-n); suffix of tree names
# $2: $seqname (-s), i.e. the PHYLIP alignment from which to make a tree
# $3: $submodel (-m) substitution model, e.g. PROTGAMMABLOSUM62 for protein alignment, or GTRGAMMA for nucleotide alignment
# $4: $outdir (-w)


# Example run:
# bash run_raxml.sh "raxml_test.txt" "aligned.ph" "PROTGAMMABLOSUM62" "/home/kcw2/data/raxml_test" &
# bash run_raxml.sh "raxml_65_66_67.txt" "/home/kcw2/data/PA3565_orthologs_65_66_67_top_evalueThreshold_1e-50.ph" "PROTGAMMABLOSUM62" "/home/kcw2/data/raxml_65_66_67" &

# need to activate the roary conda env before running 
#source /home/carmbrus/miniconda3/etc/profile.d/conda.sh
#conda activate /home/juneq/.conda/envs/roary
# source activate /home/juneq/.conda/envs/roary # borrowing someone else's conda environment that includes a RAxML installation, is this necessary?

cd /home/kcw2 # my home directory

outname=$1 # filename suffix for outputs
seqname=$2 # PHYLIP file
submodel=$3 # substitution model
outdir="${4%/}" # all output files will go here; remove trailing slash if present

mkdir -p "$outdir" # must exist prior to running RAxML

# sequential run in case parallel doesn't work
#raxmlHPC -x 2421 -f a -m "$submodel" -n $outname -p 748 -s $seqname -w $outdir -N 100

# Run RAxML with 100 bootstraps and 16 threads
# The call below is supposed to be parallel, but is sequential?
raxmlHPC-PTHREADS-SSE3 \
  -x 2421 \
  -p 748 \
  -f a \
  -m "$submodel" \
  -s "$seqname" \
  -n "$outname" \
  -T 16 \
  -N 100 \
  -w "$outdir"