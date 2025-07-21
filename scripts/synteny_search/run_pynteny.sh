#!/bin/bash

# source activate pynteny_env # synteny_wrapper.sh already activates the pynteny conda environment,
# but you should uncomment this line if running run_pynteny.sh individually.

### Inputs:
# $1: a genome database directory in which to perform the synteny searches.
# Using ${HOME} in the path is fine, but ~ seems to cause issues.

# $2: a TSV with two columns:
# Column 1 contains a synteny structure query for pynteny search
# Column 2 contains an outdir (use an absolute path!) to which the results from that pynteny search will be written

# $3: hmms_dir (the directory in which Pynteny should look for HMMs)

# $4: hmms_metadata (the metadata table associated with the HMMs in hmms_dir)

### Outputs:
# Phase 1: builds labelled peptide databases in each of the subdirs of $genome_db
# Phase 2: for each synteny structure query in column 1 of $input_file, writes the output of pynteny search
# to the corresponding outdir from column 2 of input file, in a "genomes" subdir of the outdir.

### Example runs:
# bash run_pynteny.sh -g "${HOME}/data/genome_db/" -d "${HOME}/data/synteny_input.tsv" -d "${HOME}/data/hmms/hmm_PGAP" -m "${HOME}/data/hmms/hmm_PGAP.tsv" &
# bash run_pynteny.sh -g "${HOME}/data/genbank_toy/bacteria/" -i "${HOME}/data/synteny_input_fha1.tsv" -d "${HOME}/data/hmms_fha1_search/hmms_fha1" -m "${HOME}/data/hmms_fha1_search/hmms_fha1.tsv"

### Old example runs:
# bash run_pynteny.sh "${HOME}/data/genome_db/" "${HOME}/data/synteny_input.tsv" &
# bash run_pynteny.sh "${HOME}/data/genbank_toy/bacteria/" "${HOME}/data/synteny_input_example.tsv"

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
hmms_dir=""
hmms_metadata=""

# Parse command-line flags
while getopts "g:i:h:d:m:" opt; do
  case $opt in
    g) genome_db="${OPTARG%/}/" ;;  # ensures the path to this dir ends with exactly one slash
    i) export input_file="$OPTARG" ;;
    d) export hmms_dir="$OPTARG" ;; # directory of HMMs for Pynteny to use in its synteny search
    m) export hmms_metadata="$OPTARG" ;; # metadata associated with aforementioned HMM directory
    \?) echo "Usage: $0 -g genome_db -i input_file -d hmms_dir -m hmms_metadata" >&2
        exit 1 ;;
  esac
done

# Input and output directories
#genome_db="${1%/}/" # ensures the path to this dir ends with exactly one slash
#input_file=$2

# export these as global variables
#export input_file=$2
# Pynteny wasn't able to find the following automatically, so I needed to specify them
#export hmms_dir=$3
#export hmms_metadata=$4


### Create output directories as necessary
IFS=$'\n' # lines of $input_file are newline-separated.
# Need to set up the loop this way, otherwise it only reads the first line
for next in $(cat ${input_file}); do
  outdir=$(echo $next | cut -f 2) # column 2 of $input_file contains outdirs
  mkdir -p "${outdir}/genomes" # keep outputs in a subdir so it's easier to find the output of synteny_wrapper.sh later
done


### Phase 1: Build peptide databases in parallel (only in subdirectories that don't yet have peptide databases)

build_database() {
    dir="$1" # subdirs of depth 1 from $genome_db

    if [ -d "$dir" ]; then
        labelled_db="$dir/labelled_peptides.faa" # "labelled_peptides.faa" is the labelled peptide database; it should be in $dir

        if [ ! -f "$labelled_db" ]; then # if it can't find the labelled peptide database, it will build it.
            # The code below uses for loops for convenience, but we're only expecting one annotation file per subdir of $genome_db
            # Even if there are multiple .gbff files, there will be only one labelled peptide database will be made because they'll overwrite each other
            for my_file in "$dir"/*.gbff.gz; do
                gzip -d "$my_file" # unzip so that pynteny build will work
            done

            for my_file in "$dir"/*.gbff; do
                echo "Building database from ${my_file} (PID $$)"
                pynteny build --data "$my_file" --outfile "$labelled_db"
                gzip "$my_file"
            done
        fi
    fi
}

export -f build_database

# -mindepth of 1 so it runs on $genome_db's subdirs but not itself
echo "=== PHASE 1: Building peptide databases ==="
find "$genome_db" -mindepth 1 -maxdepth 1 -type d -print0 | \
    parallel -0 -j $(nproc) --eta --progress --joblog build_db.log \
    'build_database {}'


### Phase 2: Run analyses in parallel (only in subdirectories that don't yet have peptide databases)
run_analyses() {
    dir="$1" # subdirs of depth 1 from $genome_db

    if [ -d "$dir" ]; then
        labelled_db="$dir/labelled_peptides.faa"

        if [ -f "$labelled_db" ]; then # Confirms that the labelled peptide database exists
            IFS=$'\n' # lines of $input_file are newline-separated.

            for next in $(cat ${input_file}); do # each line of $input_file represents a different synteny search query
                query=$(echo $next | cut -f 1) # column 1 of $input_file contains synteny structure queries
                outdir=$(echo $next | cut -f 2) # column 2 of $input_file contains outdirs

                # Run synteny analysis if needed
                out="${outdir}/genomes/$(basename "$dir")" # place results in the genomes subdir
                # if [ ! -f "${out}/synteny_matched.tsv" ]; then # turns out that if a database doesn't include an HMM, no synteny_matched.tsv file is written
                if [ ! -d "$out" ]; then # checks if the outdir has been made; in theory, if this exists, then pynteny search concluded successfully.
                    echo "Running ${query} analysis on $labelled_db (PID $$), writing to ${out}"
                    pynteny search \
                      --synteny_struc "${query}" \
                      --data "$labelled_db" \
                      --hmm_dir "$hmms_dir" \
                      --hmm_meta "$hmms_metadata" \
                      --outdir "$out"
                fi
            done
        else
            echo "$dir does not contain labelled_peptides.faa"
        fi
    fi
}

export -f run_analyses

# also turn these variables global so they can be used inside the function
#export input_file="$input_file"
# Pynteny wasn't able to find the following automatically, so I needed to specify them
# TODO: turn these into CLIs
#export hmms_dir="/home/kcw2/data/hmms/hmm_PGAP"
#export hmms_metadata="/home/kcw2/data/hmms/hmm_PGAP.tsv"

echo "=== PHASE 2: Running analyses ==="
find "$genome_db" -mindepth 1 -maxdepth 1 -type d -print0 | \
    parallel -0 -j "$(nproc)" --eta --progress --joblog run_analyses.log \
    'run_analyses {}'