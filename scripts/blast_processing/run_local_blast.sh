#!/bin/bash

# example run:
# $blastdir/run_local_blast.sh -p "$rundir/toydb" -d "toydb" -i "$rundir/PAO1_PA3565.fasta" -o "$rundir/PA3565_hits.blast"

print_help() {
    echo "Usage: $0 -p <db_path> -d <db_name> -i <query> -o <outfile> -m <max>"
    echo ""
    echo "  -p    Path to BLAST db"
    echo "  -d    Name of BLAST db"
    echo "  -i    Path to BLAST query FASTA"
    echo "  -o    Name of output file"
    echo "  -m    Max target seqs to return per sequence in query FASTA"
    echo "  -h    Show help message and exit"
}

# Set default value for optional flag
max=500000

# Parse options with getopts (not GNU getopts, as that causes issues on Mac systems)
while getopts "p:d:i:o:m:h" opt; do
  case $opt in
    p) db_path="$OPTARG" ;;
    d) db_name="$OPTARG" ;;
    i) query="$OPTARG" ;;
    o) outfile="$OPTARG" ;;
    m) max="$OPTARG" ;;
    h) print_help; exit 0 ;;
    *) print_help; exit 1 ;; #echo "Invalid option"; exit 1 ;;
  esac
done

# Confirm that required arguments have been provided
if [[ -z "$db_path" || -z "$db_name" || -z "$query" || -z "$outfile" ]]; then
    echo "Error: Flags -p, -d, -i and -o are required." >&2
    print_help
    exit 1
fi

# run script
cd $db_path # so that it can find the blast database
procs_to_use=$(( $(nproc) / 4 ))
mkdir -p $(dirname $outfile) # make the outdir if it doesn't exist already
scriptsdir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # get location of current script; this is where the other scripts are found too

tempBlast=$(mktemp)
echo "Running BLAST... Temp file at $tempBlast"
blastp -db "$db_name" -query "$query" -max_target_seqs "$max" -show_gis -num_threads $procs_to_use -outfmt "6 sallgi sallseqid sseq evalue salltitles" -out $tempBlast

# add organisms column
tempBlastOrgs=$(mktemp)
echo "Adding organisms... Temp file at $tempBlastOrgs"
bash "${scriptsdir}/add_organism_column.sh" "$tempBlast" "$tempBlastOrgs"

# save BLAST output as long format
echo "Expanding BLAST output file to long format..."
bash "${scriptsdir}/expand_blast_output.sh" -b "$tempBlastOrgs" -o "$outfile"