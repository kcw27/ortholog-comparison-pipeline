#!/bin/bash

print_help() {
    echo "Usage: $0 -p <db_path> -d <db_name> -i <query> -o <outfile> -m <max> -s"
    echo ""
    echo "  -g    Path to genome DB from which BLAST will be created"
    echo "  -f    (optional) Output name for FASTA containing all sequences in the genome DB"
    echo "  -d    Path to directory to which BLAST db will be saved"
    echo "  -t    Title to create BLAST db with"
    echo "  -h    Show help message and exit"
}

# Set default value for optional flags
db_fasta="comprehensive.fasta" # will be created in db_dir

# Parse options with getopts (not GNU getopts, as that causes issues on Mac systems)
while getopts "g:f:d:t:h" opt; do
  case $opt in
    g) genome_db="$OPTARG" ;; 
    f) db_fasta="$OPTARG" ;;
    d) db_dir="$OPTARG" ;;
    t) db_title="$OPTARG" ;;
    h) print_help; exit 0 ;;
  esac
done

scriptsdir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # get location of current script; this is where the other scripts are found too

mkdir -p $db_dir
cd $db_dir

echo "Producing FASTA from genome DB..."
bash "${scriptsdir}/make_comprehensive_protein_fasta.sh" $genome_db $db_fasta

echo "Producing BLAST db..."
makeblastdb -in $db_fasta -input_type fasta -dbtype prot -title $db_title -out $db_title