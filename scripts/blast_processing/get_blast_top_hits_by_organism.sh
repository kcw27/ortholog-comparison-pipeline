# example run:
# $blastdir/get_blast_top_hits_by_organism.sh -f $rundir/path1/PA3565_hits_evalueThreshold_1e-100.blast

FILES=()

print_help() {
    echo "Usage: $0 -f file1 [file2 ...]"
    echo
    echo "Filters BLAST output files to include only the top record (lowest evalue, column 4) per organism name (column 6)."
    echo "In event of a tie within an organism name, only the first record is kept."
    echo "Output files are named like so: ${file%.*}_topPerOrganism.blast"
    echo
    echo "Options:"
    echo '  -f <files...>    Required. One or more BLAST output files with format -outfmt "6 sallgi sallseqid sseq evalue salltitles", and organism column added to the end.'
    echo "  -h               Show help message and exit"
}

# Parse options
while getopts "f:h" opt; do
    case "$opt" in
        f)
            FILES=("$OPTARG")              # get the first file
            shift $((OPTIND - 1))          # shift off processed args
            FILES+=("$@")                  # remaining args are files
            break
            ;;
        h)
            print_help
            exit 0
            ;;
        *)
            print_help
            exit 1
            ;;
    esac
done

# Process each file
for file in "${FILES[@]}"; do
    output="${file%.*}_topPerOrganism.blast"
    > $output # if the file already exists, overwrite it
    
    # Sort commands: the first is to ensure that the file is sorted by organism and then by evalue (in increasing order of evalue),
    # the second is to get the top hit for each organism.
    # Since the file is sorted from smallest to biggest evalue,
    # the first hit for an organism is assumed the most significant.
    # If records are tied for lowest evalue for an organism, it only takes the first.
    cat "$file" | sort -t$'\t' -k6,6 -k4,4g | sort -u -t$'\t' -k6,6 > "$output"
    echo "Processed: $file -> $output"
done