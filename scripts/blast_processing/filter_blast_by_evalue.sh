#!/bin/bash

# Example run:
# $blastdir/filter_blast_by_evalue.sh -t "1e-100" -f $rundir/path1/PA3565_hits.blast

# Default threshold
THRESHOLD="1e-30"
FILES=()

print_help() {
    echo "Usage: $0 [-t threshold] -f file1 [file2 ...]"
    echo
    echo "Filters BLAST output files to include only rows where the e-value"
    echo "(column 4) is smaller than the provided threshold."
    echo "Output files are named like so: ${file%.*}_evalueThreshold_${THRESHOLD}.blast"
    echo
    echo "Options:"
    echo "  -t <threshold>   Optional. Default: 1e-30"
    echo '  -f <files...>    Required. One or more BLAST output files with format -outfmt "6 sallgi sallseqid sseq evalue salltitles", and organism column added to the end.'
    echo "  -h               Show help message and exit"
}

# Parse options
while getopts "t:f:h" opt; do
    case "$opt" in
        t)
            THRESHOLD="$OPTARG"
            ;;
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
    output="${file%.*}_evalueThreshold_${THRESHOLD}.blast"
    
    awk -F'\t' -v threshold="$THRESHOLD" '$4 <= threshold' "$file" > "$output" # -v is the safest way to pass a variable to awk
    echo "Processed: $file -> $output"
done