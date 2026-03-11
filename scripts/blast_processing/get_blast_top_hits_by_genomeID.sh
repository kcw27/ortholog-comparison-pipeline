#!/bin/bash

# example run:
# $blastdir/get_blast_top_hits_by_genomeID.sh -r -f $rundir/path1/PA3565_hits_filtered_metadata_forGenomeIDfilteringdemo.blast

# set default values
REMOVEHEADER="false"
GENOMECOL="1"
FILES=()

print_help() {
    echo "Usage: $0 -r -c [col] -f file1 [file2 ...]"
    echo
    echo "Filters BLAST output files to include only the top record (lowest evalue, column 4) per genome ID."
    echo "In event of a tie within a genome ID, only the first record is kept."
    echo "Output files are named like so: ${file%.*}_topPerGenome.blast"
    echo
    echo "Options:"
    echo '  -r               Optional. Include this flag to remove the first row, assumed to be the header (default: false).'
    echo '  -c <col>         Optional. Specifies the column in which genome ID is found (default: column 1, consistent with the current version of blast2gen.py)'
    echo '  -f <files...>    Required. One or more BLAST output files with format -outfmt "6 sallgi sallseqid sseq evalue salltitles", and organism column added to the end.'
    echo '                   Please ensure that all files have the same formate (i.e. having a header or not, and the column number for genome ID).'
    echo "  -h               Show help message and exit"
}

# Parse options
while getopts "rc:f:h" opt; do # notice: no colon after r or h because they don't take arguments
    case "$opt" in
        r)
            REMOVEHEADER="true"
            ;;
        c)
            GENOMECOL="$OPTARG"
            ;;
        f) # Parse all options *before* -f
            # Stop option parsing here — everything after -f is a filename
            FILES+=("$OPTARG") # get the first file
            shift $((OPTIND - 1)) # shift off processed args
            FILES+=("$@") # remaining args are files
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



# Define functions to handle the file depending on whether there's a header to remove
process_with_header() { 
  # arguments:
  # $1: input file
  # $2: output file 
  head -n 1 $1 > $2 # write header to the output file first
  
  # tail is used to remove the header
  # The first sort is to ensure that the file is sorted by genome ID and then by evalue (because the input file may not be in order),
  # the second is to get the top hit for each genome ID.
  tail -n +2 "$1" | sort -t$'\t' -k"$GENOMECOL","$GENOMECOL" -k4,4g | sort -u -t$'\t' -k"$GENOMECOL","$GENOMECOL" >> "$2" # write the top hit (lowest evalue) per genome to the file
}

process_without_header() { 
  # arguments:
  # $1: input file
  # $2: output file 
  
  # The first sort is to ensure that the file is sorted by genome ID and then by evalue (because the input file may not be in order),
  # the second is to get the top hit for each genome ID.
  cat "$1" | sort -t$'\t' -k"$GENOMECOL","$GENOMECOL" -k4,4g | sort -u -t$'\t' -k"$GENOMECOL","$GENOMECOL" >> "$2" # write the top hit (lowest evalue) per genome to the file
}

process() {
  # arguments:
  # $1: input file
  # $2: output file
  # $3: REMOVEHEADER
  
  if [[ "$3" == "true" ]]; then
    echo "As specified by user: Removing header."
    process_with_header "$1" "$2"
  else
    echo "As specified by user: the file has no header to remove."
    process_without_header "$1" "$2"
  fi
}

echo "As specified by user: Genome ID in column ${GENOMECOL}."

# Process each file
for file in "${FILES[@]}"; do
    output="${file%.*}_topPerGenome.blast"
    > $output # if the file already exists, overwrite it
    
    # Sort commands: the first is to ensure that the file is sorted by organism and then by evalue (in increasing order of evalue),
    # the second is to get the top hit for each organism.
    # Since the file is sorted from smallest to biggest evalue,
    # the first hit for an organism is assumed the most significant.
    # If records are tied for lowest evalue for an organism, it only takes the first.
    process "$file" "$output" "$REMOVEHEADER"
    echo "Processed: $file -> $output"
done