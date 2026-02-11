# Description:
# (long format- one record per row)

### Example runs:
# with organisms column, default outname:
# bash "/home/kcw2/ortholog-comparison-pipeline/scripts/blast_processing/expand_blast_output.sh" -b "${HOME}/data/testing/out/PA3565_nr_small_orgs.txt"

# without organisms column, with a specified outname: 
# bash "/home/kcw2/ortholog-comparison-pipeline/scripts/blast_processing/expand_blast_output.sh" -b "/home/kcw2/data/testing/PA3565_nr_small.txt" -o "/home/kcw2/data/testing/PA3565_nr_small_foo.txt"


print_help() {
      echo "Usage: $0 -b <blast_file> [-o <output_file>]"
      echo "  -b     Mandatory: input BLAST file with 5 or 6 columns."
      echo "-outfmt 6 with arguments sallgi sallseqid sseq evalue salltitles (may have organisms column appended)."
      echo "  -o     Optional: output filename (default: <blast>_long.blast)"
      exit 0
}


# Parse options
while getopts "b:o:" opt; do
  case $opt in
    b) blast="$OPTARG" ;;
    o) output="$OPTARG" ;;
    h) print_help; exit 0 ;;
    *) print_help; exit 1 ;;
  esac
done

# Set default filename only if not set
outfile="${output:-${blast/.*}_long.blast}"
mkdir -p $(dirname $outfile) # make the outdir if it doesn't exist already

echo "Saving to $outfile..."


# handle file differently if organisms column is or isn't present
awk -F'\t' '
function process_with_organisms() {
    split($1, gis, /;/)
    split($2, seqids, /;/)
    split($5, titles, /<>/)
    split($6, organisms, /;/)
    n = length(gis)

    if (length(seqids) != n || length(titles) != n || length(organisms) != n) {
        print "ERROR: inconsistent field counts on line " NR > "/dev/stderr"
        next
    }

    for (i = 1; i <= n; i++) {
        print gis[i] "\t" seqids[i] "\t" $3 "\t" $4 "\t" titles[i] "\t" organisms[i]
    }
}

function process_without_organisms() {
    split($1, gis, /;/)
    split($2, seqids, /;/)
    split($5, titles, /<>/)
    n = length(gis)

    if (length(seqids) != n || length(titles) != n) {
        print "ERROR: inconsistent field counts on line " NR > "/dev/stderr"
        next
    }

    for (i = 1; i <= n; i++) {
        print gis[i] "\t" seqids[i] "\t" $3 "\t" $4 "\t" titles[i]
    }
}

{
    if (NF >= 6) {
        process_with_organisms()
    } else if (NF == 5) {
        process_without_organisms()
    } else {
        print "ERROR: unexpected number of columns on line " NR > "/dev/stderr"
    }
}
' "$blast" > "$outfile"

echo "Complete!"