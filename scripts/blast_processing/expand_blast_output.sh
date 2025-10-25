# Description:
# (long format- one record per row)

### Example runs:
# with organisms column, default outname:
# bash "/home/kcw2/ortholog-comparison-pipeline/scripts/blast_processing/expand_blast_output.sh" -b "${HOME}/data/testing/out/PA3565_nr_small_orgs.txt"

# without organisms column, with a specified outname: 
# bash "/home/kcw2/ortholog-comparison-pipeline/scripts/blast_processing/expand_blast_output.sh" -b "/home/kcw2/data/testing/PA3565_nr_small.txt" -o "/home/kcw2/data/testing/PA3565_nr_small_foo.txt"

# Parse flags
while [[ $# -gt 0 ]]; do
  case "$1" in
    -b|--blast)
      blast="$2"
      shift 2
      ;;
    -o|--outfile)
      output="$2"
      shift 2
      ;;
    -h|--help)
      echo "Usage: $0 -b <blast_file> [-o <output_file>]"
      echo "  -b, --blast     Mandatory input BLAST file."
      echo "Expected to have 6 columns, e.g. -outfmt 6 with arguments sallgi sallseqid sseq evalue salltitles."
      echo "May have organism column from add_organism_column.sh"
      echo "  -o, --outfile   Optional output file (default: <blast>_long.blast)"
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      exit 1
      ;;
  esac
done

# Check mandatory argument
if [[ -z "$blast" ]]; then
  echo "Error: -b|--blast argument is required." >&2
  exit 1
fi

# Assign default outfile only if not set
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