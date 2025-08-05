# Description:
# (long format- one record per row)

# CLI:
# $1: Condensed blast file (expected to have 6 columns, e.g. -outfmt "6 sallgi sallseqid sseq evalue salltitles" and then run add_organism_column.sh)

# Example run:
# bash expand_blast_output.sh "${HOME}/data/testing/out/PA3565_nr_small_orgs.txt"
# bash expand_blast_output.sh "${HOME}/data/testing/out/PA3565_nr_orgs.txt"
# bash expand_blast_output.sh "${HOME}/data/testing/out/fha1_nr_orgs.txt"

# bash expand_blast_output.sh "${HOME}/data/blast_outputs/PA3565_nr_orgs.txt"
# bash expand_blast_output.sh "${HOME}/data/blast_outputs/fha1_nr_orgs.txt"

blast=$1
outfile="${blast/.*}_long.blast"

#> $outfile # create the output file

awk -F'\t' '{
    split($1, gis, /;/)
    split($2, seqids, /;/)
    split($5, titles, /<>/)
    split($6, organisms, /;/)
    n = length(gis)
    # columns 3 (seq) and 4 (evalue) will always have only one value per row.
    
    # make sure that all of these really have the same length
    if (length(gis) != length(seqids) || length(gis) != length(titles) || length(gis) != length(organisms)) {
#      print "gis length: " length(gis)
#      print "seqids length: " length(seqids)
#      print "titles length: " length(titles)
#      print "organisms length: " length(organisms)

      print "ERROR: inconsistent field counts on line " NR > "/dev/stderr"
      #print i ": " gis[i], seqids[i], titles[i], organisms[i] > "/dev/stderr"
      next
    }


    for (i = 1; i <= n; i++) {
        print gis[i] "\t" seqids[i] "\t" $3 "\t" $4 "\t" titles[i] "\t" organisms[i]
    }
}' "$blast" > "$outfile"



#IFS=$'\n' # $blast is newline-separated.
## Need to set up the loop this way, otherwise it only reads the first line
#for next in $(cat $blast); do
#  echo "Processing: ${next}"
#  
#  # parts of the line that may have multiple items
#  gis=$(cut -f 1 $next)
#  seqids=$(cut -f 2 $next)
#  titles=$(cut -f 5 $next)
#  organisms=$(cut -f 6 $next)
#  
#  # parts of the line that are guaranteed to have only one item
#  seq=$(cut -f 3 $next)
#  evalue=$(cut -f 4 $next)
#  
#  paste <(echo "$titles" | tr "<>" '\n') | while read title; do
#    echo $title
#  done
#done