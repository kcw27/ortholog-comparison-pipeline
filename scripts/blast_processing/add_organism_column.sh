### Description:
# Takes a protein BLAST output with -outfmt "6 sallgi sallseqid sseq evalue salltitles",
# and using the organism names in salltitles, adds a sixth column featuring organism names,
# which are separated by semicolons if there are multiple records clustered into one line.

### CLIs:
# $1: $blast, a protein BLAST output with -outfmt "6 sallgi sallseqid sseq evalue salltitles"
# $2: $outfile, the filename to which the output should be written

# Example run:
# bash add_organism_column.sh "${HOME}/data/testing/PA3565_nr_small.txt" "${HOME}/data/testing/out/PA3565_nr_small_orgs.txt"
# bash add_organism_column.sh "${HOME}/data/blast_outputs/PA3565_nr.txt" "${HOME}/data/testing/out/PA3565_nr_orgs.txt"
# bash add_organism_column.sh "${HOME}/data/blast_outputs/fha1_nr.txt" "${HOME}/data/testing/out/fha1_nr_orgs.txt"

# bash add_organism_column.sh "${HOME}/data/blast_outputs/PA3565_nr.txt" "${HOME}/data/blast_outputs/PA3565_nr_orgs.txt"
# bash add_organism_column.sh "${HOME}/data/blast_outputs/fha1_nr.txt" "${HOME}/data/blast_outputs/fha1_nr_orgs.txt"

blast=$1
outfile=$2

mkdir -p $(dirname $outfile) # make the outdir if it doesn't exist already
#> $outfile # create the output file

## awk and grep solution
#awk -F'\t' '{ cmd = "echo \"" $5 "\" | grep -oP \"(?<=\\[).*?(?=\\])\" | paste -sd \";\" -"; cmd | getline extracted; 
#  close(cmd); 
#  print $0 "\t" extracted }' "$blast" > "$outfile"
  
## pure awk solution (much faster than the awk + grep solution)
#awk -F'\t' '{
#    content = $5;
#    output = "";
#    while (match(content, /\[[^][]+\]/)) { # match everything inside square brackets
#        val = substr(content, RSTART + 1, RLENGTH - 2); # takes off the square brackets, giving us just the organism names
#        output = (output == "") ? val : output ";" val; # joins the organism names with semicolons
#        content = substr(content, RSTART + RLENGTH); # updates the content variable to advance it past the part that was already processed
#    }
#    print $0 "\t" output; # adds the organism names as the last column, column 6
#}' "$blast" > "$outfile"

# it turns out that there may be multiple pairs of square brackets within a title,
# so return only the first match found within a pair of square brackets for each title (titles are delimited by "<>")
awk -F'\t' '{
    split($5, titles, /<>/)
    output = ""

    for (i = 1; i <= length(titles); i++) {
        title = titles[i]
        if (match(title, /\[[^][]+\]/)) {
            val = substr(title, RSTART + 1, RLENGTH - 2)
            output = (output == "") ? val : output ";" val
        }
    }

    print $0 "\t" output
}' "$blast" > "$outfile"



#paste <(cut -f 5 $blast) | while read title; do
#  echo "Hello"
#  echo $title
#  #org_list=$(tr '<>' '\n' $title | grep -oP '(?<=\[).+(?=\])' | tr '\n' ';')
#  org_list=$(echo $title | grep -oP '(?<=\[).*?(?=\])' | paste -sd ';' -) # .*? for non-greedy matching of characters between brackets; replace newlines with ';' delimiters
#  
#
#  echo $org_list
#done




#IFS=$'\n' # $blast is newline-separated.
# Need to set up the loop this way, otherwise it only reads the first line
#for next in $(cat $blast); do
  #echo "Processing: ${next}"
  
  # parts of the line that may have multiple items
  #gis=$(cut -f 1 $next)
  #seqids=$(cut -f 2 $next)
  #titles=$(cut -f 5 $next)
  
  # parts of the line that are guaranteed to have only one item
  #seq=$(cut -f 3 $next)
  #evalue=$(cut -f 3 $next)
  
  #paste <(echo "$titles" | tr '<>' '\n') | while read title; do
    #echo $title
  #done
#done
