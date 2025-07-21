# CLIs: $1 is the file to filter, $2 is the output filename
# example run:
# bash get_blast_top_hits.sh ~/data/PA3565_orthologs.txt ~/data/PA3565_orthologs_top.txt
# for testing:
# bash get_blast_top_hits.sh  ~/data/results_65_67/PA3565_orthologs_synteny_filtered.tsv ~/data/test_top.tsv


# Assumption: BLAST output file has -outfmt "6 sgi sseqid sseq evalue stitle"
# Precondition: file is sorted in ascending numerical order by evalue (column 4)

# Look for an organism name in column 5 (i.e. text between square brackets []),
# write it to a new column, and if it's the first time seeing that organism,
# write it to the output file. Since the file is sorted from smallest to biggest evalue,
# the first hit for an organism is assumed the most significant.
# (Does not account for ties.)

# Attempt 1:
#awk -F'\t' '
#{
  #match($5, /\[([^\]]+)\]/, arr);
  #new_col = arr[1] ? arr[1] : "NA";
  #if (!seen[new_col]++) print $0 "\t" new_col;
#}
#' "$1" > "$2"

# Attempt 2:
# Step to sort by organism name and then by evalue (in ascending order)
#awk -F'\t' '{cmd="echo \"" $5 "\" | grep -oP \"(?<=\\[).+(?=\\])\""; cmd | getline new_col; if (new_col == "") new_col="NA"; print $0 "\t" new_col}' $1 | sort -t$'\t' -k6,6 -k4,4g > ~/temp.txt

# Step to keep only the first row seen (i.e. the record with the lowest evalue) for each organism (which is column 6)
#sort -u -t$'\t' -k6,6 ~/temp.txt > $2

# Remove the temp file
#rm ~/temp.txt

# Attempt 3:
# Optimized- cut down on overhead by not using grep within awk.
# The difference between this and the first attempt at this script is that awk is only used to add the organism column.
# Sort commands have been added. The first is to ensure that the file is sorted by organism and then by evalue (because the input file may not be in order),
# the second is to get the top hit for each organism.
awk -F'\t' '
{
    match($5, /\[(.+?)\]/, arr);
    new_col = arr[1] ? arr[1] : "NA";
    print $0 "\t" new_col;
}' "$1" | sort -t$'\t' -k6,6 -k4,4g | sort -u -t$'\t' -k6,6 > "$2"
