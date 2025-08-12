#!/usr/bin/env Rscript

### Inputs:
## args[1], blast_file: A blast output file produced with -outfmt "6 sallgi sallseqid sseq evalue salltitles" and then converted to long format using expand_blast_output.sh
# This may be the output of a protein blast, in which case it likely lacks genome IDs, but will be supplemented by the genome IDs in synteny_hits_file.
## args[2], synteny_hits_file: A synteny summary file produced by find_synteny_hits.sh, which is run by synteny_wrapper.sh
## args[3], outname: name of output file which contains the top BLAST hits per genome ID.

### Script description:
# Inner joins the data in blast_file_ and synteny_hits_file by protein ID, then filters to the top hit (i.e. lowest evalue) per genome ID.

### Example run:
# Rscript $scriptsdir/synteny_search/intersect_blast_and_synteny_with_metadata.R $datadir/PA3565_with_orgs_long.blast $datadir/synteny_outputs/results_65_67/synteny_summary.tsv $datadir/PA3565_67_topHitsPerGenome.blast
# Rscript intersect_blast_and_synteny_with_metadata.R "${HOME}/data/blast_outputs/PA3565_nr_orgs_long.txt" "${HOME}/data/results_65_67/synteny_summary.tsv" "${HOME}/data/blast_outputs/PA3565_67_topHitsPerGenome.tsv"

# for testing:
# Rscript intersect_blast_and_synteny_with_metadata.R "${HOME}/data/testing/out/PA3565_nr_small_orgs_long.txt" "${HOME}/data/results_65_67/synteny_summary.tsv" "${HOME}/data/testing/out/PA3565_67_topHitsPerGenome_small.tsv"


### Import packages
library(tidyverse)
library(glue)


### Define functions
make_dir <- function(new_dir) {
  ifelse(!dir.exists(file.path(new_dir)),
        dir.create(file.path(new_dir)),
        glue("{new_dir} directory exists"))
}

intersect_inputs <- function(blast_df, synteny_df) {
  # Performs an inner join of the blast and synteny hit data by protein ID.
  # Assumes that protein IDs are column 2 in blast_df in a form like the following example: emb|VFT34736.1|
  # Also assumes that protein IDs are column 7 in a form like the following example: VFT34736.1
  
  # Preprocess blast_df protein IDs to be consistent with synteny_df
  blast_df <- blast_df |>
    #mutate(V2 = sub(".*\\|(.*)\\|.*", "\\1", V2))
    mutate(protein_id = sub(".*\\|(.*)\\|.*", "\\1", protein_id))
  
  # now join
  #filtered_df <- inner_join(blast_df, synteny_df, by=join_by(V2 == V7))
  filtered_df <- inner_join(blast_df, synteny_df, by=join_by(protein_id))

  return(filtered_df)
}

take_top_hit_per_genome <- function(filtered_df) {
  # Takes the output of intersect_inputs.
  # Filters it to just the top hit (i.e. lowest evalue, which is column V4.x) per genome
  # Writes it to outname
#  top_df <- filtered_df |>
#    group_by(V1.y) |> # group by genome_id (V1.y)
#    arrange(V4.x) |> # sort by evalue (V4.x) in ascending order (default arrange() behavior)
#    slice_head(n=1) |> # take only the first line per genome ID
#    mutate(V2=paste(V1.y, V2, sep="-")) # diff proteins in diff genomes may have the same protein ID; concatenate genome ID and protein ID for more informative name
    
  top_df <- filtered_df |>
    group_by(genome_id.y) |> # group by genome_id
    arrange(evalue) |> # sort by evalue in ascending order (default arrange() behavior)
    slice_head(n=1) |> # take only the first line per genome ID
    mutate(protein_id=paste(genome_id.y, protein_id, sep="-")) # diff proteins in diff genomes may have the same protein ID; concatenate genome ID and protein ID for more informative name
    
  return(top_df)
}

postprocess <- function(top_df) {
  # remove genome_id.x because it's all 0's, move genome_id.y to the front, and rename to genome_id
  # Also overwrite sequence.x with sequence.y because they aren't necessarily the same; sequence.x is from BLAST and represents similar sequences that are clustered,
  # while sequence.y corresponds to the specific sequence (the BLAST hit with the lowest evalue for that genome_id) that corresponds to the metadata 
  # Additionally, overwrite organism.x with organism.y because organism.y is the specific organism that corresponds to sequence.y.
#  top_df <- top_df |>
#    subset(select = -genome_id.x) |>
#    relocate(genome_id.y) |> # default behavior: move the col to the front
#    rename(genome_id = genome_id.y) |> # new_name = old_name
    
  top_df <- top_df |>
    mutate(genome_id.x = genome_id.y, sequence.x = sequence.y, organism.x = organism.y) |>
    rename(genome_id = genome_id.x, sequence = sequence.x, organism = organism.x) |>
    subset(select = -c(genome_id.y, sequence.y, organism.y))
    # We mutate and rename, then drop, in order to ensure the new data is in the correct column positions to follow the BLAST file format
    
  return(top_df)
}


### Process CLIs
args <- commandArgs(trailingOnly = TRUE) # only get the CLIs that come after the name of the script

if (length(args) < 3) {
  stop("Please provide at least three arguments: <blast_file> <synteny_hits_file> <outname>")
}

blast_file <- args[1]
synteny_hits_file <- args[2]
outname <- args[3]

# print the arguments
cat("Blast file provided:", blast_file, "\n")
cat("Synteny hits file:", synteny_hits_file, "\n")
cat("Output file:", outname, "\n")

# make the outdir for the output file if necessary
print("Making outdir if necessary:")
make_dir(dirname(outname))

# read blast and synteny hit files as dfs
print("Reading blast and synteny hit files as dataframes:")
blast_df <- read.csv(blast_file, header=FALSE, sep="\t")
colnames(blast_df) <- c("genome_id", "protein_id", "sequence", "evalue", "protein_title", "organism")
synteny_df <- read.csv(synteny_hits_file, header=FALSE, sep="\t")
colnames(synteny_df) <- c("genome_id", "contig", "organism", "isolation_source", "titles", "locus", "protein_id", "sequence", "sequencing_technology", "assembly_method")[1:ncol(synteny_df)]

glue("Number of rows in blast_df: {nrow(blast_df)}")
glue("Number of rows in synteny_df: {nrow(synteny_df)}")

#head(blast_df$V2)
#head(synteny_df$V7)
head(blast_df$protein_id)
head(synteny_df$protein_id)

# Do an inner join to keep only information corresponding to protein IDs that intersect between the two datasets
filtered_df <- intersect_inputs(blast_df, synteny_df)
glue("Number of rows in filtered_df: {nrow(filtered_df)}")
#write.table(filtered_df, file="data/testing/out/PA3565_67_filteredIntermediateDF.tsv", sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
colnames(filtered_df)
head(filtered_df, n=1)

# Take the top hit per genome in this intersected df; write to file
top_df <- take_top_hit_per_genome(filtered_df)
glue("Number of rows in top_df: {nrow(top_df)}")

print("Postprocessing...")
top_df <- postprocess(top_df)

# write to tsv file (no row names; I've just modified it to include column names)
write.table(top_df, file=outname, sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)
print(glue("File written to {outname}"))