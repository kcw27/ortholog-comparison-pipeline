### Inputs:
## args[1], multifasta_file: path to a multifasta (sequences not aligned).
# Assume the sequence names are in the form "genomeID-proteinID".
## args[2], metadata_file: path to a metadata file
## args[3], outdir: path to output directory

### Example runs:
# Rscript /home/kcw2/ortholog-comparison-pipeline/scripts/downstream_analysis/R_scripts/process_unaligned.R "/home/kcw2/data/blast_outputs/pseudomonas_aeruginosa_ONLY_fha1_topPerGenome_completeSequences_concatGenomeProteinIDs.fasta" "/home/kcw2/data/blast_outputs/pa_fha1_top_complete_metadata.blast"
# Rscript /home/kcw2/ortholog-comparison-pipeline/scripts/downstream_analysis/R_scripts/process_unaligned.R "/home/kcw2/data/blast_outputs/pseudomonas_aeruginosa_PA3565_67_synteny_PairwiseBlastIntersected_pident99.fasta" "/home/kcw2/data/results_65_67/synteny_summary.tsv" /home/kcw2/ortholog-comparison-pipeline/scripts/downstream_analysis/R_scripts/foo


### Import libraries
library(tidyverse) # working with dataframes
library(seqinr) # to read alignments
library(glue) # fstrings
library(reticulate) # call Python functions
#library(here) # finds the project root

### Define functions
process_data <- function(multifasta_file, metadata, metadata_type, sequence_name_col, name_map_file=NULL) {
  # Takes as input the filename of a multifasta
  # as well as a metadata file (with whatever preprocessing is necessary)
  # sequence_name_col is whichever column in the metadata serves as the long sequence names (e.g. protein_id, nucleotide_id)
  # repeats_file has three columns: locus tag, start of repeat, end of repeat
  # Optional argument name_map_file can be supplied to associate long sequence names from this multifasta
  # to the short sequence names from a PHYLIP alignment produced using alignment_and_tree_wrapper.sh.
  
  data <- read.alignment(multifasta_file, format="fasta") # can still use read.alignment on non-aligned sequences
  
  seqs <- data[["seq"]] |>
    unlist() |> # convert from list to vector
    toupper() # originally in lowercase; convert to uppercase
  
  sequence_names = (data["nam"][[1]]) # need the [[1]] in order to get the actual list
  df <- data.frame(sequence_id = sequence_names, sequence = seqs) |>
    mutate(sequence_length = nchar(gsub("-", "", sequence))) # exclude gap characters from sequence length
  
  if (length(intersect(metadata$sequence_id, df$sequence_id)) == 0) { # maybe the sequence names are just protein IDs?
    if (metadata_type == "synteny_summary") {
      sequence_name_col <- "protein_id"
    } else if (metadata_type == "fetched") {
      sequence_name_col <- "subject"
    }
  } 
  
  # if the metadata has multiple records under the same sequence_name_col, take only the first so that we don't duplicate sequences in the left join.
  # Filter to keep the first row for each unique value in the 'group' column
  metadata <- metadata |>
    group_by(.data[[sequence_name_col]]) |> # since sequence_name_col is a string
    slice(1) |>
    ungroup()

  
  df <- left_join(df, metadata, by=c("sequence_id"=sequence_name_col))
  
#  # add in short_names if name_map_file was provided
#  if (!is.null(name_map_file)) {
#    name_map <- read.csv(name_map_file, header=FALSE, sep="\t", col.names=c("short_names", "sequence_id")) 
#    # as created by fasta_to_phylip.sh: two tab-separated columns, no header, short names first
#    df <- left_join(df, name_map, by="sequence_id") |>
#      relocate(short_names, .after = sequence_id) # move short_names so that it's column 2, after sequence_id
#  }
  
  return(df)
}

process_metadata <- function(metadata_file, metadata_type) {
  # loads metadata based on where it came from
  if (metadata_type == "synteny_summary") {
#    metadata <- read.csv(metadata_file, header=FALSE, sep="\t")
    
#    colnames(synteny_df) <- c("genome_id", "contig", "organism", 
#      "isolation_source", "titles", "locus", "protein_id", "sequence", 
#      "sequencing_technology", "assembly_method")[1:ncol(metadata)]

    metadata <- read.csv(metadata_file, header=TRUE, sep="\t")
      
    metadata <- metadata |>
      mutate(sequence_id = paste(genome_id, protein_id, sep = "-"))
      
  } else if (metadata_type == "fetched") {
    metadata <- read.csv(metadata_file, header=TRUE, sep="\t")
    
    metadata <- metadata |>
      mutate(sequence_id = paste(genome_id, subject, sep = "-"))
  } else {
    print(glue("Error: invalid metadata_type {metadata_type}; must be 'synteny_summary' or 'fetched'."))
  }
  
  return(metadata)
}

get_script_dir <- function() {
  # Get command line arguments
  args <- commandArgs(trailingOnly = FALSE)
  
  # Find the --file argument
  file_arg <- grep("^--file=", args, value = TRUE)
  
  # Extract the path and normalize
  if (length(file_arg) > 0) {
    script_path <- normalizePath(sub("^--file=", "", file_arg))
    return(dirname(script_path))
  } else {
    stop("Cannot determine script directory: not run via Rscript or missing --file argument.")
  }
}


### Take CLIs
args <- commandArgs(trailingOnly = TRUE) # only get the CLIs that come after the name of the script

if (length(args) < 2) {
  stop("Please provide at least two arguments: <multifasta_file> <metadata_file>. <outdir> is an optional third argument.")
}

script_dir <- get_script_dir()
#cat("Script is located in:", script_dir, "\n")
#setwd(script_dir)

multifasta_file <- args[1]
metadata_file <- args[2]
outdir <- if (length(args) >= 3) args[3] else script_dir
print(glue("Outdir: {outdir}"))
# Ensure the directory exists
dir.create(dirname(outdir), recursive = TRUE, showWarnings = FALSE)

# guess the metadata type
# Construct absolute path to metadata_processing.py
metadata_py <- normalizePath(file.path(script_dir, "..", "metadata_processing.py"))

# Confirm the file exists
if (!file.exists(metadata_py)) {
  stop(glue("Python script not found at: {metadata_py}"))
}

# Source the Python script using absolute path
metadata_module <- import_from_path("metadata_processing", path = dirname(metadata_py))
metadata_type <- metadata_module$determine_origin(metadata_file)

if (!metadata_type %in% c("synteny_summary", "fetched")) {
  stop(glue("Error: metadata_type {metadata_type} not recognized"))
}

# print the arguments
cat("Multifasta file provided:", multifasta_file, "\n")
cat("Metadata file:", metadata_file, "\n")
cat("Metadata type:", metadata_type, "\n")

### Process data
metadata <- process_metadata(metadata_file, metadata_type)

sequence_name_col <- "sequence_id"
df_raw <- process_data(multifasta_file, metadata, metadata_type, sequence_name_col)
df <- df_raw |>
  filter(!is.na(category)) |>
  filter(category != "no category")

print(glue("Benchmarking: {nrow(df_raw)} sequences in input BLAST file; {nrow(df)} were successfully categorized."))

df_raw |>
  group_by(category) |>
  summarize(n=n())
  
#df |>
#  group_by(category) |>
#  summarize(n=n())


### Call statistical tests on df
# First, source the directory of hypothesis testing scripts
#stats_dir <- here("scripts/downstream_analysis/R_scripts/hypothesis_testing") # here() finds the project root

stats_dir <- glue("{script_dir}/hypothesis_testing")
script_files <- list.files(stats_dir, pattern = "\\.R$", full.names = TRUE)
for (f in script_files) source(f)

test_type <- test_normality(df, "sequence_length", glue("{outdir}/figures"))
print(glue("Use this test type on the 'sequence_length' numerical variable: {test_type}"))

df |>
  pull(sequence_id) |>
  write.table(glue("{outdir}/categorized_ids.txt"),
    row.names = FALSE, col.names = FALSE, quote = FALSE)
print(glue("IDs of categorized sequences saved to {outdir}/categorized_ids.txt"))    
