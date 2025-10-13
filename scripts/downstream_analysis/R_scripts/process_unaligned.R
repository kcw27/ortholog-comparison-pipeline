### Inputs:
## args[1], multifasta_file: path to a multifasta (sequences not aligned).
# Assume the sequence names are in the form "genomeID-proteinID".
## args[2], metadata_file: path to a metadata file
## args[3], metadata_type: "synteny_summary" or "fetched"
## 

### Example runs:
# Rscript process_unaligned.R "/home/kcw2/data/blast_outputs/pseudomonas_aeruginosa_ONLY_fha1_topPerGenome_completeSequences_concatGenomeProteinIDs.fasta" "/home/kcw2/data/blast_outputs/pa_fha1_top_complete_metadata.blast" "fetched"
# include one for the synteny_summary metadata too

### Import libraries
library(tidyverse) # working with dataframes
library(seqinr) # to read alignments
library(glue) # fstrings

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
  df <- data.frame(locus_tags = sequence_names, sequence = seqs) |>
    mutate(sequence_length = nchar(gsub("-", "", sequence))) # exclude gap characters from sequence length
  
  if (length(intersect(metadata$locus_tags, df$locus_tags)) == 0) { # maybe the sequence names are just protein IDs?
    if (metadata_type == "synteny_summary") {
      sequence_name_col <- "protein_id"
    } else if (metadata_type == "fetched") {
      sequence_name_col <- "subject"
    }
  } 
  
  df <- left_join(df, metadata, by=c("locus_tags"=sequence_name_col))
  
  # add in short_names if name_map_file was provided
  if (!is.null(name_map_file)) {
    name_map <- read.csv(name_map_file, header=FALSE, sep="\t", col.names=c("short_names", "locus_tags")) 
    # as created by fasta_to_phylip.sh: two tab-separated columns, no header, short names first
    df <- left_join(df, name_map, by="locus_tags") |>
      relocate(short_names, .after = locus_tags) # move short_names so that it's column 2, after locus_tags
  }
  
  return(df)
}

process_metadata <- function(metadata_file, metadata_type) {
  # loads metadata based on where it came from
  if (metadata_type == "synteny_summary") {
    metadata <- read.csv(metadata_file, header=FALSE, sep="\t")
    
    colnames(synteny_df) <- c("genome_id", "contig", "organism", 
      "isolation_source", "titles", "locus", "protein_id", "sequence", 
      "sequencing_technology", "assembly_method")[1:ncol(metadata)]
      
    metadata <- metadata |>
      mutate(locus_tags = paste(genome_id, protein_id, sep = "-"))
      
  } else if (metadata_type == "fetched") {
    metadata <- read.csv(metadata_file, header=TRUE, sep="\t")
    
    metadata <- metadata |>
      mutate(locus_tags = paste(genome_id, subject, sep = "-"))
  } else {
    print(glue("Error: invalid metadata_type {metadata_type}; must be 'synteny_summary' or 'fetched'."))
  }
  
  return(metadata)
}

### Take CLIs
args <- commandArgs(trailingOnly = TRUE) # only get the CLIs that come after the name of the script

if (length(args) < 3) {
  stop("Please provide at least three arguments: <multifasta_file> <metadata_file> <metadata_type>")
}

multifasta_file <- args[1]
metadata_file <- args[2]
metadata_type <- args[3]

# print the arguments
cat("Multifasta file provided:", multifasta_file, "\n")
cat("Metadata file:", metadata_file, "\n")
cat("Metadata type:", metadata_type, "\n")

### Process data
metadata <- process_metadata(metadata_file, metadata_type)

sequence_name_col <- "locus_tags"
df_raw <- process_data(multifasta_file, metadata, metadata_type, sequence_name_col)
df <- df_raw |>
  filter(!is.na(category)) |>
  filter(category != "no category")

print(glue("Benchmarking: {nrow(df_raw)} sequences in input BLAST file; {nrow(df)} were successfully categorized."))

df_raw |>
  group_by(category) |>
  summarize(n=n())
  
df |>
  group_by(category) |>
  summarize(n=n())


### Call statistical tests on df
# First, source the directory of hypothesis testing scripts
library(here) # finds the project root
stats_dir <- here("scripts/downstream_analysis/R_scripts/hypothesis_testing")

script_files <- list.files(stats_dir, pattern = "\\.R$", full.names = TRUE)
for (f in script_files) source(f)

test_type <- test_normality(df, "sequence_length")
print(glue("Use this test type: {test_type}"))