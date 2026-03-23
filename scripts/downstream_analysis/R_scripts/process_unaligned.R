### Inputs:
## args[1], multifasta_file: path to a multifasta (sequences not aligned).
# Assume the sequence names are in the form "genomeID-proteinID".
## args[2], metadata_file: path to a metadata file
## args[3], outdir: path to output directory

### Example runs:
# Rscript /home/kcw2/ortholog-comparison-pipeline/scripts/downstream_analysis/R_scripts/process_unaligned.R "/home/kcw2/data/blast_outputs/pseudomonas_aeruginosa_ONLY_fha1_topPerGenome_completeSequences_concatGenomeProteinIDs.fasta" "/home/kcw2/data/blast_outputs/pa_fha1_top_complete_metadata.blast"
# Rscript /home/kcw2/ortholog-comparison-pipeline/scripts/downstream_analysis/R_scripts/process_unaligned.R "/home/kcw2/data/blast_outputs/pseudomonas_aeruginosa_PA3565_67_synteny_PairwiseBlastIntersected_pident99.fasta" "/home/kcw2/data/results_65_67/synteny_summary.tsv" /home/kcw2/data/testing/foo


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
  stopifnot(is.character(multifasta_file), length(multifasta_file) == 1)
  
  data <- read.alignment(multifasta_file, format="fasta") # can still use read.alignment on non-aligned sequences
  
  seqs <- data[["seq"]] |>
    unlist() |> # convert from list to vector
    toupper() # originally in lowercase; convert to uppercase
  
  sequence_names = (data["nam"][[1]]) # need the [[1]] in order to get the actual list
  df <- data.frame(sequence_id = sequence_names, sequence = seqs) |>
    mutate(sequence_length = nchar(gsub("-", "", sequence))) # exclude gap characters from sequence length
  
  # previously, the only types of sequence id were protein_id and genome_id-protein_id
#  if (length(intersect(metadata$sequence_id, df$sequence_id)) == 0) { # maybe the sequence names are just protein IDs?
#    if (metadata_type == "synteny_summary") {
#      sequence_name_col <- "protein_id"
#    } else if (metadata_type == "fetched") {
#      sequence_name_col <- "subject"
#    }
#  } 
  
  # updated approach accounts for genome_id-protein_id-locus_id sequence IDs. Counts "-" characters.
  # make whatever will serve as sequence_name_col in the metadata file
  if (metadata_type == "synteny_summary") {
    protein_col <- "protein_id"
  } else if (metadata_type == "fetched") {
    protein_col <- "subject"
  }
  
  num_dashes <- str_count(sequence_names[1], "-")
  
  if (num_dashes == 0) { # just the protein_id
    sequence_name_col <- protein_col
  } else if (num_dashes == 1) { # genome_id-protein_id
    metadata <- metadata |>
      mutate(gp_id = paste(genome_id, !!sym(protein_col), sep="-"))
    sequence_name_col <- "gp_id"
  } else { # assume there are two... if there are more, maybe they're part of the locus tag
    metadata <- metadata |>
      mutate(gpl_id = paste(genome_id, !!sym(protein_col), locus_tag, sep="-"))
    sequence_name_col <- "gpl_id"
  }
  print(glue("Sequence name col: {sequence_name_col}"))
  
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
  # Use data.table with maximum error tolerance
  if (!requireNamespace("data.table", quietly = TRUE)) {
    install.packages("data.table")
  }
  library(data.table)
  
  # Read with maximum error tolerance
  metadata <- fread(metadata_file, 
                   sep = "\t",
                   header = TRUE,
                   quote = "",
                   fill = TRUE,
                   skip = 0,
                   nThread = 1,  # Single thread for stability
                   stringsAsFactors = FALSE,
                   showProgress = FALSE)
  
  metadata <- as.data.frame(metadata)
  
  cat("Successfully read", nrow(metadata), "rows (some lines may have been skipped)\n")
  
  # Rest of your processing code remains the same...
#  if (metadata_type == "synteny_summary") {
#    metadata <- metadata |>
#      mutate(sequence_id = paste(genome_id, protein_id, sep = "-"))
#  } else if (metadata_type == "fetched") {
#    metadata <- metadata |>
#      mutate(sequence_id = paste(genome_id, subject, sep = "-"))
#  }
  
  if ("organism" %in% colnames(metadata)) {
    metadata <- metadata |>
      mutate(genus = sapply(organism, function(x) strsplit(x, " ")[[1]][1])) |>
      mutate(genus = ifelse(is.na(genus), " ", genus)) |>
      mutate(species = sapply(organism, function(x) strsplit(x, " ")[[1]][2])) |>
      mutate(species = ifelse(is.na(species), " ", species))
  }
  
  # Replace any NA ids with empty strings
  if ("genome_id" %in% colnames(metadata)) {
    metadata <- metadata |>
      mutate(genome_id = ifelse(is.na(genome_id), "", genome_id))
  }
  if ("locus_tag" %in% colnames(metadata)) {
    metadata <- metadata |>
      mutate(locus_tag = ifelse(is.na(locus_tag), "", locus_tag))
  }
  
  # protein id column is named differently depending on source
   if (metadata_type == "synteny_summary") {
     metadata <- metadata |>
       mutate(protein_id = ifelse(is.na(protein_id), "", protein_id))
   } else if (metadata_type == "fetched") {
     metadata <- metadata |>
       mutate(subject = ifelse(is.na(subject), "", subject))
   }
  
  return(metadata)
}

#get_script_dir <- function() {
#  args <- commandArgs(trailingOnly = FALSE)
#  file_arg <- "--file="
#  script_path <- sub(file_arg, "", args[grep(file_arg, args)])
#  if (length(script_path) == 0) {
#    stop("Cannot determine script path. Are you running via Rscript?")
#  }
#  normalizePath(dirname(script_path))
#}


# Function for Shiny to call
process_unaligned_shiny <- function(multifasta_file, metadata_file, outdir, script_dir, log_fn = print) {
  # Determine metadata type using the Python function
  #script_dir <- get_script_dir()
  metadata_py <- normalizePath(file.path(script_dir, "..", "metadata_processing.py"))
  
  if (!file.exists(metadata_py)) {
    stop(glue("Python script not found at: {metadata_py}"))
  }
  metadata_module <- import_from_path("metadata_processing", path = dirname(metadata_py))
  metadata_type <- metadata_module$determine_origin(metadata_file)
  
  if (!metadata_type %in% c("synteny_summary", "fetched")) {
    stop(glue("Error: metadata_type {metadata_type} not recognized"))
  }

  # Process the data
  metadata <- process_metadata(metadata_file, metadata_type)
  sequence_name_col <- "sequence_id"
  df_raw <- process_data(multifasta_file, metadata, metadata_type, sequence_name_col)
  
  if ("category" %in% names(df_raw)) {
    df <- df_raw |>
      filter(!is.na(category)) |>
      filter(category != "no category")
    
    log_fn(glue("Benchmarking: {nrow(df_raw)} sequences in input BLAST file; {nrow(df)} were successfully categorized."))
    
    log_output <- df_raw |>
      group_by(category) |>
      summarize(n = n()) |>
      capture.output()
    
    lapply(log_output, log_fn)
  } else {
    df <- df_raw
    log_fn(glue("Benchmarking: {nrow(df_raw)} sequences in input BLAST file"))
  }
   
  stats_dir <- glue("{script_dir}/hypothesis_testing")
  script_files <- list.files(stats_dir, pattern = "\\.R$", full.names = TRUE)
  for (f in script_files) source(f)
  
    
  numeric_col = ifelse("sequence_length" %in% colnames(df), "sequence_length", "sequence_length.x")
  
  # Determine which column to use for normality test
  if ("sequence_length" %in% colnames(df)) {
    numeric_col <- "sequence_length"
  } else {
    numeric_cols <- names(df)[sapply(df, is.numeric)]
    if (length(numeric_cols) > 0) {
      numeric_col <- numeric_cols[1]
      log_fn(glue("Using first numeric column '{numeric_col}' for normality test"))
    } else {
      test_type <- "non-parametric"
      log_fn("No numeric columns found - using default test type: non-parametric")
    }
  }
  
  # Only run normality test if we found a numeric column
  if (exists("numeric_col")) {
    # Capture printed output only
    test_output <- capture.output({
      invisible(test_type <- test_normality(df, numeric_col, glue("{outdir}/figures")))
    })
    lapply(test_output, log_fn)
    log_fn(glue("Use this test type on the '{numeric_col}' numerical variable: {test_type}"))
  } else {
    # test_type is already set to "non-parametric" above
    log_fn("Skipping normality test - no suitable numeric column available")
  }

  #log_fn(glue("Use this test type on the 'sequence_length' numerical variable: {test_type}"))

  # Save categorized IDs
  ids_file <- file.path(outdir, "categorized_ids.txt")
  write.table(df$sequence_id, ids_file,
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  log_fn(glue("IDs of categorized sequences saved to {outdir}/categorized_ids.txt"))
  
  #return(list(df = df, test_type = test_type))
  return(list(df = df_raw, test_type = test_type)) # return df_raw instead so that you can customize filtering later
}

# Only run as standalone script if called via Rscript with CLI arguments
if (!interactive() && length(commandArgs(trailingOnly = TRUE)) >= 2) {
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
  
  if ("category" %in% names(df_raw)) {
    df_raw |>
      group_by(category) |>
      summarize(n=n())
  }
  #df |>
  #  group_by(category) |>
  #  summarize(n=n())
  
  
  ### Call statistical tests on df
  # First, source the directory of hypothesis testing scripts
  #stats_dir <- here("scripts/downstream_analysis/R_scripts/hypothesis_testing") # here() finds the project root
  
  stats_dir <- glue("{script_dir}/hypothesis_testing")
  script_files <- list.files(stats_dir, pattern = "\\.R$", full.names = TRUE)
  for (f in script_files) source(f)
  
  #numeric_col = ifelse("sequence_length" %in% colnames(df), "sequence_length", "sequence_length.x")
  #test_type <- test_normality(df, numeric_col, glue("{outdir}/figures"))
  test_type <- test_normality(df, "sequence_length", glue("{outdir}/figures"))
  print(glue("Use this test type on the '{numeric_col}' numerical variable: {test_type}"))
  
  df |>
    pull(sequence_id) |>
    write.table(glue("{outdir}/categorized_ids.txt"),
      row.names = FALSE, col.names = FALSE, quote = FALSE)
  print(glue("IDs of categorized sequences saved to {outdir}/categorized_ids.txt"))
}