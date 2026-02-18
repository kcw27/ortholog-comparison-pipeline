# Script that takes the direct output of run_local_blast.sh (or a file with the same format)
# and produces figures to summarize the BLAST outputs.
# Note: if the optional outdir argument is provided, it will make that directory and save to a
# "blast_summary" subdirectory.

# This is intended for preliminary, exploratory visualization. To make figures with better-scaled axes,
# please manually generate them.

### Import packages
library(tidyverse)
library(glue)

theme_set(theme_classic())

### Define functions
make_dir <- function(new_dir) {
  ifelse(!dir.exists(file.path(new_dir)),
        dir.create(file.path(new_dir)),
        glue("{new_dir} directory exists"))
}

# wrapper to call all of the functions
# If you only want to use a single one, soource the script as a module in the command line
# https://linuxvox.com/blog/call-r-function-in-linux-command-line/#sourcing-external-r-functions
blast_summary <- function(blast_file, outdir) {
  # Load blast as df, providing column names
  df <- read.csv(blast_file, header=FALSE, sep="\t")
  colnames(df) <- c("genome_id", "protein_id", "sequence", "evalue", "protein_title", "organism")
  
  blast_benchmarking(df)
  genomeid_hist(df, outdir)
  seqlen_hist(df, outdir)
  evalue_hist(df, outdir)
  evalue_seqlen_plot(df, outdir)
}

# Benchmarking: print statements
blast_benchmarking <- function(df) {
  print(glue("Number of BLAST hits: ", nrow(df) ))
  print(glue("Number of unique genome IDs: ", length(unique(df$genome_id)) ))
  print(glue("Number of unique sequence IDs: ", length(unique(df$protein_id)) ))
  print(glue("Number of unique sequences: ", length(unique(df$sequence)) ))
  print(glue("Number of unique titles: ", length(unique(df$protein_title)) ))
  print(glue("Number of unique organisms: ", length(unique(df$organism)) ))
}


# Histogram: How many hits per genome are there?
genomeid_hist <- function(df, outdir) {
  fname <- file.path(outdir, "genomeid_hist.pdf")
  
  summary_df <- df |>
    group_by(genome_id) |>
    summarize(n=n()) |>
    arrange(desc(n))

  print(summary_df)

  summary_df |> ggplot() +
    geom_histogram(aes(x=n), binwidth=1) +
    labs(y="Count", x="Number of BLAST hits per genome ID", title=glue("Genome ID frequencies (n={nrow(df)} BLAST hits)"))
  
  ggsave(fname)
}


# Histogram: distribution of sequence lengths
seqlen_hist <- function(df, outdir) {
  fname <- file.path(outdir, "seqlen_hist.pdf")
  
  df |> 
    mutate(sequence_length = nchar(sequence)) |> 
    ggplot() +
    geom_histogram(aes(x=sequence_length)) +
    labs(y="Count", x="Sequence length (# amino acids)", title=glue("Sequence length frequencies (n={nrow(df)} BLAST hits)"))
  
  ggsave(fname)
}

# Histogram: Distribution of evalues
evalue_hist <- function(df, outdir) {
  fname <- file.path(outdir, "evalue_hist.pdf")
  
  df |>
    ggplot() +
    geom_histogram(aes(x=evalue)) +
    labs(y="Count", x="Evalue", title=glue("Evalue frequencies (n={nrow(df)} BLAST hits)")) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) # scales package, part of ggplot2 https://stackoverflow.com/questions/49723196/is-scales-package-included-in-tidyverse
  
  ggsave(fname)
}

# Scatter plot of evalue (y) to sequence length (x)
evalue_seqlen_plot <- function(df, outdir) {
  fname <- file.path(outdir, "evalue_seqlen_plot.pdf")
  
  df |> 
    mutate(sequence_length = nchar(sequence)) |> 
    ggplot() +
    geom_point(aes(x=sequence_length, y=evalue), alpha=0.4) +
    labs(y="Evalue", x="Sequence length (# amino acids)", title=glue("Sequence length vs evalue plot (n={nrow(df)} BLAST hits)"))
  
  ggsave(fname)
}


### Process CLIs
args <- commandArgs(trailingOnly = TRUE) # only get the CLIs that come after the name of the script

if (length(args) < 1) {
  stop('Please provide the path to a BLAST file. (Optional: outdir, if not saving in same location as BLAST file.)')
}

blast_file <- args[1]
outdir <- args[2]

# if no outdir provided, save to a subdirectory of the directory the BLAST file is in
if (is.na(outdir)) {
  outdir <- file.path(dirname(blast_file), "blast_summary")
}

make_dir(outdir) # ensure outdir exists

# get summary outputs now
blast_summary(blast_file, outdir)