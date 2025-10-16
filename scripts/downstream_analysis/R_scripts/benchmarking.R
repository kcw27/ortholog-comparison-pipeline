# Input: args[1], path to metadata file that has been acategorized
# Output: prints benchmarking

library(dplyr)
library(glue)

benchmark <- function(metadata_file) {
  df <- read.csv(metadata_file, header=TRUE, sep="\t")
  
  print(glue("There are {nrow(df)} lines in the metadata file."))
  
  n_isosource <- df |>
    filter(!isolation_source %in% c("", "NA")) |>
    nrow()
    
  print(glue("{n_isosource} rows were associated with isolation source metadata."))
  
  n_rescue <- df |> # how many didn't have isosource but did have rescue?
    filter(isolation_source %in% c("", "NA")) |>
    filter(rescued_source != "") |>
    nrow()
    
  print(glue("Among the {nrow(df) - n_isosource} rows without isolation source, {n_rescue} isolation sources were rescued from the titles column, making {n_isosource + n_rescue} rows eligible for categorization."))
  
  n_cat <- df |>
    filter(!category %in% c("", "no category")) |>
    nrow()
  n_subcat <- df |>
    filter(!subcategory %in% c("", "no subcategory")) |>
    nrow()
    
  print(glue("{n_cat} rows were successfully assigned a category, and {n_subcat} rows were successfully assigned a subcategory."))
}


### Take CLIs
args <- commandArgs(trailingOnly = TRUE) # only get the CLIs that come after the name of the script

if (length(args) < 1) {
  stop("Please provide an argument: <metadata_file> ")
}

metadata_file <- args[1]
benchmark(metadata_file)
