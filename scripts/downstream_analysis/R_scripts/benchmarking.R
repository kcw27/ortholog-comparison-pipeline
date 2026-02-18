# Input: args[1], path to metadata file that has been acategorized
# Output: prints benchmarking

library(dplyr)
library(glue)

benchmark <- function(metadata_file) {
  #df <- read.csv(metadata_file, header=TRUE, sep="\t")
  if (!requireNamespace("data.table", quietly = TRUE)) {
    install.packages("data.table")
  }
  library(data.table)
  
  # Read with maximum error tolerance
  df <- fread(metadata_file, 
                   sep = "\t",
                   header = TRUE,
                   quote = "",
                   fill = TRUE,
                   skip = 0,
                   nThread = 1,  # Single thread for stability
                   stringsAsFactors = FALSE,
                   showProgress = FALSE)
  
  df <- as.data.frame(df)
  
  cat("Successfully read", nrow(df), "rows (some lines may have been skipped)\n")
  
  print(glue("There are {nrow(df)} lines in the metadata file."))
  
  n_isosource <- df |>
    filter(!isolation_source %in% c("", "NA")) |>
    nrow()
    
  print(glue("{n_isosource} rows were associated with isolation source metadata."))
  
  if ("rescued_source" %in% names(df)) {
    n_rescue <- df |> # how many didn't have isosource but did have rescue?
      filter(isolation_source %in% c("", "NA")) |>
      filter(rescued_source != "") |>
      nrow()
      
    print(glue("Among the {nrow(df) - n_isosource} rows without isolation source, {n_rescue} isolation sources were rescued from the titles or isolation_site column, making {n_isosource + n_rescue} rows eligible for categorization."))
  }
  
  if ("category" %in% names(df) & "subcategory" %in% names(df)) {
    n_cat <- df |>
      filter(!category %in% c("", "no category")) |>
      nrow()
    n_subcat <- df |>
      filter(!subcategory %in% c("", "no subcategory")) |>
      nrow()
      
    print(glue("{n_cat} rows were successfully assigned a category, and {n_subcat} rows were successfully assigned a subcategory."))
  }
}



## Only run as standalone script if called via Rscript with CLI arguments
#if (identical(environmentName(globalenv()), "R_GlobalEnv") && !interactive()) {
#  args <- commandArgs(trailingOnly = TRUE)
#  if (length(args) < 1) {
#    stop("Please provide an argument: <metadata_file> ")
#  }
#  metadata_file <- args[1]
#  benchmark(metadata_file)
#}
