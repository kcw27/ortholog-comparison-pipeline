library(tidyverse)
library(glue)


plot_local_logo <- function(df, outdir, reference_label, sites_str, seq_colname = "sequence", width = 8, height = 6) {
  # parse sites_str, a comma-delimited string, into a sorted numeric vector
  sites <- strsplit(sites_str, ",") |> lapply(as.numeric) |> unlist() |> sort()

  # calculate sites_adj based on reference sequence
  matching_rows <- which(df$sequence_id == reference_label)
  reference_seq <- df[[seq_colname]][matching_rows[1]] # if there are multiple sequences with the same ID, use the first
  sites_adj <- calibrate_coords(reference_seq, sites)
  print(reference_label)
  print(reference_seq)
  print(sites)
  print(sites_adj)
  
  # automatically calculate window_start and window_end to flank the sites closely
  window_start <- max(0, first(sites_adj) - 3)
  window_end <- max(0, last(sites_adj) + 3)
  
  if (!requireNamespace("ggseqlogo", quietly = TRUE)) {
    install.packages("ggseqlogo")
  }
  library(ggseqlogo)
  
  ggplot() + 
    geom_logo(substring(df[[seq_colname]], first=window_start, last=window_end)) +
    scale_x_continuous(name = "Position",
                       breaks = sites_adj-window_start+1, # to correctly place the label locations 
                       labels = sites) # label as if looking at the reference genome
    # original x axis started at 1
  
  fname = glue("{outdir}/sequencelogo.pdf")
  ggsave(fname, height=height, width=width, units="in", create.dir = TRUE)
}


calibrate_coords <- function(seq, sites_of_interest) {
  # seq is a reference sequence, and sites_of_interest are the sites in the reference seq you want to annotate in the seqlogo
  # in an alignment, gap characters are inserted, so the coordinates need to be calibrated accordingly
  # at this point assume sites_of_interest is a sorted numeric vector
  non_gap = 0
  i = 1 # R is 1-indexed
  last_site <- max(sites_of_interest)
  
  sites_adj <- c()
  
  while (non_gap < last_site & i <= nchar(seq)) {
    if (substr(seq, i, i) != "-") {
      non_gap = non_gap + 1
      
      # Only add to sites_adj when we're at a non-gap position AND it matches a site of interest
      if (non_gap %in% sites_of_interest) {
        sites_adj <- append(sites_adj, i)
      }
    } 
    
    i = i + 1
  }
  
  return(sites_adj)
}