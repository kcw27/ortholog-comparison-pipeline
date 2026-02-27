library(tidyverse)
library(glue)
library(gridExtra)
library(purrr)
theme_set(theme_classic())

save_sequence_logos <- function(df, group_var = "category", sites, sites_adj, outdir, pdf_suffix, seq_colname) {
  
  if (!group_var %in% names(df)) stop(glue("Column '{group_var}' not found in df"))
  
  df <- df[!is.na(df[[group_var]]), ] # remove rows with NA's for group_var so that the split will work properly
  
  # automatically calculate window_start and window_end to flank the sites closely
  window_start <- max(0, first(sites_adj) - 3)
  window_end <- max(0, last(sites_adj) + 3)
  
  plots <- split(df, df[[group_var]]) |>
    map(~ {
      # Extract sequences for this group
      group_seqs <- .x[[seq_colname]]
      # Subset sequences to the window
      window_seqs <- substring(group_seqs, first = window_start, last = window_end)
      
      ggplot() +
        geom_logo(window_seqs) +
        scale_x_continuous(
          name = "Position",
          breaks = sites_adj - window_start + 1, # to correctly place the label locations 
          labels = sites # label as if looking at the reference genome
        ) +
        labs(title = glue("{group_var}: {unique(.x[[group_var]])[1]}; n: {nrow(.x)}"))
    })

  # also save a table of residue frequencies at each site
  split(df, df[[group_var]]) |>
    map(~ {
      get_site_frequencies(.x, seq_colname, sites, sites_adj, glue("{outdir}/{pdf_suffix}_{unique(.x[[group_var]])[1]}_freqs.tsv"))
      })

  # Split plots into groups of 4
  plot_groups <- split(plots, ceiling(seq_along(plots) / 4))

  # Create a PDF with multiple pages
  pdf(glue("{outdir}/{pdf_suffix}_sequenceLogos.pdf"), width = 12, height = 6)
  print(glue("Sequence logo saved to {outdir}/{pdf_suffix}.pdf"))

  # Loop through groups and print each set of 4 to a new page
  walk(plot_groups, ~ grid.arrange(grobs = .x, nrow = 2, ncol = 2))

  dev.off() # Close the PDF device
}


plot_local_logo <- function(df, outdir, reference_seq, sites_str, seq_colname = "sequence", width = 8, height = 6, group_var = "category") {
  # parse sites_str, a comma-delimited string, into a sorted numeric vector
  sites <- strsplit(sites_str, ",") |> lapply(as.numeric) |> unlist() |> sort()

  # calculate sites_adj based on reference sequence
#  matching_rows <- which(df$sequence_id == reference_label)
#  reference_seq <- df[[seq_colname]][matching_rows[1]] # if there are multiple sequences with the same ID, use the first
  sites_adj <- calibrate_coords(reference_seq, sites)
#  print(reference_label)
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
  
  # make the overall sequence logo
  ggplot() + 
    geom_logo(substring(df[[seq_colname]], first=window_start, last=window_end)) +
    scale_x_continuous(name = "Position",
                       breaks = sites_adj-window_start+1, # to correctly place the label locations 
                       labels = sites) + # label as if looking at the reference genome
                  labs(title = glue("Overall sequence logo; n: {nrow(df)}"))
    # original x axis started at 1
  
  fname = glue("{outdir}/overall_sequenceLogo.pdf")
  ggsave(fname, height=height, width=width, units="in", limitsize = FALSE, create.dir = TRUE)
  print(glue("Saved sequence logo to {fname}!"))
  
  # also save a table of residue frequencies at each site
  get_site_frequencies(df, seq_colname, sites, sites_adj, glue("{outdir}/overall_freqs.tsv"))
  
                                                   
  # Generate sequence logos for group_var
  save_sequence_logos(df, group_var, sites, sites_adj, outdir, glue("{group_var}"), seq_colname)
  
  
  # If group_var is "category" or "genus", can further split the data (by "subcategory" or "species" respectively)
  if (group_var == "category") { # make additional plots for subcategory
    df <- df[!is.na(df$category), ] # remove rows with NA's for group_var so that the split will work properly
    df_categories <- split(df, df$category)
  
    make_dir(glue("{outdir}/subcategories"))
    # Generate subcategory-level sequence logos
    walk(names(df_categories), ~ save_sequence_logos(df_categories[[.x]], "subcategory", sites, sites_adj,
                                                outdir, glue("subcategories/{.x}"), seq_colname
                                                ))
  } else if (group_var == "genus") { # make additional plots for species
    df <- df[!is.na(df$genus), ] # remove rows with NA's for group_var so that the split will work properly
    df_categories <- split(df, df$genus)
  
    make_dir(glue("{outdir}/species"))
    # Generate species-level histogram PDF
    walk(names(df_categories), ~ save_sequence_logos(df_categories[[.x]], "species", sites, sites_adj,
                                                outdir, glue("species/{.x}"), seq_colname
                                                ))
  }
                                                   
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

# creates directories with that name if they don't already exist
make_dir <- function(new_dir) {
  ifelse(!dir.exists(file.path(new_dir)),
        dir.create(file.path(new_dir)),
        glue("{new_dir} directory exists"))
}


get_site_frequencies <- function(df, seq_colname, sites, sites_adj, fname) {
  sites_df <- data.frame(sites, sites_adj)
  
  # Create a list to store results
  all_results <- list()
  
  for (i in 1:nrow(sites_df)) {
    site_original <- sites_df$sites[i]
    site_adjusted <- sites_df$sites_adj[i]
    
    result <- df |>
      # Extract character at the specified position
      mutate(character = substr(!!sym(seq_colname), site_adjusted, site_adjusted)) |>
      # Remove NA values
      filter(!is.na(character)) |>
      # Count frequencies
      dplyr::count(character) |> # if not explicitly specified, uses seqinr::count()
      # Calculate percentages
      mutate(
        frequency_percentage = signif(n / sum(n) * 100, digits = 3),
        site = site_original
      ) |>
      # Select and reorder columns
      select(site, character, frequency_percentage)
    
    all_results[[i]] <- result
  }
  
  # Combine all results and save to fname
  write.table(bind_rows(all_results), file = fname, sep = "\t", row.names = FALSE, quote = FALSE)
  print(glue("Saved frequencies to {fname}!"))
}