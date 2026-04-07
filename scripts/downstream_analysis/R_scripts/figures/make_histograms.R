library(tidyverse)
library(glue)
library(gridExtra)
library(purrr)
theme_set(theme_bw())

# Function to save histograms based on grouping variable
save_histograms <- function(df, group_var, outdir, pdf_suffix, column_name = "sequence_length", reference_label=NA, reference_length=NA, consistent_y_axis = TRUE, export_locus_tags = FALSE, width = 10) {
  if (!group_var %in% names(df)) stop(glue("Column '{group_var}' not found in df"))
  if (!column_name %in% names(df)) stop(glue("Column '{column_name}' not found in df"))
  
  df <- df[!is.na(df[[group_var]]), ] # remove rows with NA's for group_var so that the split will work properly

  ## To create a set of plots with consistent x and y axes, calculate the x and y limits.
  # x limits for a histogram are just the range of values
  x_limits = range(df[[column_name]])
  buffer_size = 0.10 # make buffers that are some percentage of the range
  x_limits <- x_limits + c(-buffer_size*diff(x_limits), buffer_size*diff(x_limits)) 
  # add a bit of buffer on either side so that it doesn't crop out the data at either end of the x axis
  
  if (!is.na(reference_label) & !is.na(reference_length)) {
  
    plots <- split(df, df[[group_var]]) |>
      map(~ ggplot(.x) +
            geom_histogram(aes(x = .data[[column_name]]), binwidth = width) +
            geom_vline(aes(xintercept=reference_length, colour="reference_length"), size=0.5) +
            geom_vline(aes(xintercept=mean(.data[[column_name]], na.rm = TRUE), colour="mean"), size=0.5) +
            geom_vline(aes(xintercept=median(.data[[column_name]], na.rm = TRUE), colour="median"), size=0.5) +
            labs(x = column_name, y = "Count", title = glue("{group_var}: {unique(.x[[group_var]])[1]}; n: {nrow(.x)}")) +
            scale_color_manual(name = "Values",
                               values = c(reference_length = "blue", mean = "red", median = "green"),
                               labels = c(reference_length = glue("{reference_label} = {reference_length}"),
                                          mean = glue("Mean = {round(mean(.x[[column_name]], na.rm = TRUE), 1)}"),
                                          median = glue("Median = {median(.x[[column_name]], na.rm = TRUE)}"))) +
            theme(legend.text = element_text(size = 8)) + # in case the reference_label is long
            xlim(x_limits) # notice that we set x limits here, but no y limits yet
          )
  } else {
    plots <- split(df, df[[group_var]]) |>
      map(~ ggplot(.x) +
            geom_histogram(aes(x = .data[[column_name]]), binwidth = width) +
            geom_vline(aes(xintercept=mean(.data[[column_name]], na.rm = TRUE), colour="mean"), size=0.5) +
            geom_vline(aes(xintercept=median(.data[[column_name]], na.rm = TRUE), colour="median"), size=0.5) +
            labs(x = column_name, y = "Count", title = glue("{group_var}: {unique(.x[[group_var]])[1]}; n: {nrow(.x)}")) +
            scale_color_manual(name = "Values",
                               values = c(mean = "red", median = "green"),
                               labels = c(mean = glue("Mean = {round(mean(.x[[column_name]], na.rm = TRUE), 1)}"),
                                          median = glue("Median = {median(.x[[column_name]], na.rm = TRUE)}"))) +
            xlim(x_limits) # notice that we set x limits here, but no y limits yet
          )
  }

  if (consistent_y_axis) {
    # Extract max y value from each plot using ggplot_build()
    y_max <- plots |>
      map(~ ggplot_build(.x)$data[[1]]$count) |> # what are the counts for each bin in the histograms?
      flatten_dbl() |> # flatten list of lists into a simple vector
      max() # the largest count among the histograms will be used as the upper y limit
  
    # Retroactively apply y limits to each of the plots
    # Since plots is a vector, need to use map() to add the y limits to each plot in the vector
    plots <- map(plots, ~ .x + ylim(0, y_max))
  }

  # print(plots)
  # Split plots into groups of 4
  plot_groups <- split(plots, ceiling(seq_along(plots) / 4))

  # Create a PDF with multiple pages
  pdf(glue("{outdir}/{pdf_suffix}.pdf"), width = 12, height = 6)
  print(glue("Histogram saved to {outdir}/{pdf_suffix}.pdf"))

  # Loop through groups and print each set of 4 to a new page
  walk(plot_groups, ~ grid.arrange(grobs = .x, nrow = 2, ncol = 2))

  dev.off() # Close the PDF device
  
  # Split the dataframe by group_var
  df_split <- split(df, df[[group_var]])
  
  # # Save each sub-dataframe as a CSV (for debugging, but this could also be good for supplemental figures)
  # sanitizes the group names in case they contain any character that isn't alphanumeric, a dash, or underscore
  if (export_locus_tags) {
    walk(df_split, ~ {
      group_name <- unique(.x[[group_var]])[1]
      safe_group_name <- gsub("[^A-Za-z0-9_\\-]", "_", group_name) # replace unsafe characters
      outpath <- glue("{outdir}/{pdf_suffix}_{safe_group_name}_locusTags.txt")
      write.table(.x$locus_tags, file = outpath, row.names=FALSE, col.names=FALSE, quote=FALSE)
    })
  }
}

histograms_by_source <- function(df, outdir, column_name = "sequence_length", reference_label=NA, consistent_y_axis = TRUE, group_var = "category", export_locus_tags = FALSE, width = 10) {
  make_dir(outdir)
  
  reference_length = NA # default value
  
  if (!is.na(reference_label)) {
    if ("sequence_id" %in% names(df)) {
      matching_rows <- which(df$sequence_id == reference_label)
    } else {
      # If sequence_id column doesn't exist, search through all character columns
      matching_rows <- integer(0)
      char_cols <- names(df)[sapply(df, is.character)]
      
      for (col in char_cols) {
        matches <- which(df[[col]] == reference_label)
        if (length(matches) > 0) {
          matching_rows <- matches
          print(glue("Found reference label '{reference_label}' in column '{col}'"))
          break
        }
      }
    }
    
    if (length(matching_rows) == 0) {
      print(glue("WARNING: Reference sequence ID '{reference_label}' not found in sequence_id column."))
      reference_length <- NA
    } else if (length(matching_rows) > 1) {
      print(glue("WARNING: Reference sequence ID '{reference_label}' appears {length(matching_rows)} times. Using first match."))
      reference_length <- df[[column_name]][matching_rows[1]]
      print(glue("Reference length for '{reference_label}': {reference_length}"))
    } else {
      reference_length <- df[[column_name]][matching_rows]
      print(glue("Reference length for '{reference_label}': {reference_length}"))
    }
    
  }
  
  # Generate histogram for overall data
  fname = glue("{outdir}/overall_histogram_{column_name}.pdf")
  if (!is.na(reference_label) & !is.na(reference_length)) {
    df |> ggplot() +
      geom_histogram(aes(x = df[[column_name]]), binwidth = width) +
      geom_vline(aes(xintercept=reference_length, colour="reference_length"), size=0.5) +
      geom_vline(aes(xintercept=mean(df[[column_name]], na.rm = TRUE), colour="mean"), size=0.5) +
      geom_vline(aes(xintercept=median(df[[column_name]], na.rm = TRUE), colour="median"), size=0.5) +
      labs(x = column_name, y = "Count", title = glue("{group_var}: Overall histogram; n: {nrow(df)}")) +
      scale_color_manual(name = "Values",
                           values = c(reference_length = "blue", mean = "red", median = "green"),
                           labels = c(reference_length = glue("{reference_label} = {reference_length}"),
                           mean = glue("Mean = {round(mean(df[[column_name]], na.rm = TRUE), 1)}"),
                           median = glue("Median = {median(df[[column_name]], na.rm = TRUE)}"))) +
      theme(legend.text = element_text(size = 8))
                           
      ggsave(fname, create.dir=TRUE, width = 8, height = 6)
      print(glue("Histogram saved to {fname}"))
  } else {
    df |> ggplot() +
      geom_histogram(aes(x = df[[column_name]]), binwidth = width) +
      geom_vline(aes(xintercept=mean(df[[column_name]], na.rm = TRUE), colour="mean"), size=0.5) +
      geom_vline(aes(xintercept=median(df[[column_name]], na.rm = TRUE), colour="median"), size=0.5) +
      labs(x = column_name, y = "Count", title = glue("{group_var}: Overall histogram; n: {nrow(df)}")) +
      scale_color_manual(name = "Values",
                           values = c(mean = "red", median = "green"),
                           labels = c(mean = glue("Mean = {round(mean(df[[column_name]], na.rm = TRUE), 1)}"),
                           median = glue("Median = {median(df[[column_name]], na.rm = TRUE)}")))
      ggsave(fname, create.dir=TRUE, width = 8, height = 6)
      print(glue("Histogram saved to {fname}"))
  }
    
  
  # Generate histogram PDF for group_var
  save_histograms(df, group_var, outdir, glue("{group_var}_histograms_{column_name}"), 
                  column_name = column_name, reference_label=reference_label, reference_length=reference_length, consistent_y_axis, export_locus_tags, width)
  
  if (export_locus_tags) {
    write.table(df$locus_tags, file=glue("{outdir}/{group_var}_histograms_locusTags.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
  }
  
  # If group_var is "category" or "genus", can further split the data (by "subcategory" or "species" respectively)
  if (group_var == "category") { # make additional plots for subcategory
    df <- df[!is.na(df$category), ] # remove rows with NA's for group_var so that the split will work properly
    df_categories <- split(df, df$category)
  
    make_dir(glue("{outdir}/subcategories"))
    # Generate subcategory-level histogram PDF
    walk(names(df_categories), ~ save_histograms(df_categories[[.x]], "subcategory", 
                                                outdir, glue("subcategories/{.x}_histograms_{column_name}"), 
                                                column_name = column_name, reference_label=reference_label, 
                                                reference_length=reference_length,
                                                consistent_y_axis, export_locus_tags, width))
  } else if (group_var == "genus") { # make additional plots for species
    df <- df[!is.na(df$genus), ] # remove rows with NA's for group_var so that the split will work properly
    df_categories <- split(df, df$genus)
  
    make_dir(glue("{outdir}/species"))
    # Generate species-level histogram PDF
    walk(names(df_categories), ~ save_histograms(df_categories[[.x]], "species", 
                                                outdir, glue("species/{.x}_histograms_{column_name}"), 
                                                column_name = column_name, reference_label=reference_label, 
                                                reference_length=reference_length,
                                                consistent_y_axis, export_locus_tags, width))
  }
}

# creates directories with that name if they don't already exist
make_dir <- function(new_dir) {
  ifelse(!dir.exists(file.path(new_dir)),
        dir.create(file.path(new_dir)),
        glue("{new_dir} directory exists"))
}