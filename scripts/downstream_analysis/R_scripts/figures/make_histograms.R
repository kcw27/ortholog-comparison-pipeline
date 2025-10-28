library(tidyverse)
library(glue)
library(gridExtra)
library(purrr)
theme_set(theme_bw())

# Function to save histograms based on grouping variable
save_histograms <- function(df, group_var, outdir, pdf_suffix, column_name = "sequence_length", reference_length=NA, reference_label=NA, consistent_y_axis = TRUE, export_locus_tags = FALSE, width = 10) {
  if (!group_var %in% names(df)) stop(glue("Column '{group_var}' not found in df"))
  if (!column_name %in% names(df)) stop(glue("Column '{column_name}' not found in df"))
  
  df <- df[!is.na(df[[group_var]]), ] # remove rows with NA's for group_var so that the split will work properly

  ## To create a set of plots with consistent x and y axes, calculate the x and y limits.
  # x limits for a histogram are just the range of values
  x_limits = range(df[[column_name]])
  buffer_size = 0.10 # make buffers that are some percentage of the range
  x_limits <- x_limits + c(-buffer_size*diff(x_limits), buffer_size*diff(x_limits)) 
  # add a bit of buffer on either side so that it doesn't crop out the data at either end of the x axis
  
  if (!is.na(reference_length) & !is.na(reference_label)) {
    plots <- split(df, df[[group_var]]) |>
      map(~ ggplot(.x) +
            geom_histogram(aes(x = .data[[column_name]]), binwidth = width) +
            geom_vline(aes(xintercept=reference_length, colour="reference_length"), size=0.5) +
            geom_vline(aes(xintercept=mean(.data[[column_name]]), colour="mean"), size=0.5) +
            labs(x = column_name, y = "Count", title = glue("{group_var}: {unique(.x[[group_var]])[1]}; n: {nrow(.x)}")) +
            scale_color_manual(name = "Values",
                               values = c(reference_length = "blue", mean = "red"),
                               labels = c(reference_length = glue("{reference_label} = {reference_length}"),
                                          mean = glue("Mean = {round(mean(.x[[column_name]]), 1)}"))) +
            xlim(x_limits) # notice that we set x limits here, but no y limits yet
          )
  } else {
    plots <- split(df, df[[group_var]]) |>
      map(~ ggplot(.x) +
            geom_histogram(aes(x = .data[[column_name]]), binwidth = width) +
            geom_vline(aes(xintercept=mean(.data[[column_name]]), colour="mean"), size=0.5) +
            labs(x = column_name, y = "Count", title = glue("{group_var}: {unique(.x[[group_var]])[1]}; n: {nrow(.x)}")) +
            scale_color_manual(name = "Values",
                               values = c(mean = "red"),
                               labels = c(mean = glue("Mean = {round(mean(.x[[column_name]]), 1)}"))) +
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

histograms_by_source <- function(df, outdir, column_name = "sequence_length", reference_length=NA, reference_label=NA, consistent_y_axis = TRUE, group_var = "category", export_locus_tags = FALSE, width = 10) {
  make_dir(outdir)
  # Generate histogram PDF for group_var
  save_histograms(df, group_var, outdir, glue("{group_var}_histograms"), 
                  column_name = column_name, reference_length=reference_length, reference_label=reference_label, consistent_y_axis, export_locus_tags, width)
  
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
                                                outdir, glue("subcategories/{.x}_histograms"), 
                                                column_name = column_name, reference_length=reference_length, reference_label=reference_label, 
                                                consistent_y_axis, export_locus_tags, width))
  } else if (group_var == "genus") { # make additional plots for species
    df <- df[!is.na(df$genus), ] # remove rows with NA's for group_var so that the split will work properly
    df_categories <- split(df, df$genus)
  
    make_dir(glue("{outdir}/species"))
    # Generate species-level histogram PDF
    walk(names(df_categories), ~ save_histograms(df_categories[[.x]], "species", 
                                                outdir, glue("species/{.x}_histograms"), 
                                                column_name = column_name, reference_length=reference_length, reference_label=reference_label, 
                                                consistent_y_axis, export_locus_tags, width))
  }
}

# creates directories with that name if they don't already exist
make_dir <- function(new_dir) {
  ifelse(!dir.exists(file.path(new_dir)),
        dir.create(file.path(new_dir)),
        glue("{new_dir} directory exists"))
}