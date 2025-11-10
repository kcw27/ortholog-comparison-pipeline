library(tidyverse)
library(glue)
theme_set(theme_classic())

# Function to save violin plots based on grouping variable
save_violin_plots <- function(df, group_var = "category", outdir, pdf_suffix, column_name = "sequence_length", reference_label=NA, reference_length=NA, additional = "boxplot") {
  if (!group_var %in% names(df)) stop(glue("Column '{group_var}' not found in df"))
  if (!column_name %in% names(df)) stop(glue("Column '{column_name}' not found in df"))
  
  df <- df[!is.na(df[[group_var]]), ] # remove rows with NA's for group_var 
  
  category_counts <- df |>
    group_by(!!sym(group_var)) |> # sym() converts string to symbol, and !! unquotes it
    summarise(n = n())
      
  fname = glue("{outdir}/{pdf_suffix}.pdf")
      
  # save the base plot, and add to it depending on options
  plots <- df |> ggplot() +
    geom_violin(aes(x=!!sym(group_var), y=!!sym(column_name))) +
    geom_text(data = category_counts, 
              aes(x=!!sym(group_var), y = max(df[[column_name]]) + 0.5, label = paste0("n = ", n)), 
              inherit.aes = FALSE, vjust = 0)  
  
  
  if (!is.na(reference_label) & !is.na(reference_length)) {
    # add the reference as a horizontal line
    plots <- plots + 
      geom_hline(aes(yintercept = reference_length), colour = 'blue') +
      labs(title = glue("{reference_label} {column_name} = {reference_length}"))
    
#    plots <- plots + 
#      geom_hline(aes(yintercept = reference_length, colour = "value")) +
#      scale_color_manual(name = "Values",
#          values = c("value" = "blue"),  # match factor levels as strings
#          labels = c("value" = glue("{reference_label} {column_name} = {reference_length}"))
#    )
  }
  
  if (additional == "boxplot") {
    plots <- plots + 
      geom_boxplot(aes(x=!!sym(group_var), y=!!sym(column_name)), width=0.1, color = "red")
  } else if (additional == "mean") {
    plots <- plots + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red")
  }
  # if the additional argument doesn't match either, ignore it
  
  ggsave(fname, plot = plots, create.dir=TRUE)
  print(glue("Violin plots saved to {fname}"))
}

violin_plots_by_source <- function(df, outdir, column_name = "sequence_length", reference_label=NA, group_var = "category", additional = "boxplot") {
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
  

  # Generate violin plot for group_var
  save_violin_plots(df, group_var, outdir, glue("{group_var}_violinPlots_{column_name}"), 
                  column_name = column_name, reference_label=reference_label, reference_length=reference_length, additional)
  
  
  # If group_var is "category" or "genus", can further split the data (by "subcategory" or "species" respectively)
  if (group_var == "category") { # make additional plots for subcategory
    df <- df[!is.na(df$category), ] # remove rows with NA's for group_var so that the split will work properly
    df_categories <- split(df, df$category)
  
    make_dir(glue("{outdir}/subcategories"))
    # Generate subcategory-level violin plots
    walk(names(df_categories), ~ save_violin_plots(df_categories[[.x]], "subcategory", 
                                                outdir, glue("subcategories/{.x}_violinPlots_{column_name}"), 
                                                column_name = column_name, reference_label=reference_label, 
                                                reference_length=reference_length, additional))
  } else if (group_var == "genus") { # make additional plots for species
    df <- df[!is.na(df$genus), ] # remove rows with NA's for group_var so that the split will work properly
    df_categories <- split(df, df$genus)
  
    make_dir(glue("{outdir}/species"))
    # Generate species-level violin plots
    walk(names(df_categories), ~ save_violin_plots(df_categories[[.x]], "species", 
                                                outdir, glue("species/{.x}_violinPlots_{column_name}"), 
                                                column_name = column_name, reference_label=reference_label, 
                                                reference_length=reference_length, additional))
  }
}

# creates directories with that name if they don't already exist
make_dir <- function(new_dir) {
  ifelse(!dir.exists(file.path(new_dir)),
        dir.create(file.path(new_dir)),
        glue("{new_dir} directory exists"))
}