library(rstatix)
library(moments)
library(glue)

# stratified bootstrapping approach
bootstrap_stats <- function(df, group_var = "category", data_col = "sequence_length", n = 1000, seed = 42, downsample = FALSE) {
  set.seed(seed)
  dimensions <- paste(unlist(dim(df)), collapse = ", ")
  cat(glue("\nDimensions of df: {dimensions}\n"))
  cat("\n")
  df <- df[!is.na(df[[group_var]]), ] # filter out NA's, if any
  dimensions <- paste(unlist(dim(df)), collapse = ", ")
  cat(glue("\nDimensions of df after filtering out NA's in {group_var}, if any: {dimensions}\n"))
  cat("\n")
  
  # Pre-split data by group
  grouped_data <- split(df[[data_col]], df[[group_var]])
  group_sizes <- table(df[[group_var]])
  

  # Downsample each group to the smallest group size before bootstrapping
  if (downsample) {
    cat("\nDownsampling each group to the smallest group size.\n")
    min_size <- min(group_sizes)

    df <- df |>
      group_by(!!sym(group_var)) |>
      slice_sample(n = min_size) |>
      ungroup()
  } 
  
  cat("\nGroup sizes:")
  group_sizes <- table(df[[group_var]]) # get group_sizes again in case it changed
  print(group_sizes)
  cat("\n")

  # Bootstrap loop
  replicate(n, {
    resampled <- lapply(names(grouped_data), function(g) {
      sampled_values <- sample(grouped_data[[g]], size = group_sizes[[g]], replace = TRUE)
      tibble(
        !!group_var := g,
        mean = mean(sampled_values),
        sd = sd(sampled_values),
        skewness = moments::skewness(sampled_values),
        kurtosis = moments::kurtosis(sampled_values)
      )
    })

    bind_rows(resampled)
  }, simplify = FALSE) |>
    bind_rows(.id = "replicate") |>
    mutate(replicate = as.numeric(replicate)) |>
    pivot_longer(cols = c(mean, sd, skewness, kurtosis),
                 names_to = "statistic_name",
                 values_to = "statistic")
}

# function to perform hypothesis tests
compare_bootstrapped_stats <- function(bs_df, group_var = "category", adjust = "fdr", test_type = "non-parametric") {
  for (stat_name in unique(bs_df$statistic_name)) {
    df <- bs_df |> filter(statistic_name == stat_name)
    response <- "statistic"
    predictor <- group_var
    formula <- as.formula(glue("{response} ~ {predictor}"))

    print(glue("\n=== Comparing {stat_name} across {group_var} ==="))

    if (test_type == "parametric") {
      print("Running ANOVA test.")
      test <- aov(formula, data = df)
      p <- summary(test)[[1]][["Pr(>F)"]][1]
      print(glue("ANOVA p-value: {p}"))

      if (p < 0.05) {
        cat("Conducting a post-hoc Tukey HSD test.\n")
        follow_up <- TukeyHSD(test)
        print(follow_up)
      } else {
        cat("\nNo significant difference detected between predictor levels.\n")
      }

    } else {
      print("Running Kruskal-Wallis test.")
      test <- kruskal.test(formula, data = df)
      p <- test$p.value
      print(glue("Kruskal-Wallis p-value: {p}"))
      print(test)

      if (p < 0.05) {
        cat("Conducting a post-hoc Dunn's test.\n")
        follow_up <- dunn_test(data = df,
                               formula = formula,
                               p.adjust.method = adjust)
        print(follow_up)
      } else {
        cat("\nNo significant difference detected between predictor levels.\n")
      }
    }

    cat("")  # spacer
  }
}


plot_bootstrapped_stats <- function(original_df, bootstrapped_df, outdir, group_var = "category", data_col = "sequence_length", downsample = FALSE) {
  # Compute observed statistics
  obs_stats <- original_df |>
    group_by(.data[[group_var]]) |>
    summarise(
      mean = mean(.data[[data_col]]),
      sd = sd(.data[[data_col]]),
      skewness = moments::skewness(.data[[data_col]]),
      kurtosis = moments::kurtosis(.data[[data_col]]),
      .groups = "drop"
    ) |>
    pivot_longer(cols = -all_of(group_var), names_to = "statistic_name", values_to = "statistic")

  # Get group sizes
  group_sizes <- original_df |>
    group_by(!!sym(group_var)) |>
    summarise(size = n(), .groups = "drop")

  # If downsampling, use the smallest group size for all
  if (downsample) {
    group_sizes <- group_sizes |>
      mutate(size = min(size))
  }

  # Merge sizes into obs_stats for labeling
  obs_stats <- obs_stats |>
    left_join(group_sizes, by = group_var) |>
    mutate(facet_label = paste0(.data[[group_var]], "\n(n=", size, ")"))

  # Prepare bootstrapped data with matching facet labels
  bootstrapped_df <- bootstrapped_df |>
    left_join(group_sizes, by = group_var) |>
    mutate(facet_label = paste0(.data[[group_var]], "\n(n=", size, ")"))

  # Plot
  p <- ggplot(bootstrapped_df, aes(x = statistic)) +
    geom_histogram(bins = 30, fill = "gray80", color = "black") +
    geom_vline(data = obs_stats, aes(xintercept = statistic), color = "red", linewidth = 0.8) +
    facet_grid(rows = vars(facet_label), cols = vars(statistic_name), scales = "free_x") +
    labs(title = glue("Bootstrapped distributions of {data_col} statistics by {group_var}"),
         x = "Statistic value", y = "Count") +
    theme_bw()

  outname <- glue("{outdir}/bootstrapped_histograms_{group_var}_{data_col}.pdf")
  ggsave(outname, plot = p, width = 10, height = 8)
  cat(glue("\nSaved bootstrapped stat histograms to {outname}!\n"))
}



#bootstrap_stats <- function(df, group_var = "category", data_col = "sequence_length", n = 1000, seed = 42) {
#  set.seed(seed)
#  
#  # Pre-split data by group for efficient resampling
#  grouped_data <- split(df[[data_col]], df[[group_var]])
#
#  # Get group sizes
#  group_sizes <- table(df[[group_var]])
#
#  # Bootstrap loop
#  replicate(n, {
#    # Resample within each group
#    resampled <- lapply(names(grouped_data), function(g) {
#      sampled_values <- sample(grouped_data[[g]], size = group_sizes[[g]], replace = TRUE)
#      data.frame(
#        group = g,
#        mean = mean(sampled_values),
#        sd = sd(sampled_values),
#        skewness = moments::skewness(sampled_values),
#        kurtosis = moments::kurtosis(sampled_values)
#      )
#    })
#
#    # Combine group summaries for this replicate
#    do.call(rbind, resampled)
#  }, simplify = FALSE) |>
#    bind_rows(.id = "replicate") |>
#    mutate(replicate = as.numeric(replicate)) |>
#    pivot_longer(cols = c(mean, sd, skewness, kurtosis),
#                 names_to = "statistic_name",
#                 values_to = "statistic")
#}