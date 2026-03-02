library(rstatix)
library(glue)
# allow user to select from dropdown menu: numeric columns for response, categorical columns for predictor


quantResponse_categPredictor <- function(df, response = "sequence_length", predictor = "category",
                                         test_type = "non-parametric", adjust = "fdr") {
  output <- c()
  
  summary_df <- df |>
    group_by(!!sym(predictor)) |>
    summarize(count = n())
  
  output <- capture.output(print(summary_df, width = Inf))

  
  output <- c(output, glue("Test type: {test_type}"))
  output <- c(output, glue("P-value adjustment method: {adjust}"))

  if (test_type == "parametric") {
    output <- c(output, "Running ANOVA.")
    test <- aov(as.formula(glue("{response} ~ {predictor}")), data = df)
    p <- summary(test)[[1]][["Pr(>F)"]][1]
    output <- c(output, capture.output(summary(test)))
    #output <- c(output, capture.output(print(summary(test), max = 1000)))

    if (p < 0.05) {
      output <- c(output, "Conducting a post-hoc Tukey HSD test.")
      follow_up <- test |> tukey_hsd()
      output <- c(output, capture.output(as.data.frame(follow_up))) # convert to data frame; tibbles have truncated outputs
      
      # TODO: check the format of the tukey_hsd() output table and make a sig_pairs_df similar to the one for Dunn's test below
    } else {
      output <- c(output, "No significant difference detected between predictor levels.")
    }

  } else {
    output <- c(output, "Running Kruskal-Wallis test.")
    test <- kruskal.test(df[[response]] ~ df[[predictor]])
    p <- test$p.value
    output <- c(output, capture.output(test))

    if (p < 0.05) {
      output <- c(output, "Conducting a post-hoc Dunn's test.")
      follow_up <- dunn_test(data = df,
                             formula = as.formula(glue("{response} ~ {predictor}")),
                             p.adjust.method = adjust)
                             
      output <- c(output, "Significant pairs only:")
      sig_pairs_df <- follow_up |>
        select(group1, group2, n1, n2, p.adj, p.adj.signif) |>
        filter(p.adj.signif != "ns")
      output <- c(output, capture.output(as.data.frame(sig_pairs_df)))
      
      output <- c(output, "All pairs:")
      output <- c(output, capture.output(as.data.frame(follow_up))) # convert to data frame; tibbles have truncated outputs

    } else {
      output <- c(output, "No significant difference detected between predictor levels.")
    }
  }

  # Print all collected output
  return(paste(output, collapse = "\n"))
}


