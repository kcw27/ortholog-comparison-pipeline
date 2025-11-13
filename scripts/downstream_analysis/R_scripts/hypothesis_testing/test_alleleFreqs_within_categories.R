library(dplyr)
library(tidyr)
library(rstatix)

test_alleleFreqs_within_categories <- function(df, site, residue, group_var = "category", seq_colname = "sequence", adjust = "fdr") {
  # Made a more concise version.
  # df[[seq_colname]] is assumed to contain aligned sequences
  # site is assumed to be a single int, pre-adjusted upon entry into the GUI
  adj_site <- site
  
  df |>
    group_by(!!sym(group_var)) |>
    summarize(count=n()) |>
    print()

  ### make contingency table for pairwise Fisher's test
  summary_df <- df |> 
    mutate(allele = substr(df[[seq_colname]], adj_site, adj_site)) |> 
    mutate(is_allele = (allele == residue)) |> 
    group_by(!!sym(group_var), is_allele) |> 
    summarize(n = n(), .groups = "drop") 
  
  # Ensure all combinations of group_var and is_allele exist
  summary_df <- summary_df |>
    complete(!!sym(group_var), is_allele = c(TRUE, FALSE), fill = list(n = 0))

  # get group names in the same order as in summary_df
  groups <- summary_df |> filter(is_allele == TRUE) |> pull(!!sym(group_var))
  observed <- summary_df |> filter(is_allele == TRUE) |> pull(n)
  not_observed <- summary_df |> filter(is_allele == FALSE) |> pull(n)
  
  xtab <- as.table(rbind(observed, not_observed))
  dimnames(xtab) <- list(allele_observed = c("yes", "no"), group = groups)
  print(xtab)
  
  ### NEW: Calculate overall statistics for enrichment tests
  total_sequences <- sum(xtab)
  total_allele <- sum(xtab["yes", ])
  overall_prop <- total_allele / total_sequences
  
  cat("\n=== OVERALL STATISTICS ===\n")
  cat("Total sequences:", total_sequences, "\n")
  cat("Total", residue, "alleles:", total_allele, "\n")
  cat(glue("Overall proportion: {round(overall_prop, 4)} ( {round(overall_prop * 100, 1)}% )"), "\n")
  
  cat("\n=== ENRICHMENT TESTS ===\n")
  enrichment_results <- data.frame(
    category = groups,
    total_sequences = colSums(xtab),
    observed_alleles = xtab["yes", ],
    expected_alleles = colSums(xtab) * overall_prop,
    observed_prop = xtab["yes", ] / colSums(xtab),
    p_value_fisher_2x2 = NA
  )
  
  # Perform tests for each category
  for(i in 1:nrow(enrichment_results)) {    
    # 2x2 Fisher's exact test for this category vs all others
    category_A <- enrichment_results$observed_alleles[i]
    category_nonA <- enrichment_results$total_sequences[i] - category_A
    other_A <- total_allele - category_A
    other_nonA <- (total_sequences - enrichment_results$total_sequences[i]) - other_A
    
    fisher_2x2 <- fisher.test(
      matrix(c(category_A, category_nonA, other_A, other_nonA), nrow = 2),
      alternative = "greater"
    )
    enrichment_results$p_value_fisher_2x2[i] <- fisher_2x2$p.value
  }
  
  # Adjust p-values for multiple testing
  enrichment_results$p_adj_fisher_2x2 <- p.adjust(enrichment_results$p_value_fisher_2x2, method = adjust)
  
  # Add significance flags
  enrichment_results$significant_fisher_2x2 <- enrichment_results$p_adj_fisher_2x2 < 0.05
  
  # Format and print results
  enrichment_results_formatted <- enrichment_results |>
    mutate(
      observed_prop = round(observed_prop, 4),
      expected_alleles = round(expected_alleles, 2),
      p_value_fisher_2x2 = format.pval(p_value_fisher_2x2, digits = 3),
      p_adj_fisher_2x2 = format.pval(p_adj_fisher_2x2, digits = 3)
    )
  
  row.names(enrichment_results_formatted) <- NULL
  print(enrichment_results_formatted)
  
  
  # Overall test - is there any significant variation?
  cat("\nFisher's exact test:\n")
  print(test <- fisher_test(xtab, detailed = TRUE)) #, alternative = "greater")) # since this is just to determine whether there's variation between the groups, would it even make sense to use alternative = "greater" here?
  p <- test$p
  
  if (p < 0.05) {
    cat("\np from Fisher's exact test is significant; conducting pairwise comparisons\n")
    # Pairwise tests - which specific pairs are different?
    cat("Pairwise Fisher's tests:\n")
    pairwise_results <- pairwise_fisher_test(xtab, p.adjust.method = adjust, alternative = "greater")
    print(pairwise_results)
    
    # Identify significantly different pairs
    sig_pairs <- pairwise_results |> filter(p.adj < 0.05)
    if(nrow(sig_pairs) > 0) {
      cat("\nSignificantly different category pairs (adjusted p < 0.05):\n")
      print(sig_pairs)
    } else {
      cat("\nNo significantly different category pairs after adjustment\n")
    }
  } else {
    cat("\nOverall p was not significant\n")
  }
  
  ### Summary of findings
  cat("\n=== SUMMARY ===\n")
  cat("Categories significantly enriched for", residue, "allele:\n")
  enriched_categories <- enrichment_results |> 
    filter(significant_fisher_2x2)
  if(nrow(enriched_categories) > 0) {
    for(i in 1:nrow(enriched_categories)) {
      cat("  -", enriched_categories$category[i], 
          "(observed:", enriched_categories$observed_alleles[i], "/", enriched_categories$total_sequences[i],
          "=", round(enriched_categories$observed_prop[i] * 100, 1), "%; expected:", 
          round(enriched_categories$expected_alleles[i], 1), 
          ", p-value adjusted by", adjust, "is (p =", formatC(enriched_categories$p_adj_fisher_2x2[i], format = "e", digits = 2), ") )\n")
          
      # additional print statements to facilitate data entry
      cat(glue("{enriched_categories$category[i]} (p = {formatC(enriched_categories$p_adj_fisher_2x2[i], format = 'e', digits = 2)}; {round(enriched_categories$observed_prop[i] * 100, 1)}%)"), "\n") 
      # cat(glue("expected {round(enriched_categories$expected_alleles[i] / enriched_categories$total_sequences[i] * 100, 1)}%"))
      cat(glue("expected {round(overall_prop * 100, 1)}%"), "\n\n")
    }
  } else {
    cat("  None\n")
  }
  
  # Return all results for potential further analysis
  invisible(list(
    contingency_table = xtab,
    enrichment_results = enrichment_results,
    overall_fisher_p = p,
    total_sequences = total_sequences,
    total_allele = total_allele,
    overall_proportion = overall_prop
  ))
}

#test_alleleFreqs_within_categories <- function(df, site, residue, group_var = "category", seq_colname = "sequence", adjust = "fdr") {
#  # df[[seq_colname]] is assumed to contain aligned sequences
#  # site is assumed to be a single int, pre-adjusted upon entry into the GUI
#  adj_site <- site
#
#  ### make contingency table for pairwise Fisher's test
#  summary_df <- df |> 
#    mutate(allele = substr(df[[seq_colname]], adj_site, adj_site)) |> 
#    mutate(is_allele = (allele == residue)) |> 
#    group_by(!!sym(group_var), is_allele) |> 
#    summarize(n = n(), .groups = "drop") 
#  
#  # Ensure all combinations of group_var and is_allele exist
#  summary_df <- summary_df |>
#    complete(!!sym(group_var), is_allele = c(TRUE, FALSE), fill = list(n = 0))
#
#  # get group names in the same order as in summary_df
#  groups <- summary_df |> filter(is_allele == TRUE) |> pull(!!sym(group_var))
#  observed <- summary_df |> filter(is_allele == TRUE) |> pull(n)
#  not_observed <- summary_df |> filter(is_allele == FALSE) |> pull(n)
#  
#  xtab <- as.table(rbind(observed, not_observed))
#  dimnames(xtab) <- list(allele_observed = c("yes", "no"), group = groups)
#  print(xtab)
#  
#  ### NEW: Calculate overall statistics for enrichment tests
#  total_sequences <- sum(xtab)
#  total_allele <- sum(xtab["yes", ])
#  overall_prop <- total_allele / total_sequences
#  
#  cat("\n=== OVERALL STATISTICS ===\n")
#  cat("Total sequences:", total_sequences, "\n")
#  cat("Total", residue, "alleles:", total_allele, "\n")
#  cat("Overall proportion:", round(overall_prop, 4), "(", round(overall_prop * 100, 2), "%)\n")
#  
#  ### NEW: Binomial enrichment tests for each category
#  # The binomial test asks: "If A's are randomly distributed, what's the probability of seeing 3 or more A's in a sample of 66 sequences?"
#  cat("\n=== BINOMIAL ENRICHMENT TESTS ===\n")
#  enrichment_results <- data.frame(
#    category = groups,
#    total_sequences = colSums(xtab),
#    observed_alleles = xtab["yes", ],
#    expected_alleles = colSums(xtab) * overall_prop,
#    observed_prop = xtab["yes", ] / colSums(xtab),
#    p_value_binomial = NA,
#    p_value_fisher_2x2 = NA
#  )
#  
#  # Perform binomial tests for each category
#  for(i in 1:nrow(enrichment_results)) {
#    # Binomial test for enrichment (one-sided)
#    binom_test <- binom.test(
#      x = enrichment_results$observed_alleles[i],
#      n = enrichment_results$total_sequences[i],
#      p = overall_prop,
#      alternative = "greater"
#    )
#    enrichment_results$p_value_binomial[i] <- binom_test$p.value
#    
#    # 2x2 Fisher's exact test for this category vs all others
#    category_A <- enrichment_results$observed_alleles[i]
#    category_nonA <- enrichment_results$total_sequences[i] - category_A
#    other_A <- total_allele - category_A
#    other_nonA <- (total_sequences - enrichment_results$total_sequences[i]) - other_A
#    
#    fisher_2x2 <- fisher.test(
#      matrix(c(category_A, category_nonA, other_A, other_nonA), nrow = 2),
#      alternative = "greater"
#    )
#    enrichment_results$p_value_fisher_2x2[i] <- fisher_2x2$p.value
#  }
#  
#  # Adjust p-values for multiple testing
#  enrichment_results$p_adj_binomial <- p.adjust(enrichment_results$p_value_binomial, method = adjust)
#  enrichment_results$p_adj_fisher_2x2 <- p.adjust(enrichment_results$p_value_fisher_2x2, method = adjust)
#  
#  # Add significance flags
#  enrichment_results$significant_binomial <- enrichment_results$p_adj_binomial < 0.05
#  enrichment_results$significant_fisher_2x2 <- enrichment_results$p_adj_fisher_2x2 < 0.05
#  
#  # Format and print results
#  enrichment_results_formatted <- enrichment_results |>
#    mutate(
#      observed_prop = round(observed_prop, 4),
#      expected_alleles = round(expected_alleles, 2),
#      p_value_binomial = format.pval(p_value_binomial, digits = 3),
#      p_adj_binomial = format.pval(p_adj_binomial, digits = 3),
#      p_value_fisher_2x2 = format.pval(p_value_fisher_2x2, digits = 3),
#      p_adj_fisher_2x2 = format.pval(p_adj_fisher_2x2, digits = 3)
#    )
#  
#  row.names(enrichment_results_formatted) <- NULL
#  print(enrichment_results_formatted)
#  
#  ### Original tests for overall variation
#  cat("\n=== OVERALL ASSOCIATION TESTS ===\n")
#  cat("Chi-squared test of independence:")
#  print(chisq_test(xtab))
#  
#  # Overall test - is there any significant variation?
#  cat("\nFisher's exact test:\n")
#  print(test <- fisher_test(xtab, detailed = TRUE)) #, alternative = "greater")) # since this is just to determine whether there's variation between the groups, would it even make sense to use alternative = "greater" here?
#  p <- test$p
#  
#  if (p < 0.05) {
#    cat("\np from Fisher's exact test is significant; conducting pairwise comparisons\n")
#    # Pairwise tests - which specific pairs are different?
#    cat("Pairwise Fisher's tests:\n")
#    pairwise_results <- pairwise_fisher_test(xtab, p.adjust.method = adjust, alternative = "greater")
#    print(pairwise_results)
#    
#    # Identify significantly different pairs
#    sig_pairs <- pairwise_results |> filter(p.adj < 0.05)
#    if(nrow(sig_pairs) > 0) {
#      cat("\nSignificantly different category pairs (adjusted p < 0.05):\n")
#      print(sig_pairs)
#    } else {
#      cat("\nNo significantly different category pairs after adjustment\n")
#    }
#  } else {
#    cat("\nOverall p was not significant\n")
#  }
#  
#  ### Summary of findings
#  cat("\n=== SUMMARY ===\n")
#  cat("Categories significantly enriched for", residue, "allele:\n")
#  enriched_categories <- enrichment_results |> 
#    filter(significant_binomial | significant_fisher_2x2)
#  if(nrow(enriched_categories) > 0) {
#    for(i in 1:nrow(enriched_categories)) {
#      cat("  -", enriched_categories$category[i], 
#          "(observed:", enriched_categories$observed_alleles[i], "/", enriched_categories$total_sequences[i],
#          "=", round(enriched_categories$observed_prop[i] * 100, 1), "%; expected:", 
#          round(enriched_categories$expected_alleles[i], 1), ")\n")
#    }
#  } else {
#    cat("  None\n")
#  }
#  
##  # Code recommended by Copilot: the pre- and post- adjustment pvals look like p_value_binomial and p_adj_binomial from the earlier binomial test,
##  # so this code is redundant.
##  # Total number of "yes" alleles
##  total_yes <- sum(xtab["yes", ])
##  # Total number of samples
##  total_samples <- colSums(xtab)
##  
##  # Expected proportion of "A"
##  expected_p <- total_yes / sum(total_samples)
##  
##  # Binomial tests per group
##  pvals <- sapply(1:5, function(i) {
##    binom.test(xtab["yes", i], total_samples[i], p = expected_p, alternative = "greater")$p.value
##  })
##  names(pvals) <- colnames(xtab)
##  print(pvals)
##  print(p.adjust(pvals, method = adjust)) 
#
#  
#  # Return all results for potential further analysis
#  invisible(list(
#    contingency_table = xtab,
#    enrichment_results = enrichment_results,
#    overall_fisher_p = p,
#    total_sequences = total_sequences,
#    total_allele = total_allele,
#    overall_proportion = overall_prop
#  ))
#}


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