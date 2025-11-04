library(tidyverse)
library(stringdist)
library(glue)
theme_set(theme_bw())

get_pairwise_dists <- function(df, group_var = "category", seq_colname = "sequence") {
  df_categories <- split(df, df[[group_var]])

  # Collect partial dataframes
  all_dists <- map_dfr(names(df_categories), function(g) {
    # pull out the set of sequences corresponding to the current group
    seqs <- df_categories[[g]][[seq_colname]]
    
    # get pairwise distances within the category
    mtx <- stringdistmatrix(a = seqs, b = seqs, method = "osa")
    pairwise_dists <- mtx[upper.tri(mtx)]
    
    tibble(pairwise_dist = pairwise_dists, group := g)
  })

  return(all_dists)
}

compare_pairwise_dists <- function(dist_df, test_type = "non-parametric", adjust = "fdr") {  
  if (test_type == "parametric") {
    print("Running ANOVA test.")
    test <- aov(pairwise_dist ~ group, data = dist_df)
    p <- summary(test)[[1]][["Pr(>F)"]][1]
    print(glue("ANOVA p-value: {signif(p, 4)}"))

    if (p < 0.05) {
      print("Conducting a post-hoc Tukey HSD test.")
      follow_up <- TukeyHSD(test)
      print(follow_up)
    } else {
      print("No significant difference detected between predictor levels.")
    }

  } else {
    print("Running Kruskal-Wallis test.")
    test <- kruskal.test(pairwise_dist ~ group, data = dist_df)
    p <- test$p.value
    print(glue("Kruskal-Wallis p-value: {signif(p, 4)}"))
    print(test)

    if (p < 0.05) {
      print("Conducting a post-hoc Dunn's test.")
      follow_up <- dunn_test(data = dist_df,
                             formula = pairwise_dist ~ group,
                             p.adjust.method = adjust)
      print(follow_up)
    } else {
      print("No significant difference detected between predictor levels.")
    }
  }
}

plot_pairwise_dists <- function(dist_df, outdir, group_var = "category") {
  # Compute group-wise means
  mean_lines <- dist_df |>
    group_by(group) |>
    summarise(mean_dist = mean(pairwise_dist), .groups = "drop")
  
  # Plot with histogram and mean overlay
  p <- ggplot(dist_df, aes(x = pairwise_dist)) +
    geom_histogram(bins = 30, color = "white") + # change the outline with the color argument for legibility
    geom_vline(data = mean_lines, aes(xintercept = mean_dist, color = "Mean pairwise distance")) +
    scale_color_manual(name = "", values = c("Mean pairwise distance" = "red")) +
    facet_wrap(vars(group), scales = "free_y") +
    labs(title = glue("Pairwise sequence distances by {group_var}"),
         x = "Pairwise distance", y = "Count") +
    theme_bw()
    
  # save figure
  outname = glue("{outdir}/pairwiseSeqDists_within_{group_var}.pdf")
  ggsave()
}

#get_pairwise_dists <- function(df, group_var = "category", seq_colname = "sequence") {
#  df_categories <- split(df, df[[group_var]])
#  
#  walk(names(df_categories), ~ {
#    # pull out the set of sequences corresponding to the current group
#    indices <- which(df[[group_var]] == .x)
#    seqs <- df[indices, seq_colname]
#
#    # get pairwise distances within the category
#    mtx <- stringdistmatrix(a = seqs, b = seqs, method = "osa")
#    pairwise_dists <- mtx[upper.tri(mtx)]
#    
#    partial_df <- data.frame(pairwise_dists, .x)
#    colnames(partial_df) = c("pairwise_dist", group_var)
#    #print(head(partial_df))
#    #print(dim(partial_df))
#  })
#}