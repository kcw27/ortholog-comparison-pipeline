library(rstatix)
library(glue)
# allow user to select from dropdown menu: numeric columns for response, categorical columns for predictor


quantResponse_categPredictor <-function(df, response="sequence_length", predictor="category", test_type = "non-parametric", adjust = "fdr") {
  # response: string specifying name of numeric column
  # predictor: string specifying name of categorical column
  print(glue("Test type: {test_type}"))
  print(glue("P-value adjustment method: {adjust}"))
  
  if (test_type == "parametric") {
    # do ANOVA aov(), and if significant, TukeyHSD() follow-up
    print("Running ANOVA.")
    test <- test <- aov(as.formula(glue("{response} ~ {predictor}")), data = df)
    p <- summary(test)[[1]][["Pr(>F)"]][1]
    print(summary(test))

    if (p < 0.05) {
      print("Conducting a post-hoc Tukey HSD test.")
      follow_up <- test |> tukey_hsd()
      print(follow_up)
    } else {
      print("No significant difference detected between predictor levels.")
    }
    
  } else { # have to use non-parametric tests; less power
    # first obtain p.val from kruskal, then if p<0.05, do dunn follow-up
    print("Running Kruskal-Wallis test.")
    test <- kruskal.test(df[[response]]~df[[predictor]])
    p <- test$p.val
    print(test)
    
    if (p < 0.05) {
      print("Conducting a post-hoc Dunn's test.")
      follow_up <- dunn_test(data = df, formula = as.formula(glue("{response} ~ {predictor}")), p.adjust.method = adjust)
      print(follow_up)
    } else {
      print("No significant difference detected between predictor levels.")
    }
  }
}

