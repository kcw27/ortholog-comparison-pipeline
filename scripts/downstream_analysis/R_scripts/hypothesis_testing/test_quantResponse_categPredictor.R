
# allow user to select from dropdown menu: numeric columns for response, categorical columns for predictor


quantResponse_categPredictor <-function(df, response="sequence_length", predictor="category", test_type = "non-parametric") {
  # response: string specifying name of numeric column
  # predictor: string specifying name of categorical column
  if (test_type == "parametric") {
    # do ANOVA aov(), and if significant, TukeyHSD() follow-up
  } else { # have to use non-parametric tests; less power
    # first obtain p.val from kruskal, then if p<0.05, do dunn follow-up
    # print(kruskal.test(sequence_length~category, data=df[df$category != "no category",]))
  }
}

