library(ggplot2)
library(glue)

test_normality <- function(df, colname, fig_dir) {
  outfile = glue("{fig_dir}/qqplot_{colname}.pdf")
  
  # Ensure the directory exists
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)

#  pdf(outfile, width = 8, height = 6)
#  qqnorm(df[[colname]], main = paste("QQ Plot of", colname))
#  qqline(df[[colname]])
#  dev.off()
#  print(glue("QQ plot saved to {outfile}"))
  
  pngfile = glue("{fig_dir}/qqplot_{colname}.png")
  png(pngfile, width = 800, height = 600)
  qqnorm(df[[colname]], main = paste("QQ Plot of", colname))
  qqline(df[[colname]])
  dev.off()
  print(glue("QQ plot saved to {pngfile}"))

  set.seed(42)
  if (length(df[[colname]]) > 5000) { # shapiro.test() doesn't work on datasets larger than 5000 values
    subset <- sample(df[[colname]], size = 5000)
    normality_test <- shapiro.test(subset)
  } else {
    normality_test <- shapiro.test(df[[colname]])
  }
  
  print(normality_test)
  
  if (normality_test$p.val < 0.05) {
    print(glue("{colname} is probably not normal; p={signif(normality_test$p.val, 3)}"))
    test_type = "non-parametric" 
  } else {
    print(glue("{colname} is probably normal; p={signif(normality_test$p.val, 3)}"))
    test_type = "parametric"
  }
  return(test_type)
}
