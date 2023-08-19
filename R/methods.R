#' @method print jointCalib
#' @exportS3Method
print.jointCalib <- function(x,...) {
  cat("Weights calibrated using: ", x$method, " with", x$backend, "backend.\n", sep="")
  cat("Summary statistics for g-weights:\n")

  print(summary(x$g))

  cat("Totals and precision (abs diff: ", sum(abs(x$diff)), ")\n", sep="")

  print(cbind(totals=x$totals, precision=x$diff))

  #cat("Quantiles and estimates based on g-weights:\n")

  #print(cbind(true = quants_known$x,
  #            est = weightedQuantile(df_resp$x, result1$g*df_resp$d, probs)))

  invisible(x)
}

