#' @title Print method for jointCalib class
#' @param x object of jointCalib class
#' @param ... Additional optional arguments
#' @method print jointCalib
#' @exportS3Method
print.jointCalib <- function(x,...) {
  cat("Weights calibrated using: ", x$method, " (backend: ", x$backend,")\n", sep="")
  cat("Summary statistics for g-weights:\n")

  print(summary(x$g))

  cat("Totals and precision (abs diff: ", sum(abs(x$diff)), ")\n", sep="")

  print(cbind(totals=x$totals, precision=x$diff))
  invisible(x)
}

