#' @import mathjaxr
NULL
#' @title Weights calibration using empirical likelihood function
#'
#' @author Maciej Beręsewicz based on Zhang, Han and Wu (2022)
#'
#' \loadmathjax
#' @description
#' \code{calib_el} performs calibration using empirical likelihood method. The function is taken from Wu (2005)
#'
#' @param X matrix of variables for calibration of quantiles and totals
#' @param d initial d-weights for calibration (e.g. design-weights)
#' @param totals vector of totals (where 1 element is population size)
#' @param maxit a numeric value giving the maximum number of iterations
#' @param tol the desired accuracy for the iterative procedure
#'
#' @references
#' Wu, C. (2005). Algorithms and R codes for the pseudo empirical likelihood method in survey sampling. Survey Methodology, 31(2), 239..
#'
#' @return Returns a vector of empirical likelihood weights that sums to N
#'
#' @export
calib_el <- function(X, d, totals, maxit=50, tol=1e-8) {

  n_col <- NCOL(X[, -1])
  n_row <- NROW(X)
  N <- totals[1]
  mu <- totals[-1]/N
  U <- X[, -1] - rep(1, n_row) %*% t(mu)
  lambda <- numeric(n_col)
  dif <- 1
  k <- 0
  while (dif > tol & k <= maxit) {
    D1 <- t(U) %*% as.vector(d / (1 + U %*% lambda))
    DD <- -t(U) %*% (as.vector(d / (1 + U %*% lambda) ^ 2) * U)
    D2 <- solve(DD, D1, tol = 1e-40)
    dif <- max(abs(D2))
    rule <- 1
    while (rule > 0) {
      rule <- 0
      if (min(1 + t(lambda - D2) %*% t(U)) <= 0)
        rule <- rule + 1
      if (rule > 0)
        D2 <- D2 / 2
    }
    lambda <- lambda - D2
    k <- k + 1
  }
  if (k >= maxit)
    stop("Maximum number of iterations reached")

  gweight <- c(d / (1 + U %*% lambda) / n_row)
  gweight <- gweight/sum(gweight)*N/d

  return(gweight)
}
