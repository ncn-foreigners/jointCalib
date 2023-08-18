#' @import mathjaxr
NULL
#' @title Weights calibration using empirical likelihood function
#'
#' @author Maciej Beręsewicz based on Wu (2005) and Zhang, Han and Wu (2022)
#'
#' \loadmathjax
#' @description
#' \code{calib_el} performs calibration using empirical likelihood method. The function is taken from Wu (2005), if algorithm has problem with convergence codes from Zhang, Han and Wu (2022) using \code{constrOptim} is used.
#'
#' @param X matrix of variables for calibration of quantiles and totals
#' @param d initial d-weights for calibration (e.g. design-weights)
#' @param totals vector of totals (where 1 element is population size)
#' @param maxit a numeric value giving the maximum number of iterations
#' @param tol the desired accuracy for the iterative procedure
#'
#' @references
#' Wu, C. (2005). Algorithms and R codes for the pseudo empirical likelihood method in survey sampling. Survey Methodology, 31(2), 239 (code is taken from \url{https://sas.uwaterloo.ca/~cbwu/Rcodes/LagrangeM2.txt}).
#' Zhang, S., Han, P., and Wu, C. (2023) Calibration Techniques Encompassing Survey Sampling, Missing Data Analysis and Causal Inference. International Statistical Review, 91: 165–192. https://doi.org/10.1111/insr.12518 (code is taken from Supplementary Materials).
#'
#' @return Returns a vector of empirical likelihood weights that sums to N
#'
#' @export
calib_el <- function(X, d, totals, maxit=50, tol=1e-8) {
  ## add constrOptim if something is not working
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
    DD <- -t(U) %*% (as.vector(d / (1 + U %*% lambda)^2) * U)
    D2 <- try(MASS::ginv(DD, tol = 1e-40) %*% D1)
    if (class(D2)[1] == "try-error") break
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
  if (k >= maxit | class(D2)[1] == "try-error") {

    message("Maximum iteration exceeded or MASS::ginv resulted in error. Changing to stats::constrOptim.")

    ll_fun <- function(d, rho, g_hat) {
      -sum(d*log(1 + g_hat %*% rho))
    }
    ll_grad <- function(d, rho, g_hat) {
      -colSums(d*g_hat / c(1 + g_hat %*% rho))
    }

    rho_hat <-
      stats::constrOptim(
        theta = rep(0, n_col),
        f = ll_fun,
        grad = ll_grad,
        ui = U,
        ci = rep(1 / n_row - 1, n_row),
        g_hat = U,
        d = d
      )
    lambda <- rho_hat$par
  }

  gweight <- c(d / (1 + U %*% lambda) / n_row)
  gweight <- gweight/sum(gweight)*N/d

  return(gweight)
}
