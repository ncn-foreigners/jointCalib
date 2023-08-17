#' @import mathjaxr
NULL
#' @title Weights calibration using empirical likelihood function
#'
#' @author Maciej Beręsewicz based on Zhang, Han and Wu (2022)
#'
#' \loadmathjax
#' @description
#' \code{calib_el} performs calibration using empirical likelihood method. The function is taken from the supplementary materials from the publication Zhang, Han and Wu (2022)
#'
#' @param X matrix of variables for calibration of quantiles
#' @param d initial d-weights for calibration (e.g. design-weights)
#' @param totals vector of totals
#' @param N population size
#' @param ... arguments passed to \code{stats::constrOptim}
#'
#' @references
#' Zhang, S., Han, P., and Wu, C. (2023) Calibration Techniques Encompassing Survey Sampling, Missing Data Analysis and Causal Inference. International Statistical Review, 91: 165–192. https://doi.org/10.1111/insr.12518.
#'
#' @return Returns a vector of empirical likelihood weights that sums to 1
#'
#' @export
calib_el <- function(X, d, totals, N, ...) {

  n_col <- NCOL(X)
  n_row <- NROW(X)

  Xs <- X - do.call("rbind", replicate(n_row, totals/N, simplify = F))

  log_lik <- function(rho, g_hat) {
    -sum(log(1 + g_hat %*% rho))
  }
  grad <- function(rho, g_hat) {
    -colSums(g_hat / c(1 + g_hat %*% rho))
  }

  rho.hat <-
    stats::constrOptim(
      theta = rep(0, n_col),
      f = log_lik,
      grad = grad,
      ui = Xs,
      ci = rep(1 / n_row - 1, n_row),
      g_hat = Xs,
      ...
    )$par

  weight <- c(1 / (1 + Xs %*% rho.hat) / n_row)
  weight <- weight / sum(weight)

  return(weight)
}
