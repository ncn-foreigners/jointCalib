#' @import mathjaxr
NULL
#' @title An internal function for calibration of weights using empirical likelihood method
#'
#' @author Maciej Beręsewicz based on Wu (2005) and Zhang, Han and Wu (2022)
#'
#' \loadmathjax
#' @description
#' \code{calib_el} performs calibration using empirical likelihood (EL) method. The function is taken from Wu (2005), if algorithm has problem with convergence codes from Zhang, Han and Wu (2022) using \code{constrOptim} is used.
#'
#' In (pseudo) EL the following (pseudo) EL function is maximized
#'
#' \mjsdeqn{\sum_{i \in r} d_i\log(p_i),}
#'
#' under the following constraint
#'
#' \mjsdeqn{\sum_{i \in r} p_i = 1,}
#'
#' with constraints on quantiles (with notation as in Harms and Duchesne (2006))
#'
#' \mjsdeqn{\sum_{i \in r} p_i(a_{i} - \alpha/N) = 0,}
#'
#' where \mjseqn{a_{i}} is created using \code{joint_calib_create_matrix} function, and possibly means
#'
#' \mjsdeqn{\sum_{i \in r} p_i(x_{i} - \mu_{x}) = 0,}
#'
#' where \mjseqn{\mu_{x}} is known population mean of X. For simplicity of notation we assume only one quantile and one mean is known. This can be generalized to multiple quantiles and means.
#'
#' @param X matrix of variables for calibration of quantiles and totals (first column should be intercept),
#' @param d initial d-weights for calibration (e.g. design-weights),
#' @param totals vector of totals (where 1 element is the population size),
#' @param maxit a numeric value giving the maximum number of iterations,
#' @param tol the desired accuracy for the iterative procedure,
#' @param eps the desired accuracy for computing the Moore-Penrose generalized inverse (see [MASS::ginv()]),
#' @param ... arguments passed to [stats::optim] via [stats::constrOptim].
#'
#' @references
#' Wu, C. (2005). Algorithms and R codes for the pseudo empirical likelihood method in survey sampling. Survey Methodology, 31(2), 239 (code is taken from \url{https://sas.uwaterloo.ca/~cbwu/Rcodes/LagrangeM2.txt}).
#'
#' Zhang, S., Han, P., and Wu, C. (2023) Calibration Techniques Encompassing Survey Sampling, Missing Data Analysis and Causal Inference. International Statistical Review, 91: 165–192. https://doi.org/10.1111/insr.12518 (code is taken from Supplementary Materials).
#'
#' @return Returns a vector of empirical likelihood g-weights
#'
#' @examples
#' ## generate data based on Haziza and Lesage (2016)
#' set.seed(123)
#' N <- 1000
#' x <- runif(N, 0, 80)
#' y <- exp(-0.1 + 0.1*x) + rnorm(N, 0, 300)
#' p <- rbinom(N, 1, prob = exp(-0.2 - 0.014*x))
#' totals_known <- c(N=N, x=sum(x))
#' df <- data.frame(x, y, p)
#' df_resp <- df[df$p == 1, ]
#' df_resp$d <- N/nrow(df_resp)
#' res <- calib_el(X = model.matrix(~x, df_resp),
#'                 d = df_resp$d,
#'                 totals = totals_known)
#' data.frame(known = totals_known, estimated=colSums(res*df_resp$d*model.matrix(~x, df_resp)))
#'
#' @export
calib_el <- function(X, d, totals, maxit=50, tol=1e-8, eps = .Machine$double.eps, ...) {
  n_col <- NCOL(X[, -1])
  n_row <- NROW(X)
  N <- totals[1]
  mu <- totals[-1]/N
  U <- X[, -1] - rep(1, n_row) %*% t(mu)
  lambda <- numeric(n_col)
  dif <- 1
  k <- 0
  while (dif > tol & k <= maxit) {
    p <- (1 + U %*% lambda)
    D1 <- t(U) %*% as.vector(d / p)
    DD <- -t(U) %*% (as.vector(d / p^2) * U)
    D2 <- try(MASS::ginv(DD, tol = eps) %*% D1)
    if (class(D2)[1] == "try-error") break
    dif <- max(abs(D2))
    ################################################
    ## here we use the same stoping rule as in sampling::calib
    ## tr <- crossprod(X, d/p)
    ## dif <- max(abs(tr - totals)/totals)
    ################################################
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
        method = "BFGS",
        ui = U,
        ci = rep(1 / n_row - 1, n_row),
        g_hat = U,
        d = d,
        control = list(maxit = maxit, ...)
      )
    lambda <- rho_hat$par
  }

  gweight <- as.numeric(1/(1 + U %*% lambda))
  return(gweight)
}
