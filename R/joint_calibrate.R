#' @title Function to calibrate totals and quantiles jointly
#' @author Maciej Beręsewicz
#'
#' @description
#' \code{joint_calibrate} allows to specify matrix of X variables for calibration of totals (intercept should not be included) and matrix of X_q variables for calibration of quantiles.
#'
#' @param X matrix of variables for calibration of totals
#' @param X_q matrix of variables for calibration of quantiles
#' @param d initial weights for calibration (e.g. design weights)
#' @param N population size for calibration of quantiles
#' @param pop_totals a vector of population totals for \code{X}
#' @param pop_quantiles a vector of population quantiles for \code{X_q}
#' @param backend specify an R package to perform calibration
#' @param method specify distance function for calibration
#' @param ... arguments passed either to \code{sampling::calib} or \code{laeken::calibWeights}
#'
#' @references
#'
#' Deville, J. C., and Särndal, C. E. (1992). Calibration estimators in survey sampling.
#' Journal of the American statistical Association, 87(418), 376-382.
#'
#' Harms, T. and Duchesne, P. (2006). On calibration estimation for quantiles.
#' Survey Methodology, 32(1), 37.
#'
#' @returns Returns a list with containing:\cr
#' \itemize{
#' \item{\code{w} -- final weight}
#' \item{\code{Xs} -- matrix used for calibration (i.e. Intercept, X and X_q transformed for calibration of quantiles)}
#' \item{\code{totals} -- a vector of totals (i.e. \code{N}, \code{pop_totals} and quantiles)}
#' \item{\code{diff} -- difference between \code{colSums(Xs*w)} and \code{totals}}
#' }
#'
#' @examples
#' \dontrun{
#'
#' }
joint_calibrate <-
function(X = NULL,
         X_q = NULL,
         d = NULL,
         N = NULL,
         pop_totals = NULL,
         pop_quantiles = NULL,
         backend = c("sampling", "laeken"),
         method = c("raking", "linear", "logit"),
         ...) {

  stopifnot("X and pop_totals have different dimensions" = ncol(X) == NROW(pop_totals))
  stopifnot("X_q and pop_quantiles have different dimensions"= ncol(X_q) == length(pop_quantiles))
  stopifnot("At least one element of pop_quantiles is empty (length of 0)" = all(lengths(pop_quantiles) > 0))
  stopifnot("Ony `sampling` and `laeken` are possible backends" = backend %in% c("sampling", "laeken"))
  stopifnot("Ony `raking`, `linear` and `logit` are possible" = method %in% c("linear", "raking", "logit"))


  ## quantiles
  totals_q_vec <- unlist(pop_quantiles)
  quantiles <- names(totals_q_vec)
  if (all(grepl("%", quantiles))) {
    quantiles <- as.numeric(gsub("%", "", quantiles))/100
  }
  ## pop_quantiles
  T_mat <- c(N, quantiles, pop_totals)
  A <- joint_calib_create_matrix(X_q, N, pop_quantiles)
  X <- cbind(1, A, X)

  if (backend == "sampling") {
    w_res <- sampling::calib(Xs = X, d = d, total = T_mat, method = method, ...)
  }
  if (backend == "laeken") {
    w_res <- laeken::calibWeights(X=X, d= d, totals = T_mat, method = method, ...)
  }
  w <- w_res*d
  return(list(w=w, Xs = X, totals = T_mat, diff = colSums(X*w) - T_mat))
}
