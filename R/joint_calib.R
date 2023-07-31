#' @title Function to calibrate totals and quantiles jointly
#' @author Maciej Beręsewicz
#'
#' @description
#' \code{joint_calib} allows to specify matrix of X variables for calibration of totals (intercept should not be included) and matrix of X_q variables for calibration of quantiles.
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
#' Haziza, D., and Lesage, É. (2016). A discussion of weighting procedures for unit nonresponse. Journal of Official Statistics, 32(1), 129-145.
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
#' ## generate data based on Haziza and Lesage (2016)
#' set.seed(123)
#' N <- 1000
#' x <- runif(N, 0, 80)
#' y <- exp(-0.1 + 0.1*x) + rnorm(N, 0, 300)
#' p <- rbinom(N, 1, prob = exp(-0.2 - 0.014*x))
#' probs <- seq(0.1, 0.9, 0.1)
#' quants_known <- quantile(x, probs)
#' totals_known <- sum(x)
#' df <- data.frame(x, y, p)
#' df_resp <- df[df$p == 1, ]
#' df_resp$d <- N/nrow(df_resp)
#' y_quant_true <- quantile(y, probs)
#' ## standard calibration for comparison
#' result0 <- sampling::calib(Xs = cbind(1, df_resp$x),
#'                            d = df_resp$d,
#'                            total = c(N, totals_known),
#'                            method = "linear")
#' y_quant_hat0 <- laeken::weightedQuantile(x = df_resp$y, probs = probs, weights = result0*df_resp$d)
#' ## example 1: calibrate only quantiles (deciles)
#' result1 <- joint_calib(X_q = as.matrix(df_resp$x),
#'                       d = df_resp$d,
#'                       N = N,
#'                       pop_quantiles = list(quants_known),
#'                       method = "linear",
#'                       backend = "sampling")
#' ## estimate quantiles
#' y_quant_hat1 <- laeken::weightedQuantile(x = df_resp$y, probs = probs, weights = result1$w)
#'
#' ## compare with known
#' data.frame(standard = y_quant_hat0, est=y_quant_hat1, true=y_quant_true)
#'
#' ## example 2: calibrate with quantiles (deciles) and totals
#'
#' result2 <- joint_calib(X_q = as.matrix(df_resp$x),
#'                        X = as.matrix(df_resp$x),
#'                        d = df_resp$d,
#'                        N = N,
#'                        pop_quantiles = list(quants_known),
#'                        pop_totals = totals_known,
#'                        method = "linear",
#'                        backend = "sampling")
#' ## estimate quantiles
#' y_quant_hat2 <- laeken::weightedQuantile(x = df_resp$y, probs = probs, weights = result2$w)
#'
#' ## compare with known
#' data.frame(standard = y_quant_hat0, est1=y_quant_hat1, est2=y_quant_hat2, true=y_quant_true)
#' }
#'
#' @export
joint_calib <-
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
  stopifnot("X contains constant" = all(base::apply(X, 2, stats::sd) > 0))
  stopifnot("X_q contains constant" = all(base::apply(X_q, 2, stats::sd) > 0))
  stopifnot("Ony `sampling` and `laeken` are possible backends" = backend %in% c("sampling", "laeken"))
  stopifnot("Ony `raking`, `linear` and `logit` are possible" = method %in% c("linear", "raking", "logit"))

  ## processing quantiles
  totals_q_vec <- unlist(pop_quantiles)
  quantiles <- names(totals_q_vec)
  if (all(grepl("%", quantiles))) {
    quantiles <- as.numeric(gsub("%", "", quantiles))/100
  }
  ## create pop_totals
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
