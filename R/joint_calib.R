#' @title Function to calibrate totals and quantiles jointly
#' @author Maciej Beręsewicz
#'
#' @description
#' \code{joint_calib} allows to specify matrix of X variables for calibration of totals (intercept should not be included) and matrix of X_q variables for calibration of quantiles.
#'
#' @param formula_totals a formula with variables for calibration of totals
#' @param formula_quantiles a formula with variables for calibration of quantiles
#' @param data a data.frame with variables
#' @param design a svydesign object (currently not supported)
#' @param dweights initial d-weights for calibration (e.g. design-weights)
#' @param N population size for calibration of quantiles
#' @param pop_totals a named vector of population totals for \code{X}, should be provided exactly as in `survey` package (see `survey::calibrate`)
#' @param pop_quantiles a named list of of population quantiles for \code{X_q}
#' @param subset a formula for subset of data
#' @param backend specify an R package to perform calibration
#' @param method specify distance function for calibration
#' @param control a list of control parameters (currently not supported)
#' @param ... arguments passed either to \code{sampling::calib}, \code{laeken::calibWeights} or \code{survey::calibrate}
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
#' \item{\code{totals} -- a vector of totals (i.e. \code{N}, \code{pop_totals} and \code{pop_quantiles})}
#' \item{\code{diff} -- difference between \code{colSums(Xs*w)} and \code{totals}}
#' }
#'
#' @examples
#' \donttest{
#' ## generate data based on Haziza and Lesage (2016)
#' set.seed(123)
#' N <- 1000
#' x <- runif(N, 0, 80)
#' y <- exp(-0.1 + 0.1*x) + rnorm(N, 0, 300)
#' p <- rbinom(N, 1, prob = exp(-0.2 - 0.014*x))
#' probs <- seq(0.1, 0.9, 0.1)
#' quants_known <- list(x=quantile(x, probs))
#' totals_known <- c(x=sum(x))
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
#' result1 <- joint_calib(formula_quantiles = ~x,
#'                        data = df_resp,
#'                        dweights = df_resp$d,
#'                        N = N,
#'                        pop_quantiles = quants_known,
#'                        method = "linear",
#'                        backend = "sampling")
#' ## estimate quantiles
#' y_quant_hat1 <- laeken::weightedQuantile(x = df_resp$y, probs = probs, weights = result1$w)
#'
#' ## compare with known
#' data.frame(standard = y_quant_hat0, est=y_quant_hat1, true=y_quant_true)
#'
#' ## example 2: calibrate with quantiles (deciles) and totals
#'
#' result2 <- joint_calib(formula_totals = ~x,
#'                        formula_quantiles = ~x,
#'                        data = df_resp,
#'                        dweights = df_resp$d,
#'                        N = N,
#'                        pop_quantiles = quants_known,
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
#' @seealso
#' [sampling::calib()] -- for standard calibration.
#' [laeken::calibWeights()] -- for standard calibration.
#' [survey::calibrate()] -- for standard calibration.
#'
#' @export
joint_calib <-
function(formula_totals = NULL,
         formula_quantiles = NULL,
         data = NULL,
         design = NULL, ## TBA
         dweights = NULL,
         N = NULL,
         pop_totals = NULL,
         pop_quantiles = NULL,
         subset = NULL,
         control = NULL,
         backend = c("sampling", "laeken", "survey"),
         method = c("raking", "linear", "logit"),
         ...) {

  ## processing
  if (is.null(formula_quantiles)) {
    stop("Parameter formula_quantiles is required. If you would like to use standard calibration we suggest using survey, sampling or laeken package.")
  }

  if (missing(backend)) backend <- "sampling"
  if (missing(method)) method <- "linear"


  stopifnot("Ony `sampling` and `laeken` are possible backends" = backend %in% c("sampling", "laeken", "survey"))
  stopifnot("Ony `raking`, `linear` and `logit` are possible" = method %in% c("linear", "raking", "logit"))

  subset <- parse(text = deparse(substitute(subset)))

  if (!is.logical(subset)) {subset <- eval(subset, data)}
  if (is.null(subset)) {subset <- TRUE}

  data <- base::subset(data, subset = subset)

  if (!is.data.frame(data)) {
    data <- data.frame(data)
  }

  ## parse formulas

  if (!is.null(formula_totals)) {
    X <- stats::model.matrix(formula_totals, data)
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]

    stopifnot("X and pop_totals have different dimensions" = ncol(X) == NROW(pop_totals))
    stopifnot("X and pop_totals have different names"=all.equal(colnames(X), names(pop_totals)))
    stopifnot("X contains constant" = all(base::apply(X, 2, stats::sd) > 0))
  } else {
    X <- NULL
  }

   X_q <- stats::model.matrix(formula_quantiles, data)
   X_q <- X_q[, colnames(X_q) != "(Intercept)", drop = FALSE]
   stopifnot("pop_quantiles contains unnamed elements" = all(sapply(names(pop_quantiles), nchar) > 0))
   stopifnot("X_q and pop_quantiles have different dimensions" = ncol(X_q) == length(pop_quantiles))
   stopifnot("X_q and pop_quantiles have different names" = all.equal(colnames(X_q), names(pop_quantiles)))
   stopifnot("At least one element of pop_quantiles is empty (length of 0)" = all(lengths(pop_quantiles) > 0))
   stopifnot("X_q contains constant" = all(base::apply(X_q, 2, stats::sd) > 0))


  if (!is.null(dweights)) {
    dweights <- as.numeric(dweights)
  } else {
    dweights <- rep(1, nrow(X_q))
  }

  ## processing quantiles
  pop_quantiles <- pop_quantiles[colnames(X_q)]
  names(pop_quantiles) <- NULL
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
    w_res <- sampling::calib(Xs = X, d = dweights, total = T_mat, method = method, ...)
  }
  if (backend == "laeken") {
    w_res <- laeken::calibWeights(X=X, d= dweights, totals = T_mat, method = method, ...)
  }
  w <- w_res*dweights
  return(list(w=w, Xs = X, totals = T_mat, diff = colSums(X*w) - T_mat))
}
