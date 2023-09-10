#' @title Function for the balancing control to treatment group using \code{joint_calib}
#' @author Maciej Beręsewicz
#'
#' @description
#' \code{joint_calib_att} allows balancing control to treatment group on quantiles or means and quantiles. It provides a user-friendly interface that allows to specify variables and quantiles to be balanced.
#'
#' @param formula_means a formula with variables to be balanced at means,
#' @param formula_quantiles a formula with variables to be balanced at quantiles,
#' @param treatment a formula with treatment indicator
#' @param data a data.frame with variables,
#' @param probs a vector or a named list with quantiles to be balanced,
#' @param ... other parameters passed to \code{joint_calib} function.
#'
#' @references reference
#'
#' @returns Returns a list with containing:\cr
#' \itemize{
#' \item{\code{g} -- g-weight that sums up to sample size,}
#' \item{\code{Xs} -- matrix used for calibration (i.e. Intercept, X and X_q transformed for calibration of quantiles),}
#' \item{\code{totals} -- a vector of totals (i.e. \code{N}, \code{pop_totals} and \code{pop_quantiles}),}
#' \item{\code{method} -- selected method,}
#' \item{\code{backend} -- selected backend.}
#' }
#'
#' @examples
#' # generate data as in hbal package
#' set.seed(123)
#' N <- 1500
#' X1 <- rnorm(N)
#' X2 <- rnorm(N)
#' X3 <- rbinom(N, size = 1, prob = .5)
#' X1X3 <- X1*X3
#' D_star <- 0.5*X1 + 0.3*X2 + 0.2*X1*X2 - 0.5*X1*X3 - 1
#' D <- ifelse(D_star > rnorm(N), 1, 0) # Treatment indicator
#' y <- 0.5*D + X1 + X2 + X2*X3 + rnorm(N) # Outcome
#' dat <- data.frame(D = D, X1 = X1, X2 = X2, X3 = X3, X1X3=X1X3, Y = y)
#' head(dat)
#'
#' results <- joint_calib_att(
#' formula_means = ~ X1 + X2 + X3,
#' formula_quantiles = ~ X1 + X2,
#' treatment = ~ D,
#' data = dat,
#' method = "raking"
#' )
#'
#' results
#'
#' ## interaction for means
#' results2 <- joint_calib_att(
#' formula_means = ~ X1 + X2 + X3 + X1*X3,
#' formula_quantiles = ~ X1 + X2,
#' treatment = ~ D,
#' data = dat,
#' method = "raking"
#' )
#'
#' results2
#'
#' ## probs as a named list with varying probs
#' results3 <- joint_calib_att(
#' formula_means = ~ X1 + X2 + X3 + X1*X3,
#' formula_quantiles = ~ X1 + X2,
#' treatment = ~ D,
#' data = dat,
#' method = "raking",
#' probs = list(X1 = 0.5, X2 = c(0.25, 0.75))
#' )
#'
#' results3
#'
#' @export
joint_calib_att <-
  function(formula_means = NULL,
           formula_quantiles = NULL,
           treatment = NULL,
           data,
           probs = c(0.25, 0.5, 0.75),
           ... ) {

    ## checks
    stopifnot("`probs` should be a vector or a named list" = inherits(probs, "numeric") | inherits(probs, "list"))
    ## split data by treatment
    data_splitted <- split(data, treatment)
    N_size <- sapply(data_splitted, nrow)

    ## calculate totals
    if (is.null(formula_means)) {
      pop_totals <- NULL
    } else {
      X <- stats::model.matrix(formula_means, data_splitted[["1"]])
      X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
      pop_totals <- colSums(X)
    }

    ## calculate quantiles
    X_q <- stats::model.matrix(formula_quantiles, data_splitted[["1"]])
    X_q <- as.data.frame(X_q[, colnames(X_q) != "(Intercept)", drop = FALSE])

    ## check quantiles
    ## list of quantiles
    if (inherits(probs, "list")) {
      pop_quantiles <- lapply(colnames(X_q), FUN = function(x) stats::quantile(X_q[, x], probs = probs[[x]]))
      names(pop_quantiles) <- colnames(X_q)
    } else {
      pop_quantiles <- lapply(X_q, stats::quantile, probs = probs)
    }


    result <- joint_calib(formula_totals = formula_means,
                          formula_quantiles = formula_quantiles,
                          data = data_splitted[["0"]],
                          dweights = rep(1, unname(N_size["0"])),
                          N = unname(N_size["1"]),
                          pop_totals = pop_totals,
                          pop_quantiles = pop_quantiles,
                          ...)

    return(result)
}

