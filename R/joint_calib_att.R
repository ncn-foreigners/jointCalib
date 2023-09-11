#' @title Function for balancing the control and treatment group using \code{joint_calib}
#' @author Maciej Beręsewicz
#'
#' @description
#' \code{joint_calib_att} allows quantile or mean and quantile balancing of control and treatment groups. It provides a user-friendly interface to specify the variables and quantiles to be balanced.
#'
#' @param formula_means a formula with variables to be balanced at means,
#' @param formula_quantiles a formula with variables to be balanced at quantiles,
#' @param treatment a formula with treatment indicator
#' @param data a data.frame with variables,
#' @param probs a vector or a named list with quantiles to be balanced (default is \code{c(0.25, 0.5, 0.75)}),
#' @param ... other parameters passed to \code{joint_calib} function.
#'
#' @references
#'
#' Greifer N (2023). WeightIt: Weighting for Covariate Balance in Observational Studies.
#' R package version 0.14.2, <https://CRAN.R-project.org/package=WeightIt>.
#'
#' Greifer N (2023). cobalt: Covariate Balance Tables and Plots.
#' R package version 4.5.1, <https://CRAN.R-project.org/package=cobalt>.
#'
#' Ho, D., Imai, K., King, G., & Stuart, E. A. (2011).
#' MatchIt: Nonparametric Preprocessing for Parametric Causal Inference.
#' Journal of Statistical Software, 42(8), 1–28. <https://doi.org/10.18637/jss.v042.i08>
#'
#' Xu, Y., & Yang, E. (2023). Hierarchically Regularized Entropy Balancing.
#' Political Analysis, 31(3), 457-464. <https://doi.org/10.1017/pan.2022.12>
#'
#' @returns Returns a list with containing:\cr
#' \itemize{
#' \item{\code{g} -- g-weight that sums up to treatment group size,}
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
#' ## probs for interactions
#' results4 <- joint_calib_att(
#' formula_means = ~ X1 + X2 + X3,
#' formula_quantiles = ~ X1 + X2 + X1:X3,
#' treatment = ~ D,
#' data = dat,
#' probs = list(X1=0.5, X2 = c(0.25, 0.5), `X1:X3` = 0.75),
#' method = "raking"
#' )
#'
#' results4
#'
#' @export
joint_calib_att <-
  function(formula_means = NULL,
           formula_quantiles = NULL,
           treatment = NULL,
           data,
           probs = c(0.25, 0.5, 0.75),
           ...) {

    ## checks
    stopifnot("`formula_quantiles` should be provided" = !is.null(formula_quantiles))
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

    ## check for unique values
    check_uniques <- sapply(pop_quantiles, FUN=function(x) length(unique(x)))
    n_quantiles <- lengths(pop_quantiles)

    if (sum(n_quantiles-check_uniques) > 0) {
      locate_error <- abs(check_uniques-n_quantiles) > 0
      stop(
        paste0(
        "Non-unique values when calculating quantiles for the following variables: ",
        paste(names(which(locate_error)), collapse = ", "),
        ". Adjust `probs` for quantiles to avoid repeated reference values."
        )
      )
    }

    result <- joint_calib(formula_totals = formula_means,
                          formula_quantiles = formula_quantiles,
                          data = data_splitted[["0"]],
                          dweights = rep(1, unname(N_size["0"])),
                          N = unname(N_size["1"]),
                          pop_totals = pop_totals,
                          pop_quantiles = pop_quantiles,
                          control = control_calib(el_att = TRUE),
                          ...)

    return(result)
}

