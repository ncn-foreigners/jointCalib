#' @title Function to balance the covariate distributions of a control and treatment group using \code{joint_calib}
#' @author Maciej Beręsewicz
#'
#' @description
#'
#' \code{joint_calib_att} allows quantile or mean and quantile balancing of the covariate distributions of the control and treatment groups. It provides a user-friendly interface for specifying the variables and quantiles to be balanced.
#' \code{joint_calib_att} uses \code{joint_calib} function, so the user can apply different methods to find the weights that balance the control and treatment groups. For more details see [jointCalib::joint_calib()] and Beręsewicz and Szymkowiak (2023) working paper.
#'
#' @param formula_means a formula with variables to be balanced at means,
#' @param formula_quantiles a formula with variables to be balanced at quantiles,
#' @param treatment a formula with a treatment indicator,
#' @param data a data.frame with variables,
#' @param probs a vector or a named list of quantiles to be balanced (default is \code{c(0.25, 0.5, 0.75)}),
#' @param ... other parameters passed to \code{joint_calib} function.
#'
#' @references
#'
#' Beręsewicz, M. and Szymkowiak, M. (2023) A note on joint calibration estimators for totals and quantiles
#' Arxiv preprint <https://arxiv.org/abs/2308.13281>
#'
#' Greifer N (2023). WeightIt: Weighting for Covariate Balance in Observational Studies.
#' R package version 0.14.2, <https://CRAN.R-project.org/package=WeightIt>.
#'
#' Greifer N (2023). cobalt: Covariate Balance Tables and Plots.
#' R package version 4.5.1, <https://CRAN.R-project.org/package=cobalt>.
#'
#' Ho, D., Imai, K., King, G., and Stuart, E. A. (2011).
#' MatchIt: Nonparametric Preprocessing for Parametric Causal Inference.
#' Journal of Statistical Software, 42(8), 1–28. <https://doi.org/10.18637/jss.v042.i08>
#'
#' Xu, Y., and Yang, E. (2023). Hierarchically Regularized Entropy Balancing.
#' Political Analysis, 31(3), 457-464. <https://doi.org/10.1017/pan.2022.12>
#'
#' @returns Returns a list with containing:\cr
#' \itemize{
#' \item{\code{g} -- g-weight that sums up to treatment group size,}
#' \item{\code{Xs} -- matrix used for balancing (i.e. Intercept, X based on \code{formula_means} and X_q transformed for balancing of quantiles based on \code{formula_quantiles} and \code{probs}),}
#' \item{\code{totals} -- a vector of treatment reference size (\code{N}), means (\code{pop_totals}) and order of quantiles (based on \code{formula_quantiles} and \code{probs})}.
#' \item{\code{method} -- selected method,}
#' \item{\code{backend} -- selected backend.}
#' }
#'
#'
#' @examples
#'
#' ## generate data as in the hbal package
#' set.seed(123)
#' N <- 1500
#' X1 <- rnorm(N)
#' X2 <- rnorm(N)
#' X3 <- rbinom(N, size = 1, prob = .5)
#' X1X3 <- X1*X3
#' D_star <- 0.5*X1 + 0.3*X2 + 0.2*X1*X2 - 0.5*X1*X3 -1
#' D <- ifelse(D_star > rnorm(N), 1, 0)
#' y <- 0.5*D + X1 + X2 + X2*X3 + rnorm(N)
#' dat <- data.frame(D = D, X1 = X1, X2 = X2, X3 = X3, X1X3=X1X3, Y = y)
#' head(dat)
#'
#' ## Balancing means of X1, X2 and X3 and quartiles (0.25, 0.5, 0.75) of X1 and X2
#' ## sampling::raking is used
#' results <- joint_calib_att(
#' formula_means = ~ X1 + X2 + X3,
#' formula_quantiles = ~ X1 + X2,
#' treatment = ~ D,
#' data = dat,
#' method = "raking"
#' )
#'
#' ## Results are presented with summary statistics of balance weights (g-weights)
#' ## and information on the accuracy of reproducing reference treatment distributions
#' results
#'
#' ## An interaction between X1 and X2 is added to means
#' results2 <- joint_calib_att(
#' formula_means = ~ X1 + X2 + X3 + X1*X3,
#' formula_quantiles = ~ X1 + X2,
#' treatment = ~ D,
#' data = dat,
#' method = "raking"
#' )
#'
#' ## Results with interaction are presented below
#' results2
#'
#' ## As noted in the documentation, the probs argument can be a named list of different orders
#' ## In this example, we specify that X1 should be balanced at the mean,
#' ## while X2 should be balanced at Q1 and Q3
#' results3 <- joint_calib_att(
#' formula_means = ~ X1 + X2 + X3 + X1*X3,
#' formula_quantiles = ~ X1 + X2,
#' treatment = ~ D,
#' data = dat,
#' method = "raking",
#' probs = list(X1 = 0.5, X2 = c(0.25, 0.75))
#' )
#'
#' ## Results with different orders are presented below
#' results3
#'
#' ## Finally, we specify an order of quantile for the interaction
#' results4 <- joint_calib_att(
#' formula_means = ~ X1 + X2 + X3,
#' formula_quantiles = ~ X1 + X2 + X1:X3,
#' treatment = ~ D,
#' data = dat,
#' probs = list(X1=0.5, X2 = c(0.25, 0.5), `X1:X3` = 0.75),
#' method = "raking"
#' )
#'
#' ## Results with Q3 balancing for interaction are presented below
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
    ## take into account that treatment may be boolean
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

