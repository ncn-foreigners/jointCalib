#' @title Function for the covariate balancing propensity score using \code{CBPS}
#' @author Maciej Beręsewicz
#'
#' @description
#' \code{joint_calib_cbps} allows quantile or mean and quantile balancing using covariate balancing propensity score via [CBPS::CBPS()] function.
#'
#' @param formula_means a formula with variables to be balanced at means,
#' @param formula_quantiles a formula with variables to be balanced at quantiles,
#' @param treatment a formula with treatment indicator,
#' @param data a data.frame with variables,
#' @param probs a vector or a named list with quantiles to be balanced (default is \code{c(0.25, 0.5, 0.75)}),
#' @param control a control list of parameters for creation of X_q matrix,
#' @param standardize Default is TRUE, which normalizes weights to sum to 1 within each treatment group,
#' @param method method passed to \code{CBPS} function,
#' @param ... other parameters passed to \code{CBPS} function.
#'
#' @references
#'
#' Imai, K., & Ratkovic, M. (2014). Covariate balancing propensity score.
#' Journal of the Royal Statistical Society Series B: Statistical Methodology, 76(1), 243-263.
#'
#' Fong C, Ratkovic M, Imai K (2022). CBPS: Covariate Balancing Propensity Score.
#' R package version 0.23, <https://CRAN.R-project.org/package=CBPS>.
#'
#' @returns Returns an \code{CBPS} object.
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
#' dat <- data.frame(D = D, X1 = X1, X2 = X2, X3 = X3, X1X3 = X1X3, Y = y)
#' head(dat)
#'
#' result <- joint_calib_cbps(formula_means = ~ X1 + X2 + X3,
#'                            formula_quantiles = ~ X1 + X2,
#'                            treatment = ~ D,
#'                            data = dat)
#'
#' result
#'
#' # calculate ATE
#' w_1 <- dat$D/fitted(result)
#' w_1 <- w_1/mean(w_1)
#' w_0 <- (1-dat$D)/(1-fitted(result))
#' w_0 <- w_0/mean(w_0)
#' mean((w_1-w_0)*dat$Y)
#'
#' # compare with standard CBPS
#'
#' result2 <- CBPS::CBPS(D ~ X1 + X2 + X3, data = dat, method = "exact", standardize = FALSE, ATT = 0)
#'
#' w_1a <- dat$D/fitted(result2)
#' w_1a <- w_1a/mean(w_1a)
#' w_0a <- (1-dat$D)/(1-fitted(result2))
#' w_0a <- w_0a/mean(w_0a)
#' mean((w_1a-w_0a)*dat$Y)
#'
#' @export
joint_calib_cbps <-
  function(formula_means = NULL,
           formula_quantiles = NULL,
           treatment = NULL,
           data,
           probs = c(0.25, 0.5, 0.75),
           control = control_calib(),
           standardize = FALSE,
           method = "exact",
           ...) {

    ## checks
    stopifnot("`formula_quantiles` should be provided" = !is.null(formula_quantiles))
    stopifnot("`probs` should be a vector or a named list" = inherits(probs, "numeric") | inherits(probs, "list"))

    ## split data by treatment
    treat <- stats::model.frame(treatment, data)
    id_treatment <- which(treat == 1)
    id_control <- which(treat == 0)

    ## calculate totals
    if (is.null(formula_means)) {
      X <- NULL
    } else {
      X <- stats::model.matrix(formula_means, data)
      X <- as.data.frame(X[, colnames(X) != "(Intercept)", drop = FALSE])
    }

    ## calculate quantiles
    X_q <- stats::model.matrix(formula_quantiles, data)
    X_q <- as.data.frame(X_q[, colnames(X_q) != "(Intercept)", drop = FALSE])

    ## check quantiles
    ## list of quantiles
    if (inherits(probs, "list")) {
      pop_quantiles <- lapply(colnames(X_q), FUN = function(x) stats::quantile(X_q[id_treatment, x], probs = probs[[x]]))
      names(pop_quantiles) <- colnames(X_q)
    } else {
      pop_quantiles <- lapply(X_q[id_treatment, ], stats::quantile, probs = probs)
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

    A_treat <- joint_calib_create_matrix(X_q = X_q[id_treatment,], N = NROW(X_q[id_treatment, ]),
                                         pop_quantiles = pop_quantiles, control = control)
    A_control <- joint_calib_create_matrix(X_q = X_q[id_control,], N = NROW(X_q[id_treatment, ]),
                                           pop_quantiles = pop_quantiles, control = control)

    A <- matrix(0, nrow = NROW(X_q), ncol = NCOL(A_control))
    A[id_treatment, ] <- A_treat
    A[id_control, ] <- A_control

    colnames(A) <- gsub("%", "", names(unlist(pop_quantiles)))

    Xs <- as.data.frame(cbind(treat=treat[,1], X, A))

    result <- CBPS::CBPS(formula = treat ~ ., data = Xs,
                         ATT = 0,
                         standardize = standardize,
                         method = method, ...)

    return(result)
}
