#' @title Function to balance the covariate distributions using covariate balancing propensity score \code{CBPS}
#' @author Maciej Beręsewicz
#'
#' @description
#'
#' \code{joint_calib_cbps} allows quantile or mean and quantile balancing of the covariate distributions of the control and treatment groups using the covariate balancing propensity score method (Imai & Ratkovic (2014)). CBPS::CBPS()] and [CBPS::hdCBPS()] are used a backend for estimating the parameters.
#' This function works in a similar way to the [jointCalib::joint_calib_att()] function, i.e. the user can specify variables for the balancing means as well as the quantiles.
#'
#' @param formula_means a formula with variables to be balanced at means,
#' @param formula_quantiles a formula with variables to be balanced at quantiles,
#' @param treatment a formula with a treatment indicator,
#' @param data a data.frame with variables,
#' @param probs a vector or a named list of quantiles to be balanced (default is \code{c(0.25, 0.5, 0.75)}),
#' @param control a control list of parameters for creation of X_q matrix based on \code{formula_quantiles} and \code{probs} (see [jointCalib::joint_calib_create_matrix()]),
#' @param standardize default is FALSE, which normalizes weights to sum to 1 within each treatment group (passed to \code{CBPS()} function),
#' @param method default is "exact". Choose "over" to fit an over-identified model that combines the propensity score and covariate balancing conditions; choose "exact" to fit a model that only contains the covariate balancing conditions (passed to \code{CBPS()} function)
#' @param variable_selection default is FALSE. Set to TRUE to select high dimension CBPS via [CBPS::hdCBPS()],
#' @param target specify target (y) variable for \code{hdCBPS} function,
#' @param ... other parameters passed to \code{CBPS} or \code{hdCBPS} functions.
#'
#' @references
#'
#' Imai, K., and Ratkovic, M. (2014). Covariate balancing propensity score.
#' Journal of the Royal Statistical Society Series B: Statistical Methodology, 76(1), 243-263.
#'
#' Fong C, Ratkovic M, and Imai K (2022). CBPS: Covariate Balancing Propensity Score.
#' R package version 0.23, <https://CRAN.R-project.org/package=CBPS>.
#'
#' @returns Returns a \code{CBPS} or a \code{list} object as a result of the \code{hdCBPS} function.
#'
#' @examples
#'
#' ## generate data as in the hbal package (see [hbal::hbal()])
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
#' ## Balancing means of X1, X2 and X3 and quartiles (0.25, 0.5, 0.75) of X1 and X2.
#' result <- joint_calib_cbps(formula_means = ~ X1 + X2 + X3,
#'                            formula_quantiles = ~ X1 + X2,
#'                            treatment = ~ D,
#'                            data = dat)
#'
#' ## CBPS output is presented
#' result
#'
#' ## calculate ATE by hand
#' w_1 <- dat$D/fitted(result)
#' w_1 <- w_1/mean(w_1)
#' w_0 <- (1-dat$D)/(1-fitted(result))
#' w_0 <- w_0/mean(w_0)
#' mean((w_1-w_0)*dat$Y)
#'
#' ## Compare with standard CBPS using only means
#' result2 <- CBPS::CBPS(D ~ X1 + X2 + X3, data = dat, method = "exact", standardize = FALSE, ATT = 0)
#'
#' ## calculate ATE by hand
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
           variable_selection = FALSE,
           target = NULL,
           ...) {

    ## checks
    stopifnot("`formula_quantiles` should be provided" = !is.null(formula_quantiles))
    stopifnot("`probs` should be a vector or a named list" = inherits(probs, "numeric") | inherits(probs, "list"))
    stopifnot("If `variable_selection == TRUE` then `target` should be provided" = !(variable_selection == TRUE & is.null(target)))

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
    names(Xs)[1] <- names(treat)

    if (variable_selection) {
      cat("Variable selection via `CBPS::hdCBPS` started...")
      y <- stats::model.frame(target, data)
      result <- CBPS::hdCBPS(formula = stats::as.formula(paste(names(Xs)[1], "~ .")),
                             data = Xs,
                             y = y[, 1],
                             ATT = 0,
                             ...)
    } else {
      result <- CBPS::CBPS(formula = stats::as.formula(paste(names(Xs)[1], "~ .")),
                           data = Xs,
                           ATT = 0,
                           standardize = standardize,
                           method = method, ...)
    }

    return(result)
}
