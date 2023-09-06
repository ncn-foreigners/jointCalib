#' @title Function for the joint calibration of totals and quantiles
#' @author Maciej Beręsewicz
#'
#' @description
#' \code{joint_calib} allows joint calibration of totals and quantiles. It provides a user-friendly interface that includes the specification of variables in formula notation, a vector of population totals, a list of quantiles, and a variety of backends and methods.
#'
#' @param formula_totals a formula with variables to calibrate the totals,
#' @param formula_quantiles a formula with variables for quantile calibration,
#' @param data a data.frame with variables,
#' @param dweights initial d-weights for calibration (e.g. design weights),
#' @param N population size for calibration of quantiles,
#' @param pop_totals a named vector of population totals for \code{formula_totals}. Should be provided exactly as in `survey` package (see `survey::calibrate`),
#' @param pop_quantiles a named list of population quantiles for \code{formula_quantiles} or an \code{newsvyquantile} class object (from \code{survey::svyquantile} function),
#' @param subset a formula for subset of data,
#' @param backend specify an R package to perform the calibration. Only \code{sampling}, \code{laeken}, \code{survey}, \code{ebal} or \code{base} are allowed,
#' @param method specify method (i.e. distance function) for the calibration. Only \code{raking}, \code{linear}, \code{logit}, \code{sinh}, \code{truncated}, \code{el} (empirical likelihood), \code{eb} (entropy balancing) are allowed,
#' @param bounds a numeric vector of length two giving bounds for the g-weights,
#' @param maxit a numeric value representing the maximum number of iterations,
#' @param tol the desired accuracy for the iterative procedure (for \code{sampling}, \code{laeken}, \code{ebal}, \code{el}) or tolerance in matching population total for \code{survey::grake} (see help for [survey::grake])
#' @param eps the desired accuracy for computing the Moore-Penrose generalized inverse (see [MASS::ginv()])
#' @param control a list of control parameters (currently only for \code{joint_calib_create_matrix})
#' @param ... arguments passed either to \code{sampling::calib}, \code{laeken::calibWeights}, \code{survey::calibrate} or \code{optim::constrOptim}
#'
#' @references
#'
#' Beręsewicz,  M., and Szymkowiak, M. (2023). A note on joint calibration estimators for totals and quantiles
#' Arxiv preprint <https://arxiv.org/abs/2308.13281>
#'
#' Deville, J. C., and Särndal, C. E. (1992). Calibration estimators in survey sampling.
#' Journal of the American statistical Association, 87(418), 376-382.
#'
#' Harms, T. and Duchesne, P. (2006). On calibration estimation for quantiles.
#' Survey Methodology, 32(1), 37.
#'
#' Wu, C. (2005) Algorithms and R codes for the pseudo empirical likelihood method in survey sampling,
#' Survey Methodology, 31(2), 239.
#'
#' Zhang, S., Han, P., and Wu, C. (2023) Calibration Techniques Encompassing Survey Sampling,
#' Missing Data Analysis and Causal Inference, International Statistical Review 91, 165–192.
#'
#' Haziza, D., and Lesage, É. (2016). A discussion of weighting procedures for unit nonresponse.
#' Journal of Official Statistics, 32(1), 129-145.
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
#'
#' y_quant_hat0 <- laeken::weightedQuantile(x = df_resp$y,
#'                                          probs = probs,
#'                                          weights = result0*df_resp$d)
#' x_quant_hat0 <- laeken::weightedQuantile(x = df_resp$x,
#'                                          probs = probs,
#'                                          weights = result0*df_resp$d)
#'
#' ## example 1: calibrate only quantiles (deciles)
#' result1 <- joint_calib(formula_quantiles = ~x,
#'                        data = df_resp,
#'                        dweights = df_resp$d,
#'                        N = N,
#'                        pop_quantiles = quants_known,
#'                        method = "linear",
#'                        backend = "sampling")
#' ## estimate quantiles
#' y_quant_hat1 <- laeken::weightedQuantile(x = df_resp$y,
#'                                          probs = probs,
#'                                          weights = result1$g*df_resp$d)
#' x_quant_hat1 <- laeken::weightedQuantile(x = df_resp$x,
#'                                          probs = probs,
#'                                          weights = result1$g*df_resp$d)
#'
#' ## compare with known
#' data.frame(standard = y_quant_hat0, est=y_quant_hat1, true=y_quant_true)
#'
#' ## example 2: calibrate with quantiles (deciles) and totals
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
#' y_quant_hat2 <- laeken::weightedQuantile(x = df_resp$y,
#'                                          probs = probs,
#'                                          weights = result2$g*df_resp$d)
#' x_quant_hat2 <- laeken::weightedQuantile(x = df_resp$x,
#'                                          probs = probs,
#'                                          weights = result2$g*df_resp$d)
#'
#' ## compare with known
#' data.frame(standard = y_quant_hat0, est1=y_quant_hat1,
#'            est2=y_quant_hat2, true=y_quant_true)
#'
#' ## example 3: calibrate wigh quantiles (deciles) and totals with
#' ## hyperbolic sinus (sinh) and survey package
#'
#' result3 <- joint_calib(formula_totals = ~x,
#'                        formula_quantiles = ~x,
#'                        data = df_resp,
#'                        dweights = df_resp$d,
#'                        N = N,
#'                        pop_quantiles = quants_known,
#'                        pop_totals = totals_known,
#'                        method = "sinh",
#'                        backend = "survey")
#'
#' ## estimate quantiles
#' y_quant_hat3 <- laeken::weightedQuantile(x = df_resp$y,
#'                                          probs = probs,
#'                                          weights = result3$g*df_resp$d)
#' x_quant_hat3 <- laeken::weightedQuantile(x = df_resp$x,
#'                                          probs = probs,
#'                                          weights = result3$g*df_resp$d)
#'
#' ## example 4: calibrate wigh quantiles (deciles) and totals with ebal package
#' result4 <- joint_calib(formula_totals = ~x,
#'                        formula_quantiles = ~x,
#'                        data = df_resp,
#'                        dweights = df_resp$d,
#'                        N = N,
#'                        pop_quantiles = quants_known,
#'                        pop_totals = totals_known,
#'                        method = "eb",
#'                        backend = "ebal")
#'
#' ## estimate quantiles
#' y_quant_hat4 <- laeken::weightedQuantile(x = df_resp$y,
#'                                          probs = probs,
#'                                          weights = result4$g*df_resp$d)
#' x_quant_hat4 <- laeken::weightedQuantile(x = df_resp$x,
#'                                          probs = probs,
#'                                          weights = result4$g*df_resp$d)
#'
#' ## compare with known
#' data.frame(standard = y_quant_hat0,
#'            est1=y_quant_hat1,
#'            est2=y_quant_hat2,
#'            est3=y_quant_hat3,
#'            est4=y_quant_hat4,
#'            true=y_quant_true)
#' ## compare with known X
#' data.frame(standard = x_quant_hat0,
#'            est1=x_quant_hat1,
#'            est2=x_quant_hat2,
#'            est3=x_quant_hat3,
#'            est4=x_quant_hat4,
#'            true = quants_known$x)
#'
#'
#' @seealso
#' [sampling::calib()] -- for standard calibration.
#'
#' [laeken::calibWeights()] -- for standard calibration.
#'
#' [survey::calibrate()] -- for standard and more advanced calibration.
#'
#' [ebal::ebalance()] -- for standard entropy balancing.
#'
#'
#' @export
joint_calib <-
function(formula_totals = NULL,
         formula_quantiles = NULL,
         data = NULL,
         dweights = NULL,
         N = NULL,
         pop_totals = NULL,
         pop_quantiles = NULL,
         subset = NULL,
         backend = c("sampling", "laeken", "survey", "ebal", "base"),
         method = c("raking", "linear", "logit", "sinh", "truncated", "el", "eb"),
         bounds = c(0, 10),
         maxit = 50,
         tol = 1e-8,
         eps = .Machine$double.eps,
         control = control_calib(),
         ...) {

  ## processing
  if (is.null(formula_quantiles)) {
    stop("The `formula_quantiles` parameter is required. If you want to use
         standard calibration, we recommend using the `survey`, `sampling`, `laeken` or `ebal` packages.")
  }

  stopifnot("The `pop_quantiles` parameter should be of `list` or `newsvyquantile` class" =
              inherits(pop_quantiles, "newsvyquantile") | inherits(pop_quantiles, "list"))

  if (missing(backend)) backend <- "sampling"
  if (missing(method)) method <- "linear"

  if (method == "eb") backend <- "ebal"
  if (method == "el") backend <- "base"

  stopifnot("Only `survey`, `sampling`, `laeken`, `ebal` and `base` are possible backends" =
              backend %in% c("sampling", "laeken", "survey", "ebal", "base"))
  stopifnot("Only `raking`, `linear`, logit`, `sinh`, `truncated`, `el` and `eb` are possible" =
              method %in% c("linear", "raking", "logit", "sinh", "truncated", "el", "eb"))

  stopifnot("Method `sinh` is only possible with `survey`" = !(method == "sinh" & backend != "survey"))
  stopifnot("Method `truncated` is only possible with `sampling`" = !(method == "truncated" & backend != "sampling"))

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
    stopifnot("`X` and `pop_totals` have different dimensions" = ncol(X) == NROW(pop_totals))
    stopifnot("`X` and `pop_totals` have different names"=all.equal(colnames(X), names(pop_totals)))
    stopifnot("`X` contains constant" = all(base::apply(X, 2, stats::sd) > 0))
  } else {
    X <- NULL
  }

   X_q <- stats::model.matrix(formula_quantiles, data)
   X_q <- X_q[, colnames(X_q) != "(Intercept)", drop = FALSE]

   stopifnot("`pop_quantiles` contains unnamed elements" = all(sapply(names(pop_quantiles), nchar) > 0))
   stopifnot("`formula_quantiles` and `pop_quantiles` have different dimensions" = ncol(X_q) == length(pop_quantiles))
   stopifnot("`formula_quantiles` and `pop_quantiles` have different names" = all.equal(colnames(X_q), names(pop_quantiles)))
   stopifnot("At least one element of `pop_quantiles` is empty (length of 0)" = all(lengths(pop_quantiles) > 0))
   stopifnot("`formula_quantiles` contains constant" = all(base::apply(X_q, 2, stats::sd) > 0))


  if (!is.null(dweights)) {
    dweights <- as.numeric(dweights)
  } else {
    dweights <- rep(1, nrow(X_q))
  }

  ## processing quantiles [more robust is needed]
  if (inherits(pop_quantiles, "list")) {
    pop_quantiles <- pop_quantiles[colnames(X_q)]
    names(pop_quantiles) <- NULL
    totals_q_vec <- unlist(pop_quantiles)
    quantiles <- names(totals_q_vec)
    if (all(grepl("%", quantiles))) {
      quantiles <- as.numeric(gsub("%", "", quantiles))/100
    }
  }
  if (inherits(pop_quantiles, "newsvyquantile")) {
    pop_quantiles <- lapply(pop_quantiles, FUN=function(x) x[,1])
    pop_quantiles <- pop_quantiles[colnames(X_q)]
    names(pop_quantiles) <- NULL
    totals_q_vec <- unlist(pop_quantiles)
    quantiles <- names(totals_q_vec)
    quantiles <- as.numeric(quantiles)
  }

  ## create pop_totals -- scaled totals are needed
  #T_mat <- c(N/N, quantiles*N, pop_totals/N)
  T_mat <- c(N, quantiles, pop_totals)

  A <- joint_calib_create_matrix(X_q, N, pop_quantiles,
                                 control = control)
  #X <- cbind(1/N, A*N, X/N)
  X <- cbind(1, A, X)

  ## change to large switch
  gweights <- switch(backend,
                     "sampling" = {
                       if (method %in% c("logit", "truncated")) {
                         sampling::calib(Xs = X,
                                         d = dweights,
                                         total = T_mat,
                                         method = method,
                                         bounds = bounds,
                                         max_iter = maxit,
                                         ...)
                       } else {
                         sampling::calib(Xs = X,
                                         d = dweights,
                                         total = T_mat,
                                         method = method,
                                         max_iter = maxit,
                                         ...)
                       }
                       },
                     "laeken" = laeken::calibWeights(X = X,
                                                     d = dweights,
                                                     totals = T_mat,
                                                     method = method,
                                                     bounds = bounds,
                                                     maxit = maxit,
                                                     tol = tol,
                                                     eps = .Machine$double.eps,
                                                     ...),
                     "survey" = survey::grake(mm = X,
                                              ww = dweights,
                                              calfun = switch(method,
                                                              "linear" = survey::cal.linear,
                                                              "raking" = survey::cal.raking,
                                                              "logit" = survey::cal.logit,
                                                              "sinh" = survey::cal.sinh),
                                              population = T_mat,
                                              bounds = list(lower = bounds[1], upper = bounds[2]),
                                              epsilon = tol,
                                              maxit = maxit,
                                              verbose = FALSE,
                                              variance = NULL),
                     "ebal" =  ebal::eb(tr.total = T_mat,
                                        co.x = X,
                                        coefs = c(log(T_mat[1]/NROW(X)), rep(0, (NCOL(X) - 1))),
                                        base.weight = dweights,
                                        max.iterations = maxit,
                                        constraint.tolerance = control$ebal_constraint_tolerance,
                                        print.level = control$ebal_print_level
                                        )$Weights.ebal/dweights,
                     "base" = calib_el(X = X,
                                        d = dweights,
                                        totals = T_mat,
                                        tol = tol,
                                        maxit = maxit,
                                        eps = eps))

  gweights <- as.numeric(gweights)

  return(
    structure(
      list(g = gweights,
           Xs = X,
           totals = c(N, quantiles, pop_totals),
           diff = colSums(X * dweights * gweights) - T_mat,
           method = method,
           backend = backend),
      class = "jointCalib"
      )
  )
}

