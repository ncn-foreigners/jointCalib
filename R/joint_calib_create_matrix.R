#' @import mathjaxr
NULL
#' @title An internal function to create an A matrix for calibration of quantiles
#' @author Maciej Beręsewicz
#'
#' \loadmathjax
#' @description
#' \code{joint_calib_create_matrix} is function that creates an \mjseqn{A = [a_{ij}]} matrix for calibration of quantiles. Function allows to create matrix using \code{logistic} interpolation (using \code{stats::plogis}, default) or \code{linear} (as in Harms and Duchesne (2006), i.e. slightly modified Heavyside function).
#'
#' In case of \code{logistic} interpolation elements of \mjseqn{A} are created as follows
#'
#' \mjsdeqn{a_{i j} = \frac{1}{(1 + \exp\left(-2l\left(x_{ij}-Q_{x_j, \alpha}\right)\right))N},}
#'
#' where \mjseqn{x_{ij}} is the \mjseqn{i}th row of the auxiliary variable \mjseqn{X_j}, \mjseqn{N} is the population size, \mjseqn{Q_{x_j, \alpha}} is the known population \mjseqn{\alpha}th quantile, and \mjseqn{l} is set to -1000 (by default).
#'
#' In case of  \code{linear} interpolation elements of \mjseqn{A} are created as follows
#'
#' \mjsdeqn{a_{i j}=
#' \begin{cases}
#' N^{-1}, &  x_{i j} \leqslant L_{x_{j}, r} \left(Q_{x_j, \alpha}\right), \cr
#' N^{-1} \beta_{x_{j}, r}\left(Q_{x_j, \alpha}\right), & x_{i j}=U_{x_{j}, r}\left(Q_{x_j, \alpha}\right), \cr
#' 0, & x_{i j}>U_{x_{j}, r} \left(Q_{x_j, \alpha}\right),
#' \end{cases}}
#'
#' \mjseqn{i=1,...,r}, \mjseqn{j=1,...,k}, where \mjseqn{r} is the set of respondents, \mjseqn{k} is the auxiliary variable index and
#'
#' \mjsdeqn{L_{x_{j}, r}(t) = \max \left\lbrace\left\lbrace{x_{i j}}, i \in s \mid x_{i j} \leqslant t\right\rbrace \cup \lbrace-\infty\rbrace \right\rbrace,}
#' \mjsdeqn{U_{x_{j}, r}(t) = \min \left\lbrace\left\lbrace{x_{i j}}, i \in s \mid x_{i j}>t\right\rbrace \cup \lbrace\infty\rbrace \right\rbrace,}
#' \mjsdeqn{\beta_{x_{j}, r}(t) = \frac{t-L_{x_{j}, s}(t)}{U_{x_{j}, s}(t)-L_{x_{j}, s}(t)},}
#'
#' \mjseqn{i=1,...,r}, \mjseqn{j=1,...,k}, \mjseqn{t \in \mathbb{R}}.
#'
#' @param X_q matrix of variables for calibration of quantiles,
#' @param N population size for calibration of quantiles,
#' @param pop_quantiles a vector of population quantiles for \code{X_q},
#' @param control a control parameter for creation of \code{X_q} matrix.
#'
#' @references
#'
#' Harms, T. and Duchesne, P. (2006). On calibration estimation for quantiles.
#' Survey Methodology, 32(1), 37.
#'
#' @return Return matrix A
#'
#' @examples
#' # Create matrix for one variable and 3 quantiles
#' set.seed(123)
#' N <- 1000
#' x <- as.matrix(rnorm(N))
#' quants <- list(quantile(x, c(0.25,0.5,0.75)))
#' A <- joint_calib_create_matrix(x, N, quants)
#' head(A)
#' colSums(A)
#'
#' # Create matrix with linear interpolation
#' A <- joint_calib_create_matrix(x, N, quants, control_calib(interpolation="linear"))
#' head(A)
#' colSums(A)
#'
#' # Create matrix for two variables and different number of quantiles
#'
#' set.seed(123)
#' x1 <- rnorm(N)
#' x2 <- rchisq(N, 1)
#' x <- cbind(x1, x2)
#' quants <- list(quantile(x1, 0.5), quantile(x2, c(0.1, 0.75, 0.9)))
#' B <- joint_calib_create_matrix(x, N, quants)
#' head(B)
#' colSums(B)
#' @export
joint_calib_create_matrix <-
function(X_q, N, pop_quantiles, control = control_calib()) {

  interpolation <-   control$interpolation
  logit_const <- control$logit_const

  A <- switch(interpolation,
              linear = {
                A <- list()
                for (k in 1:NCOL(X_q)) {
                  totals_q_k  <-   pop_quantiles[[k]]
                  A_q <- matrix(0, nrow=NROW(X_q), ncol = NROW(pop_quantiles[[k]]))
                  x_sorted <- sort(X_q[, k])
                  for (i in 1:NROW(totals_q_k)) {
                    poz <- which(x_sorted <= totals_q_k[i])
                    n_poz <- NROW(poz)
                    L <- x_sorted[n_poz]
                    U <- x_sorted[n_poz + 1]
                    B <- (totals_q_k[i] - L) / (U - L)
                    A_q[X_q[, k] < L, i] <- 1/N
                    A_q[X_q[, k] == U, i] <- B/N
                  }
                  A[[k]] <- A_q
                }
                A <- Reduce("cbind", A)
                A
              },
              logit = {
                A <- lapply(1:length(pop_quantiles), FUN = function(x) {
                  sapply(unname(pop_quantiles[[x]]), FUN = function(y) stats::plogis(logit_const*(X_q[, x]-y))/N)
                })

                A <- Reduce("cbind", A)
                A
              })

  return(A)
}
