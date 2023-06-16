calib_quantiles_create_matrix <- function(x, N, pop_quantiles) {
  
  stopifnot("one x allowed"=NCOL(x) == 1)
  
  A_q <- matrix(0, nrow=NROW(x), ncol = NROW(pop_quantiles))
  
  x_sorted <- sort(x)
  
  for (i in 1:NROW(pop_quantiles)) {
    poz <- which(x_sorted <= pop_quantiles[i])
    n_poz <- NROW(poz)
    L <- x_sorted[n_poz]
    U <- x_sorted[n_poz + 1]
    B <- (pop_quantiles[i] - L) / (U - L)
    A_q[x < L, i] <- 1/N
    A_q[x == U, i] <- B/N
  }
  
  return(A_q)
}

calib_quantiles <- function(X, X_q, d, N, totals, totals_q, 
                            backend = c("sampling", "laeken"), 
                            method = c("raking", "linear", "logit"), 
                            ...) {
  
  
  ## quantiles
  quantiles <- names(totals_q)
  if (all(grepl("%", quantiles))) {
    quantiles <- as.numeric(gsub("%", "", quantiles))/100
  }
  ## totals_q
  T_mat <- c(N, quantiles, totals)
  A <- calib_quantiles_create_matrix(X_q, N, totals_q)
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