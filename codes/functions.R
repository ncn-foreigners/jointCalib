calib_quantiles_create_matrix <- function(X_q, N, totals_q) {
  
  stopifnot("dim of X_q is not equal to length of totals_q"= NCOL(X_q) == length(totals_q))
  
  A <- list()
  for (k in 1:NCOL(X_q)) {
  
    totals_q_k  <-   totals_q[[k]]
    A_q <- matrix(0, nrow=NROW(X_q), ncol = NROW(totals_q[[k]]))
    
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
  return(A)
}

calib_quantiles <- function(X=NULL, X_q, d, N, totals=NULL, totals_q, 
                            backend = c("sampling", "laeken"), 
                            method = c("raking", "linear", "logit"), 
                            ...) {
  
  #stopifnot("If Z is provided then backend should be `sampling`" = NROW(Z) > 0 & backend == "sampling") 
  
  ## quantiles
  totals_q_vec <- unlist(totals_q)
  quantiles <- names(totals_q_vec)
  if (all(grepl("%", quantiles))) {
    quantiles <- as.numeric(gsub("%", "", quantiles))/100
  } else {
    quantiles <- as.numeric(quantiles)
  }
  ## totals_q
  if (is.null(X)) {
    T_mat <- c(N, quantiles)
    A <- calib_quantiles_create_matrix(X_q, N, totals_q)
    X <- cbind(1, A)
  } else {
    T_mat <- c(N, quantiles, totals)
    A <- calib_quantiles_create_matrix(X_q, N, totals_q)
    X <- cbind(1, A, X)  
  }
  
  
  if (backend == "sampling") {
    w_res <- sampling::calib(Xs = X, d = d, total = T_mat, method = method, ...)  
  }
  if (backend == "laeken") {
    w_res <- laeken::calibWeights(X=X, d= d, totals = T_mat, method = method, ...)
  }
  w <- w_res*d
  return(list(w=w, Xs = X, totals = T_mat, diff = colSums(X*w) - T_mat))
}