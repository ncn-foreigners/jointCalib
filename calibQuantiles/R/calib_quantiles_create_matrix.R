calib_quantiles_create_matrix <-
function(X_q, N, totals_q) {
  
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
