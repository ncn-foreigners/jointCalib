calib_quantiles <-
function(X, X_q, d, N, totals, totals_q, 
                            backend = c("sampling", "laeken"), 
                            method = c("raking", "linear", "logit"), 
                            ...) {
  
  #stopifnot("If Z is provided then backend should be `sampling`" = NROW(Z) > 0 & backend == "sampling") 
  
  ## quantiles
  totals_q_vec <- unlist(totals_q)
  quantiles <- names(totals_q_vec)
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
