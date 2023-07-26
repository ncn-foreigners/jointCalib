library(sampling)
library(laeken)
library(survey)
library(data.table)
library(ggplot2)
library(scales)
library(Rcpp)
library(nonprobsvy)
library(xtable)
source("codes/functions.R") ## to be replaced by package

sourceCpp("codes/syspps.cpp")


seed_for_sim <- 2023-07-14
set.seed(seed_for_sim+1)
n_reps <- 1000


# generate data according to the paper ------------------------------------

N <- 20000
n_1 <- 1000
n_2 <- 1000
z1 <- rbinom(N, 1, 0.5)
z2 <- runif(N,0,2)
z3 <- rexp(N,1)
z4 <- rchisq(N, 4)
x1 <- z1
x2 <- z2 + 0.3*z1
x3 <- z3 + 0.2*(x1 + x2)
x4 <- z4 + 0.1*(x1 + x2 + x3)
e <- rnorm(N)

etol <- 1e-8
sigma03 <- stats::uniroot(f = function(s) cor(2+x1+x2+x3+x4, 2+x1+x2+x3+x4+s*e) - 0.3, c(0,20), tol = etol)$root
sigma05 <- stats::uniroot(f = function(s) cor(2+x1+x2+x3+x4, 2+x1+x2+x3+x4+s*e) - 0.5, c(0,20), tol = etol)$root
sigma08 <- stats::uniroot(f = function(s) cor(2+x1+x2+x3+x4, 2+x1+x2+x3+x4+s*e) - 0.8, c(0,20), tol = etol)$root

y1 <- 2+x1+x2+x3+x4+sigma03*e
y2 <- 2+x1+x2+x3+x4+sigma05*e
y3 <- 2+x1+x2+x3+x4+sigma08*e

theta_A1 <- stats::uniroot(f = function(s) sum(plogis(s + 0.1*x1+0.2*x2+0.3*x3+0.2*x4))-n_1, c(-20,0), tol = etol)$root

pA_1 <- plogis(theta_A1 + 0.1*x1+0.2*x2+0.3*x3+0.2*x4)
w_pA_1 <- 1/pA_1

s_B <- stats::uniroot(f = function(s) max(s+x3)/min(s+x3) - 50, c(0, 100), tol = etol)$root
pB_1 <- inclusionprobabilities(s_B + x3, n_1)
w_pB_1 <- 1/pB_1

pop_data <- data.frame(x1,x2,x3,x4,y1,y2,y3,pA_1,w_pA_1,pB_1,w_pB_1) |> setDT()


# known population quantities ---------------------------------------------

x_totals <- with(pop_data, c(sum(x1), sum(x2), sum(x3), sum(x4)))
p_quantiles_est <- seq(0.25, 0.75, 0.25) ## for estimation
#p_quantiles <- seq(0.1, 0.9, 0.1) ## for calibration
p_quantiles <- p_quantiles_est ## for calibration
x2_q <- with(pop_data, quantile(x2, p_quantiles))
x3_q <- with(pop_data, quantile(x3, p_quantiles))
x4_q <- with(pop_data, quantile(x4, p_quantiles))
k_quants <- NROW(p_quantiles)-1
k_quants_est <- NROW(p_quantiles_est)-1


# known population y ------------------------------------------------------

pop_true_vals <- pop_data[, lapply(.SD, 
                                   function(x) data.frame(m=mean(x), 
                                                          q25=quantile(x, 0.25), 
                                                          q50=quantile(x, 0.5), 
                                                          q75=quantile(x, 0.75))), 
                          .SDcols = c("y1", "y2", "y3")] |> 
  melt(value.name = "true") |> 
  {\(x) x[, c("y", "estimator") := tstrsplit(variable, "\\.")][, .(y, estimator, true)]}()

pop_true_vals[, estimator2:= factor(estimator, c("m", "q25", "q50", "q75"), c("Mean", "Q25", "Q50", "Q75"))]



# simulation --------------------------------------------------------------


results_y3_sample <- 
  results_y2_sample <- 
  results_y1_sample <- matrix(data=0, 
                              nrow = n_reps, 
                              ncol = (1+(k_quants_est+1))*4 + 1 + 3 + (k_quants_est+1)*3)



for (r in 1:n_reps) {
  set.seed(r)
  if (r %% 50 == 0) cat("Iteration: ", r, "\n")
  
  ## nonprob sample
  sample_nonprob <- pop_data[which(sampling::UPpoisson(pop_data$pA_1)==1),]
  sample_nonprob$w_naive <- N/n_1
  
  ## probability sample
  sample_prob <- pop_data[syspps_cpp(pop_data$pB_1, n_1), ] 
  sample_prob_svy <- svydesign(ids=~1, probs =~pB_1, data = sample_prob)
  
  ## estimated population size, totals and quantiles
  N_hat <- sum(weights(sample_prob_svy))
  x_totals <- svytotal( ~ x1 + x2 + x3 + x4, sample_prob_svy) |> as.numeric()
  q_est <- svyquantile( ~ x2 + x3 + x4, sample_prob_svy, p_quantiles)
  
  ## IPW
  ipw_y1 <- nonprob(selection = ~ x1 + x2 + x3 + x4,
                    target = ~ y1,
                    svydesign = sample_prob_svy,
                    data = sample_nonprob)
  
  ipw_y2 <- nonprob(selection = ~ x1 + x2 + x3 + x4,
                    target = ~ y2,
                    svydesign = sample_prob_svy,
                    data = sample_nonprob)
  
  ipw_y3 <- nonprob(selection = ~ x1 + x2 + x3 + x4,
                    target = ~ y3,
                    svydesign = sample_prob_svy,
                    data = sample_nonprob)
  
  ## MI based on model
  mi_y1 <- nonprob(outcome = y1 ~ x1 + x2 + x3 + x4,
                   svydesign = sample_prob_svy,
                   data = sample_nonprob)
  
  mi_y2 <- nonprob(outcome = y2 ~ x1 + x2 + x3 + x4,
                   svydesign = sample_prob_svy,
                   data = sample_nonprob)
  
  mi_y3 <- nonprob(outcome = y3 ~ x1 + x2 + x3 + x4,
                   svydesign = sample_prob_svy,
                   data = sample_nonprob)

  
  ## DR
  dr_y1 <- nonprob(selection = ~ x1 + x2 + x3 + x4,
                   outcome = y1 ~ x1 + x2 + x3 + x4,
                   svydesign = sample_prob_svy,
                   data = sample_nonprob)
  
  dr_y2 <- nonprob(selection = ~ x1 + x2 + x3 + x4,
                   outcome = y2 ~ x1 + x2 + x3 + x4,
                   svydesign = sample_prob_svy,
                   data = sample_nonprob)
  
  dr_y3 <- nonprob(selection = ~ x1 + x2 + x3 + x4,
                   outcome = y3 ~ x1 + x2 + x3 + x4,
                   svydesign = sample_prob_svy,
                   data = sample_nonprob)
  
  ## calibration with quantiles only
  w_res <- calib_quantiles(X_q = with(sample_nonprob, cbind(x2, x3, x4)),  
                           d = sample_nonprob$w_naive,  
                           N=N_hat,  
                           totals_q = list(q_est$x2[, 1], q_est$x3[, 1], q_est$x4[, 1]), 
                           method = "raking",
                           backend = "sampling")
  
  ## calibration with quantiles and totals
  w_res2 <- calib_quantiles(X_q = with(sample_nonprob, cbind(x2, x3, x4)),
                            X =  with(sample_nonprob, cbind(x1, x2, x3, x4)),
                            d = sample_nonprob$w_naive,
                            N=N_hat,
                            totals = x_totals,
                            totals_q = list(q_est$x2[, 1], q_est$x3[, 1], q_est$x4[, 1]),
                            method = "raking",
                            backend = "sampling")
  
  ## calibration with totals only
  w_res_st <- calib(Xs = with(sample_nonprob, cbind(x1, x2, x3, x4)), 
                    d = sample_nonprob$w_naive, 
                    total = x_totals, 
                    method = "raking")
  
  ## weights
  sample_nonprob$w_cal_q1 <- w_res$w
  sample_nonprob$w_cal_q2 <- w_res2$w
  sample_nonprob$w_cal_s <- sample_nonprob$w_naive*w_res_st
  sample_nonprob$w_ipw <- ipw_y1$weights
  sample_nonprob$w_dr <- dr_y1$weights
  
  results_y1_sample[r, 1] <- r
  results_y1_sample[r, 2] <- with(sample_nonprob, mean(y1))
  results_y1_sample[r, 3] <- with(sample_nonprob, weighted.mean(y1, w_cal_s))
  results_y1_sample[r, 4] <- with(sample_nonprob, weighted.mean(y1, w_cal_q1))
  results_y1_sample[r, 5] <- with(sample_nonprob, weighted.mean(y1, w_cal_q2))
  results_y1_sample[r, 6] <- ipw_y1$output$mean
  results_y1_sample[r, 7] <- mi_y1$output$mean
  results_y1_sample[r, 8] <- dr_y1$output$mean
  results_y1_sample[r, 9:(9+k_quants_est)] <- with(sample_nonprob, quantile(y1, p_quantiles_est))
  results_y1_sample[r, (9+1*k_quants_est+1):(9+2*k_quants_est+1)] <- with(sample_nonprob, weightedQuantile(y1, w_ipw,p_quantiles_est))
  results_y1_sample[r, (9+2*k_quants_est+2):(9+3*k_quants_est+2)] <- quantile(mi_y1$outcome$glm$fitted.values, p_quantiles_est)
  results_y1_sample[r, (9+3*k_quants_est+3):(9+4*k_quants_est+3)] <- with(sample_nonprob, weightedQuantile(y1, w_dr,p_quantiles_est))
  results_y1_sample[r, (9+4*k_quants_est+4):(9+5*k_quants_est+4)] <- with(sample_nonprob, weightedQuantile(y1, w_cal_s, p_quantiles_est))
  results_y1_sample[r, (9+5*k_quants_est+5):(9+6*k_quants_est+5)] <- with(sample_nonprob, weightedQuantile(y1, w_cal_q1, p_quantiles_est))
  results_y1_sample[r, (9+6*k_quants_est+6):(9+7*k_quants_est+6)] <- with(sample_nonprob, weightedQuantile(y1, w_cal_q2, p_quantiles_est))
  
  results_y2_sample[r, 1] <- r
  results_y2_sample[r, 2] <- with(sample_nonprob, mean(y2))
  results_y2_sample[r, 3] <- with(sample_nonprob, weighted.mean(y2, w_cal_s))
  results_y2_sample[r, 4] <- with(sample_nonprob, weighted.mean(y2, w_cal_q1))
  results_y2_sample[r, 5] <- with(sample_nonprob, weighted.mean(y2, w_cal_q2))
  results_y2_sample[r, 6] <- ipw_y2$output$mean
  results_y2_sample[r, 7] <- mi_y2$output$mean
  results_y2_sample[r, 8] <- dr_y2$output$mean
  results_y2_sample[r, 9:(9+k_quants_est)] <- with(sample_nonprob, quantile(y2, p_quantiles_est))
  results_y2_sample[r, (9+1*k_quants_est+1):(9+2*k_quants_est+1)] <- with(sample_nonprob, weightedQuantile(y2, w_ipw,p_quantiles_est))
  results_y2_sample[r, (9+2*k_quants_est+2):(9+3*k_quants_est+2)] <- quantile(mi_y2$outcome$glm$fitted.values, p_quantiles_est)
  results_y2_sample[r, (9+3*k_quants_est+3):(9+4*k_quants_est+3)] <- with(sample_nonprob, weightedQuantile(y2, w_dr,p_quantiles_est))
  results_y2_sample[r, (9+4*k_quants_est+4):(9+5*k_quants_est+4)] <- with(sample_nonprob, weightedQuantile(y2, w_cal_s, p_quantiles_est))
  results_y2_sample[r, (9+5*k_quants_est+5):(9+6*k_quants_est+5)] <- with(sample_nonprob, weightedQuantile(y2, w_cal_q1, p_quantiles_est))
  results_y2_sample[r, (9+6*k_quants_est+6):(9+7*k_quants_est+6)] <- with(sample_nonprob, weightedQuantile(y2, w_cal_q2, p_quantiles_est))
  
  results_y3_sample[r, 1] <- r
  results_y3_sample[r, 2] <- with(sample_nonprob, mean(y3))
  results_y3_sample[r, 3] <- with(sample_nonprob, weighted.mean(y3, w_cal_s))
  results_y3_sample[r, 4] <- with(sample_nonprob, weighted.mean(y3, w_cal_q1))
  results_y3_sample[r, 5] <- with(sample_nonprob, weighted.mean(y3, w_cal_q2))
  results_y3_sample[r, 6] <- ipw_y3$output$mean
  results_y3_sample[r, 7] <- mi_y3$output$mean
  results_y3_sample[r, 8] <- dr_y3$output$mean
  results_y3_sample[r, 9:(9+k_quants_est)] <- with(sample_nonprob, quantile(y3, p_quantiles_est))
  results_y3_sample[r, (9+1*k_quants_est+1):(9+2*k_quants_est+1)] <- with(sample_nonprob, weightedQuantile(y3, w_ipw,p_quantiles_est))
  results_y3_sample[r, (9+2*k_quants_est+2):(9+3*k_quants_est+2)] <- quantile(mi_y3$outcome$glm$fitted.values, p_quantiles_est)
  results_y3_sample[r, (9+3*k_quants_est+3):(9+4*k_quants_est+3)] <- with(sample_nonprob, weightedQuantile(y3, w_dr,p_quantiles_est))
  results_y3_sample[r, (9+4*k_quants_est+4):(9+5*k_quants_est+4)] <- with(sample_nonprob, weightedQuantile(y3, w_cal_s, p_quantiles_est))
  results_y3_sample[r, (9+5*k_quants_est+5):(9+6*k_quants_est+5)] <- with(sample_nonprob, weightedQuantile(y3, w_cal_q1, p_quantiles_est))
  results_y3_sample[r, (9+6*k_quants_est+6):(9+7*k_quants_est+6)] <- with(sample_nonprob, weightedQuantile(y3, w_cal_q2, p_quantiles_est))
  
}


results_sample <- rbind(results_y1_sample, results_y2_sample, results_y3_sample) |> as.data.frame() |> setDT()

setnames(results_sample, names(results_sample), 
         c("iter", "m_naive", "m_cal", "m_qcal1", "m_qcal2",
           "m_ipw", "m_mi", "m_dr",
           paste(paste0("q", rep(c(25, 50, 75), times = 4)), 
                 rep(c("naive", "ipw", "mi", "dr", "cal", "qcal1", "qcal2"), each = 3), sep="_")))

results_sample[, y:=rep(c("y1", "y2", "y3"), each = n_reps)]
results_sample <- melt(results_sample, id.vars = c("iter", "y"), value.name = "value", variable.name = "estimator") 
results_sample[, c("estimator", "type") := tstrsplit(estimator, "_")]
results_sample <- results_sample[pop_true_vals, on = c("y", "estimator")]
results_sample[, type:= factor(type, c("naive", "ipw", "mi", "dr", "cal", "qcal1", "qcal2"),
                               c("Naive", "IPW", "MI", "DR", "Cal", "QCal1", "QCal2"))]
results_sample[, estimator:= factor(estimator, c("m", "q25", "q50", "q75"), c("Mean", "Q25", "Q50", "Q75"))]


saveRDS(results_sample, file = "results/sim-chen2020-jasa.rds")
fwrite(x = results_sample, file = "results/sim-chen2020-jasa.csv.gz")

