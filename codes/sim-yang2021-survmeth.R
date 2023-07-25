library(sampling)
library(laeken)
library(survey)
library(data.table)
library(ggplot2)
library(scales)
library(Rcpp)
library(xtable)
library(glue)
library(nonprobsvy)
source("codes/functions.R") ## to be replaced by package


# generate data ------------------------------------------------------------

seed_for_sim <- 2023-07-14
set.seed(seed_for_sim+1)

n_reps <- 1000 ## number of simulations for the paper

N <- 100000
n <- 1000
x1 <- rnorm(N,1,1)
x2 <- rexp(N,1)
alp <- rnorm(N)
epsilon <- rnorm(N)
y11 <- 1 + x1 + x2 + alp + epsilon
y12 <- 0.5*(x1-1.5)^2 + x2^2 + alp + epsilon
y21 <- rbinom(N,1,plogis(1 + x1 + x2 + alp))
y22 <- rbinom(N,1,plogis(0.5*(x1-1.5)^2 + x2^2 + alp))
p1 <- plogis(x2)
p2 <- plogis(-3+(x1-1.5)^2+(x2-2)^2)
pop_data <- data.frame(x1,x2,y11,y12,y21,y22,p1,p2) |> setDT()
p_quantiles <- seq(0.25, 0.75, 0.25)
#p_quantiles <- seq(0.1, 0.9, 0.1)


# get true values ---------------------------------------------------------


pop_true_vals <- pop_data[, lapply(.SD, 
                                   function(x) data.frame(mean=mean(x), 
                                                          q25=quantile(x, 0.25), 
                                                          q50=quantile(x, 0.5), 
                                                          q75=quantile(x, 0.75))), 
                          .SDcols = c("y11", "y12", "y21", "y22")] |> 
  melt(value.name = "true") |> 
  {\(x) x[, c("y", "type") := tstrsplit(variable, "\\.")][, .(y, type, true)]}()


# main simulation ---------------------------------------------------------


results <- list()

for (r in 1:n_reps) {
  set.seed(r)
  if (r %% 10 == 0) cat("iteration: ", r, "\n")
  sample_prob <- pop_data[sample(1:N, n),]
  sample_prob$w <- N/n
  sample_prob_svy <- svydesign(ids=~1, weights = ~w, data = sample_prob)
  q_est <- svyquantile( ~ x1 + x2, sample_prob_svy, p_quantiles)
  x_totals <- svytotal( ~ x1 + x2, sample_prob_svy) |> as.numeric()
  
  sample_bd1 <- pop_data[rbinom(N,1,pop_data$p1)==1, ]
  sample_bd1$w_naive <- N/nrow(sample_bd1)
  sample_bd2 <- pop_data[rbinom(N,1,pop_data$p2)==1, ]
  sample_bd2$w_naive <- N/nrow(sample_bd2)
  
  ##########################################
  ## linear inclusion probability (BD1)
  ##########################################
  ## calibration with quantiles only (QCAL1)
  w_res <- calib_quantiles(X_q = with(sample_bd1, cbind(x1,x2)),
                           d = sample_bd1$w_naive,
                           N = N,
                           totals_q = list(q_est$x1[, 1], q_est$x2[, 1]),
                           method = "raking",
                           backend = "sampling")
  
  ## calibration with quantiles and totals (QCAL2)
  w_res2 <- calib_quantiles(X_q = with(sample_bd1, cbind(x1,x2)),
                            X =  with(sample_bd1, cbind(x1,x2)),
                            d = sample_bd1$w_naive,
                            N = N,
                            totals = x_totals,
                            totals_q = list(q_est$x1[, 1], q_est$x2[, 1]),
                            method = "raking",
                            backend = "sampling")
  ## calibration with totals only (CAL)
  w_res_st <- calib(Xs = with(sample_bd1, cbind(x1,x2)),
                    d = sample_bd1$w_naive,
                    total = x_totals,
                    method = "raking")
  ## ipw 
  ipw <- nonprob(selection = ~ x1 + x2,
                 target = ~ y11,
                 svydesign = sample_prob_svy,
                 data = sample_bd1)
  
  ## mi for y11
  mi_y11 <- nonprob(outcome = y11 ~ x1 + x2,
                    svydesign = sample_prob_svy,
                    data = sample_bd1,
                    method_outcome = "nn")
  
  ## mi for y12
  mi_y12 <- nonprob(outcome = y12 ~ x1 + x2,
                    svydesign = sample_prob_svy,
                    data = sample_bd1,
                    method_outcome = "nn")
  ## mi for y21
  mi_y21 <- nonprob(outcome = y21 ~ x1 + x2,
                    svydesign = sample_prob_svy,
                    data = sample_bd1,
                    method_outcome = "nn")
  ## mi for y22
  mi_y22 <- nonprob(outcome = y22 ~ x1 + x2,
                    svydesign = sample_prob_svy,
                    data = sample_bd1,
                    method_outcome = "nn")
  
  ## for quantile estimation based on mi
  mi_y11_vals <- apply(mi_y11$outcome$model_rand$nn.idx, 1, FUN=\(x) mean(sample_bd1$y11[x]))
  mi_y12_vals <- apply(mi_y12$outcome$model_rand$nn.idx, 1, FUN=\(x) mean(sample_bd1$y12[x]))
  
  ## dr for y11
  dr_y11 <- nonprob(selection = ~ x1 + x2,
                    outcome = y11 ~ x1 + x2,
                    svydesign = sample_prob_svy,
                    data = sample_bd1)
  ## dr for y12
  dr_y12 <- nonprob(selection = ~ x1 + x2,
                    outcome = y12 ~ x1 + x2, 
                    svydesign = sample_prob_svy,
                    data = sample_bd1)
  
  ## dr for y21
  dr_y21 <- nonprob(selection = ~ x1 + x2,
                    outcome = y21 ~ x1 + x2,
                    svydesign = sample_prob_svy,
                    data = sample_bd1,
                    family_outcome = "binomial")
  ## dr for y22
  dr_y22 <- nonprob(selection = ~ x1 + x2,
                    outcome = y22 ~ x1 + x2,
                    svydesign = sample_prob_svy,
                    data = sample_bd1,
                    family_outcome = "binomial")
  
  ############################################
  ## non-linear inclusion probability
  ############################################
  ## calibration with quantiles only (QCAL1)
  w_res_2 <- calib_quantiles(X_q = with(sample_bd2, cbind(x1, x2)),
                             d = sample_bd2$w_naive,
                             N = N,
                             totals_q = list(q_est$x1[, 1], q_est$x2[, 1]),
                             method = "raking",
                             backend = "sampling")
  
  ## calibration with quantiles only (QCAL2)
  w_res2_2 <- calib_quantiles(X_q = with(sample_bd2, cbind(x1, x2)),
                              X =  with(sample_bd2, cbind(x1, x2)),
                              d = sample_bd2$w_naive,
                              N = N,
                              totals = x_totals,
                              totals_q = list(q_est$x1[, 1], q_est$x2[, 1]),
                              method = "raking",
                              backend = "sampling")
  
  ## calibration with quantiles only (CAL)
  w_res_st_2 <- calib(Xs = with(sample_bd2, cbind(x1, x2)),
                      d = sample_bd2$w_naive,
                      total = x_totals,
                      method = "raking")
  
  ## ipw
  ipw_2 <- nonprob(selection = ~ x1 + x2,
                   target = ~ y11,
                   svydesign = sample_prob_svy,
                   data = sample_bd2)
  ## mi for y11
  mi_y11_2 <- nonprob(outcome = y11 ~ x1 + x2,
                      svydesign = sample_prob_svy,
                      data = sample_bd2,
                      method_outcome = "nn")
  ## mi for y12
  mi_y12_2 <- nonprob(outcome = y12 ~ x1 + x2,
                      svydesign = sample_prob_svy,
                      data = sample_bd2,
                      method_outcome = "nn")
  ## mi for y21
  mi_y21_2 <- nonprob(outcome = y21 ~ x1 + x2,
                      svydesign = sample_prob_svy,
                      data = sample_bd2,
                      method_outcome = "nn")
  ## mi for y22
  mi_y22_2 <- nonprob(outcome = y22 ~ x1 + x2,
                      svydesign = sample_prob_svy,
                      data = sample_bd2,
                      method_outcome = "nn")
  
  ## for quantiles based on mi
  mi_y11_vals_2 <- apply(mi_y11_2$outcome$model_rand$nn.idx, 1, FUN=\(x) mean(sample_bd2$y11[x]))
  mi_y12_vals_2 <- apply(mi_y12_2$outcome$model_rand$nn.idx, 1, FUN=\(x) mean(sample_bd2$y12[x]))
  
  ## dr for y11
  dr_y11_2 <- nonprob(selection = ~ x1 + x2,
                      outcome = y11 ~ x1 + x2,
                      svydesign = sample_prob_svy,
                      data = sample_bd2)
  ## dr for y12
  dr_y12_2 <- nonprob(selection = ~ x1 + x2,
                      outcome = y12 ~ x1 + x2,
                      svydesign = sample_prob_svy,
                      data = sample_bd2)
  ## dr for y21
  dr_y21_2 <- nonprob(selection = ~ x1 + x2,
                      outcome = y21 ~ x1 + x2,
                      svydesign = sample_prob_svy,
                      data = sample_bd2,
                      family_outcome = "binomial")
  ## dr for y22
  dr_y22_2 <- nonprob(selection = ~ x1 + x2,
                      outcome = y22 ~ x1 + x2,
                      svydesign = sample_prob_svy,
                      data = sample_bd2,
                      family_outcome = "binomial")
  ## add weights
  sample_bd1[, ":="(w_cal_q1 =w_res$w,  ## QCAL1 weights
                    w_cal_q2 =w_res2$w,  ## QCAL2 weights
                    w_cal_s = w_naive*w_res_st, ## CAL weights
                    w_ipw = ipw$weights,  ## IPW weighs
                    w_dr=dr_y11$weights)] ## DR (IPW) weights
  sample_bd2[, ":="(w_cal_q1 =w_res_2$w,  ## QCAL1 weights
                    w_cal_q2 = w_res2_2$w,  ## QCAL2 weights
                    w_cal_s = w_naive*w_res_st_2, ## CAL weights
                    w_ipw = ipw_2$weights,  ## IPW weighs
                    w_dr=dr_y11_2$weights)] ## DR (IPW) weights
  ################################################################
  ## summary of results
  ################################################################
  
  ## results for means
  res1 <- sample_bd1[, .(r=r, data="bd1", type = "mean", 
                         y11_naive=mean(y11), y11_cal=weighted.mean(y11, w_cal_s), 
                         y11_qcal1=weighted.mean(y11, w_cal_q1),  y11_qcal2=weighted.mean(y11, w_cal_q2),
                         y11_ipw=weighted.mean(y11, w_ipw), y11_dr=dr_y11$output$mean, y11_mi=mi_y11$output$mean, 
                         
                         y12_naive=mean(y12), y12_cal=weighted.mean(y12, w_cal_s),
                         y12_qcal1=weighted.mean(y12, w_cal_q1), y12_qcal2=weighted.mean(y12, w_cal_q2),
                         y12_ipw=weighted.mean(y12, w_ipw), y12_dr=dr_y12$output$mean, y12_mi=mi_y12$output$mean, 
                         
                         y21_naive=mean(y21), y21_cal=weighted.mean(y21, w_cal_s),
                         y21_qcal1=weighted.mean(y21, w_cal_q1), y21_qcal2=weighted.mean(y21, w_cal_q2),
                         y21_mi=mi_y21$output$mean, y21_ipw=weighted.mean(y21, w_ipw), y21_dr=dr_y21$output$mean,
                         
                         y22_naive=mean(y22), y22_cal=weighted.mean(y22, w_cal_s),
                         y22_qcal1=weighted.mean(y22, w_cal_q1), y22_qcal2=weighted.mean(y22, w_cal_q2),
                         y22_mi=mi_y22$output$mean,y22_ipw=weighted.mean(y22, w_ipw), y22_dr=dr_y22$output$mean)]
  
  
  res2 <- sample_bd2[, .(r=r, data="bd2", type = "mean",
                         y11_naive=mean(y11), y11_cal=weighted.mean(y11, w_cal_s),
                         y11_qcal1=weighted.mean(y11, w_cal_q1), y11_qcal2=weighted.mean(y11, w_cal_q2),
                         y11_ipw=weighted.mean(y11, w_ipw), y11_dr=dr_y11_2$output$mean, y11_mi=mi_y11_2$output$mean, 
                         
                         y12_naive=mean(y12), y12_cal=weighted.mean(y12, w_cal_s),
                         y12_qcal1=weighted.mean(y12, w_cal_q1), y12_qcal2=weighted.mean(y12, w_cal_q2),
                         y12_ipw=weighted.mean(y12, w_ipw), y12_dr=dr_y12_2$output$mean, y12_mi=mi_y12_2$output$mean, 
                         
                         y21_naive=mean(y21), y21_cal=weighted.mean(y21, w_cal_s),
                         y21_qcal1=weighted.mean(y21, w_cal_q1), y21_qcal2=weighted.mean(y21, w_cal_q2),
                         y21_mi=mi_y21_2$output$mean, y21_ipw=weighted.mean(y21, w_ipw), y21_dr=dr_y21_2$output$mean, 
                         
                         y22_naive=mean(y22), y22_cal=weighted.mean(y22, w_cal_s),
                         y22_qcal1=weighted.mean(y22, w_cal_q1), y22_qcal2=weighted.mean(y22, w_cal_q2),
                         y22_mi=mi_y22_2$output$mean, y22_ipw=weighted.mean(y22, w_ipw), y22_dr=dr_y22_2$output$mean)]
  
  ## results for quartiles
  res3 <- sample_bd1[, .(r=r, data="bd1", type = "q25",
                         y11_naive=quantile(y11,0.25), y11_cal=weightedQuantile(y11, w_cal_s,0.25),
                         y11_qcal1=weightedQuantile(y11, w_cal_q1,0.25), y11_qcal2=weightedQuantile(y11, w_cal_q2,0.25),
                         y11_ipw=weightedQuantile(y11, w_ipw,0.25), y11_dr=weightedQuantile(y11, w_dr,0.25),
                         y11_mi = quantile(mi_y11_vals, 0.25),
                         y12_naive=quantile(y12,0.25), y12_cal=weightedQuantile(y12, w_cal_s,0.25),
                         y12_qcal1=weightedQuantile(y12, w_cal_q1,0.25), y12_qcal2=weightedQuantile(y12, w_cal_q2,0.25),
                         y12_ipw=weightedQuantile(y12, w_ipw,0.25), y12_dr=weightedQuantile(y12, w_dr,0.25),
                         y12_mi = quantile(mi_y12_vals, 0.25))]
  
  res4 <- sample_bd2[, .(r=r, data="bd2", type = "q25",
                         y11_naive=quantile(y11,0.25), y11_cal=weightedQuantile(y11, w_cal_s,0.25),
                         y11_qcal1=weightedQuantile(y11, w_cal_q1,0.25), y11_qcal2=weightedQuantile(y11, w_cal_q2,0.25),
                         y11_ipw=weightedQuantile(y11, w_ipw,0.25), y11_dr=weightedQuantile(y11, w_dr,0.25),
                         y11_mi = quantile(mi_y11_vals_2, 0.25),
                         y12_naive=quantile(y12,0.25), y12_cal=weightedQuantile(y12, w_cal_s,0.25),
                         y12_qcal1=weightedQuantile(y12, w_cal_q1,0.25), y12_qcal2=weightedQuantile(y12, w_cal_q2,0.25),
                         y12_ipw=weightedQuantile(y12, w_ipw,0.25), y12_dr=weightedQuantile(y12, w_dr,0.25),
                         y12_mi = quantile(mi_y12_vals_2, 0.25))]
  
  res5 <- sample_bd1[, .(r=r, data="bd1", type = "q50",
                         y11_naive=quantile(y11,0.50), y11_cal=weightedQuantile(y11, w_cal_s,0.50),
                         y11_qcal1=weightedQuantile(y11, w_cal_q1,0.50), y11_qcal2=weightedQuantile(y11, w_cal_q2,0.50),
                         y11_ipw=weightedQuantile(y11, w_ipw,0.50), y11_dr=weightedQuantile(y11, w_dr,0.50),
                         y11_mi = quantile(mi_y11_vals, 0.50),
                         y12_naive=quantile(y12,0.50), y12_cal=weightedQuantile(y12, w_cal_s,0.50),
                         y12_qcal1=weightedQuantile(y12, w_cal_q1,0.50), y12_qcal2=weightedQuantile(y12, w_cal_q2,0.50),
                         y12_ipw=weightedQuantile(y12, w_ipw,0.50), y12_dr=weightedQuantile(y12, w_dr,0.50),
                         y12_mi = quantile(mi_y12_vals, 0.50))]
  
  res6 <- sample_bd2[, .(r=r, data="bd2", type = "q50",
                         y11_naive=quantile(y11,0.50), y11_cal=weightedQuantile(y11, w_cal_s,0.50),
                         y11_qcal1=weightedQuantile(y11, w_cal_q1,0.50), y11_qcal2=weightedQuantile(y11, w_cal_q2,0.50),
                         y11_ipw=weightedQuantile(y11, w_ipw,0.50), y11_dr=weightedQuantile(y11, w_dr,0.50),
                         y11_mi = quantile(mi_y11_vals_2, 0.50),
                         y12_naive=quantile(y12,0.50), y12_cal=weightedQuantile(y12, w_cal_s,0.50),
                         y12_qcal1=weightedQuantile(y12, w_cal_q1,0.50), y12_qcal2=weightedQuantile(y12, w_cal_q2,0.50),
                         y12_ipw=weightedQuantile(y12, w_ipw,0.50), y12_dr=weightedQuantile(y12, w_dr,0.50),
                         y12_mi = quantile(mi_y12_vals_2, 0.50))]
  
  
  res7 <- sample_bd1[, .(r=r, data="bd1", type = "q75",
                         y11_naive=quantile(y11,0.75), y11_cal=weightedQuantile(y11, w_cal_s,0.75),
                         y11_qcal1=weightedQuantile(y11, w_cal_q1,0.75), y11_qcal2=weightedQuantile(y11, w_cal_q2,0.75),
                         y11_ipw=weightedQuantile(y11, w_ipw,0.75), y11_dr=weightedQuantile(y11, w_dr,0.75),
                         y11_mi = quantile(mi_y11_vals, 0.75),
                         y12_naive=quantile(y12,0.75), y12_cal=weightedQuantile(y12, w_cal_s,0.75),
                         y12_qcal1=weightedQuantile(y12, w_cal_q1,0.75), y12_qcal2=weightedQuantile(y12, w_cal_q2,0.75),
                         y12_ipw=weightedQuantile(y12, w_ipw,0.75), y12_dr=weightedQuantile(y12, w_dr,0.75),
                         y12_mi = quantile(mi_y12_vals, 0.75))]
  
  res8 <- sample_bd2[, .(r=r, data="bd2", type = "q75",
                         y11_naive=quantile(y11,0.75), y11_cal=weightedQuantile(y11, w_cal_s,0.75),
                         y11_qcal1=weightedQuantile(y11, w_cal_q1,0.75), y11_qcal2=weightedQuantile(y11, w_cal_q2,0.75),
                         y11_ipw=weightedQuantile(y11, w_ipw,0.75), y11_dr=weightedQuantile(y11, w_dr,0.75),
                         y11_mi = quantile(mi_y11_vals_2, 0.75),
                         y12_naive=quantile(y12,0.75), y12_cal=weightedQuantile(y12, w_cal_s,0.75),
                         y12_qcal1=weightedQuantile(y12, w_cal_q1,0.75), y12_qcal2=weightedQuantile(y12, w_cal_q2,0.75),
                         y12_ipw=weightedQuantile(y12, w_ipw,0.75), y12_dr=weightedQuantile(y12, w_dr,0.75),
                         y12_mi = quantile(mi_y12_vals_2, 0.75))]
  
  
  
  results[[r]] <- rbind(res1,res2, res3, res4, res5, res6, res7, res8, fill = T)
}


# save results  -----------------------------------------------------------

results_df <- rbindlist(results) |> 
  melt(id.vars = c("r", "data", "type"), value.name = "value", variable.name = "estimator") |> 
  na.omit() |> 
  {\(x) x[, c("y", "estimator"):=tstrsplit(estimator, "_")]}() |> 
  {\(x) x[, estimator:=factor(estimator, c("naive", "cal", "ipw", "mi", "dr", "qcal1", "qcal2"))]}()

results_df[pop_true_vals, true := i.true, on = c("y", "type")]
results_df[, estimator:=factor(estimator, c("naive", "cal", "ipw", "mi", "dr", "qcal1", "qcal2"),
                               c("Naive", "Cal", "IPW", "MI", "DR", "QCal1", "QCal2"))]

saveRDS(object = results_df, file = "results/sim-yang2021-svymeth.rds")

fwrite(x = results_df, file = "results/sim-yang2021-svymeth.csv.gz")

