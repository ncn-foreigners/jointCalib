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


generate_calibrate <- function(data, totals, N, method, quants) {
  ## p1
  sample_nonresp1 <- data[rbinom(N, 1, p1)==1, ]
  sample_nonresp1[, d := N/.N]
  
  w0 <- calib(Xs = cbind(1, as.matrix(sample_nonresp1$x)),
              d = sample_nonresp1$d,
              total = c(N, totals),
              method = method)*sample_nonresp1$d
  
  w1 <- calib_quantiles(X_q = as.matrix(sample_nonresp1$x),
                        d = sample_nonresp1$d,
                        N = N,
                        totals_q = list(quantile(df$x, quants)),
                        method = method,
                        backend = "sampling")
  
  w2 <- calib_quantiles(X_q = as.matrix(sample_nonresp1$x),
                        X = as.matrix(sample_nonresp1$x),
                        d = sample_nonresp1$d,
                        N = N,
                        totals = totals,
                        totals_q = list(quantile(df$x, quants)),
                        method = method,
                        backend = "sampling")
  
  res1 <- rbind(
    sample_nonresp1[, lapply(.SD, weighted.mean, w = w0), .SDcols = patterns("y")][, ":="(w="w0", p="p1", m = method)],
    sample_nonresp1[, lapply(.SD, weighted.mean, w = w1$w), .SDcols = patterns("y")][, ":="(w="w1", p="p1", m = method)],
    sample_nonresp1[, lapply(.SD, weighted.mean, w = w2$w), .SDcols = patterns("y")][, ":="(w="w2", p="p1", m = method)])
  
  ## p2
  sample_nonresp1 <- data[rbinom(N, 1, p2)==1, ]
  sample_nonresp1[, d := N/.N]
  
  w0 <- calib(Xs = cbind(1, as.matrix(sample_nonresp1$x)),
              d = sample_nonresp1$d,
              total = c(N, totals),
              method = method)*sample_nonresp1$d
  
  w1 <- calib_quantiles(X_q = as.matrix(sample_nonresp1$x),
                        d = sample_nonresp1$d,
                        N = N,
                        totals_q = list(quantile(df$x, quants)),
                        method = method,
                        backend = "sampling")
  
  w2 <- calib_quantiles(X_q = as.matrix(sample_nonresp1$x),
                        X = as.matrix(sample_nonresp1$x),
                        d = sample_nonresp1$d,
                        N = N,
                        totals = totals,
                        totals_q = list(quantile(df$x, quants)),
                        method = method,
                        backend = "sampling")
  
  res2 <- rbind(
    sample_nonresp1[, lapply(.SD, weighted.mean, w = w0), .SDcols = patterns("y")][, ":="(w="w0", p="p2", m = method)],
    sample_nonresp1[, lapply(.SD, weighted.mean, w = w1$w), .SDcols = patterns("y")][, ":="(w="w1", p="p2", m = method)],
    sample_nonresp1[, lapply(.SD, weighted.mean, w = w2$w), .SDcols = patterns("y")][, ":="(w="w2", p="p2", m = method)])
  
  ## p3
  sample_nonresp1 <- data[rbinom(N, 1, p3)==1, ]
  sample_nonresp1[, d := N/.N]
  
  w0 <- calib(Xs = cbind(1, as.matrix(sample_nonresp1$x)),
              d = sample_nonresp1$d,
              total = c(N, totals),
              method = method)*sample_nonresp1$d
  
  w1 <- calib_quantiles(X_q = as.matrix(sample_nonresp1$x),
                        d = sample_nonresp1$d,
                        N = N,
                        totals_q = list(quantile(df$x, quants)),
                        method = method,
                        backend = "sampling")
  
  w2 <- calib_quantiles(X_q = as.matrix(sample_nonresp1$x),
                        X = as.matrix(sample_nonresp1$x),
                        d = sample_nonresp1$d,
                        N = N,
                        totals = totals,
                        totals_q = list(quantile(df$x, quants)),
                        method = method,
                        backend = "sampling")
  
  res3 <- rbind(
    sample_nonresp1[, lapply(.SD, weighted.mean, w = w0), .SDcols = patterns("y")][, ":="(w="w0", p="p3", m = method)],
    sample_nonresp1[, lapply(.SD, weighted.mean, w = w1$w), .SDcols = patterns("y")][, ":="(w="w1", p="p3", m = method)],
    sample_nonresp1[, lapply(.SD, weighted.mean, w = w2$w), .SDcols = patterns("y")][, ":="(w="w2", p="p3", m = method)])
  
  ## p4
  sample_nonresp1 <- data[rbinom(N, 1, p4)==1, ]
  sample_nonresp1[, d := N/.N]
  
  w0 <- calib(Xs = cbind(1, as.matrix(sample_nonresp1$x)),
              d = sample_nonresp1$d,
              total = c(N, totals),
              method = method)*sample_nonresp1$d
  
  w1 <- calib_quantiles(X_q = as.matrix(sample_nonresp1$x),
                        d = sample_nonresp1$d,
                        N = N,
                        totals_q = list(quantile(df$x, quants)),
                        method = method,
                        backend = "sampling")
  
  w2 <- calib_quantiles(X_q = as.matrix(sample_nonresp1$x),
                        X = as.matrix(sample_nonresp1$x),
                        d = sample_nonresp1$d,
                        N = N,
                        totals = totals,
                        totals_q = list(quantile(df$x, quants)),
                        method = method,
                        backend = "sampling")
  
  res4 <- rbind(
    sample_nonresp1[, lapply(.SD, weighted.mean, w = w0), .SDcols = patterns("y")][, ":="(w="w0", p="p4", m = method)],
    sample_nonresp1[, lapply(.SD, weighted.mean, w = w1$w), .SDcols = patterns("y")][, ":="(w="w1", p="p4", m = method)],
    sample_nonresp1[, lapply(.SD, weighted.mean, w = w2$w), .SDcols = patterns("y")][, ":="(w="w2", p="p4", m = method)])
  
  return(rbind(res1, res2, res3, res4))
  
}

seed_for_sim <- 2023-07-14
set.seed(seed_for_sim+1)
N <- 1000
x <- runif(N, 0, 80)
y1 <- 1000+10*x+rnorm(N, 0, 300)
y2 <- exp(-0.1 + 0.1*x) + rnorm(N, 0, 300)
y3 <- rbinom(N, 1, prob = 1/(exp(-0.5*(x-55))+1))
y4 <- 1300-(x-40)^2 + rnorm(N, 0, 300)
p1 <- (1.2+0.024*x)^-1
p2 <- exp(-0.2 - 0.014*x)
p3 <- 0.2+0.6*(1 + exp(-5 + x/8))^-1
p4 <- 0.07+0.45*(x/40-1)^2+0.0025*x
df <- data.frame(x,y1,y2,y3,y4,p1,p2,p3,p4) |> setDT()
quants <- seq(0.1, 0.9, 0.1)

results <- list()
R <- 1000
for (r in 1:R) {
  set.seed(r)
  if (r %% 100 == 0) cat("iteration:", r, "\n")
  results[[r]] <- rbind(
    generate_calibrate(data=df,totals=sum(x), N=N, method = "linear", quants = quants),
    generate_calibrate(data=df,totals=sum(x), N=N, method = "raking", quants = quants),
    generate_calibrate(data=df,totals=sum(x), N=N, method = "logit", quants = quants)
    )
}

pop_true <- df[, lapply(.SD, mean), .SDcols = patterns("y")] |> melt(value.var = "value", variable.name = "y")
pop_quants <- df[, lapply(.SD, quantile, p=0.5), 
                 .SDcols = patterns("y")] |> melt(value.var = "value", variable.name = "y")

results_all <- rbindlist(results, idcol = "r") |> 
  melt(id.vars = c("r", "w", "p", "m"), value.var = "value", variable.name = "y") 

results_all[pop_true, true:=i.value, on = "y"]
results_all[pop_true, q50:=i.value, on = "y"]

saveRDS(results_all, file = "results/sim-haziza2016-jos.rds")
fwrite(x = results_all, file = "results/sim-haziza2016-jos.csv.gz")

