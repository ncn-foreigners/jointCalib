library(survey)
library(laeken)
library(sampling)

source("functions.R")

## one simulation

set.seed(20230602)
N <- 100000
n_a <- 500
n_b <- 10000
n_b1 <- 0.7*n_b
n_b2 <- 0.3*n_b
x <- rchisq(N, 1)#rnorm(N, 2, 1)
x2 <- rbinom(N, 1, 0.7)
e <- rnorm(N)
y1 <- 1 + 2*x2 + 2*x + e
y2 <- 3 + x + 2*e
y3 <- 2.5 + 0.5*x^2 + e
strata <- x < median(x) ## MAR
#strata <- y1 < median(y1) ## NMAR
pop <- data.frame(x, x2, y1, y2, y3, strata)
sample_a <- pop[sample(1:N, n_a),]
sample_a$w_a <- N/n_a
pop1 <- subset(pop, strata == TRUE)
pop2 <- subset(pop, strata == FALSE)
sample_b <- rbind(pop1[sample(1:nrow(pop1), n_b1), ],
                  pop2[sample(1:nrow(pop2), n_b2), ])
sample_b$w_b <- N/n_b
p_quantiles <- seq(0.1, 0.9, 0.1)
tot_qualtiles <- quantile(x, p_quantiles)

## calibrate weights
w_res <- calib_quantiles(X_q = sample_b$x,  ## do odtworzenia rzędów kwantyli dla x
                         X = cbind(sample_b$x, sample_b$x2), ## to do totali
                         d = sample_b$w_b,  ## waga wejściowa
                         N=N,  ## wielkosc populacji
                         totals = c(sum(x), sum(x2)), ## wartosci globalne dla totali
                         totals_q = tot_qualtiles, ## kwantyle z populacji
                         method = "linear",
                         backend = "sampling")


sample_b$w_cal <- w_res$w
## before calibration
przed <- svydesign(ids= ~1, weights = ~w_b, data = sample_b)
## after calibration with quantiles
po <- svydesign(ids= ~1, weights = ~w_cal, data = sample_b)
## after calibration to totals
po_kal <- calibrate(design = przed, 
                    formula = ~x + x2, 
                    population = c(`(Intercept)`= N, x = sum(x), x2 = sum(x2)), calfun = cal.raking)

## results
results <- data.frame(
  before = svyquantile(~x, design = przed, quantiles = p_quantiles)$x[,1],
  calib = svyquantile(~x, design = po, quantiles = p_quantiles)$x[,1],
  cal_standard = svyquantile(~x, design = po_kal, quantiles = p_quantiles)$x[,1],
  known = tot_qualtiles)

results

## simulation study

B <- 500
results_sim <- matrix(0, nrow = B, ncol = 6)
for (b in 1:B) {
  set.seed(b)
  
  sample_b <- rbind(pop1[sample(1:nrow(pop1), n_b1), ],
                    pop2[sample(1:nrow(pop2), n_b2), ])
  sample_b$w_b <- N/n_b
  #sample_b <- sample_b[order(sample_b$x),]
  
  w_res <- calib_quantiles(X_q = sample_b$x, 
                           X = cbind(sample_b$x, sample_b$x2),
                           d = sample_b$w_b, 
                           N=N, 
                           totals = c(sum(x), sum(x2)),
                           totals_q = tot_qualtiles,
                           method = "linear",
                           backend = "sampling")
  
  sample_b$w_cal <- w_res$w
  przed <- svydesign(ids= ~1, weights = ~w_b, data = sample_b)
  po <- svydesign(ids= ~1, weights = ~w_cal, data = sample_b)
  po_kal <- calibrate(design = przed, formula = ~x, 
                      population = c(`(Intercept)`= N, x = sum(x)), 
                      calfun = cal.raking)
  
  results_sim[b,1] <- svymean(~y1, przed)[1]
  results_sim[b,2] <- svymean(~y1, po)[1]
  results_sim[b,3] <- svymean(~y1, po_kal)[1]
  
  results_sim[b,4] <- svyquantile(~y1, przed, quantiles = 0.5)$y1[1]
  results_sim[b,5] <- svyquantile(~y1, po, quantiles = 0.5)$y1[1]
  results_sim[b,6] <- svyquantile(~y1, po_kal, quantiles = 0.5)$y1[1]
}

## rmse dla mediany
sqrt((apply(results_sim[, 4:6], 2, mean)-median(y1))^2 + apply(results_sim[, 4:6], 2, var))
## rmse dla sredniej
sqrt((apply(results_sim[, 1:3], 2, mean)-mean(y1))^2 + apply(results_sim[, 1:3], 2, var))

colnames(results_sim) <- c("mean\nnaive", "mean\nkal-quant", "mean\nkal-tot", 
                           "median\nnaive", "median\nkal-quant", "median\nkal-tot")
boxplot(results_sim)
abline(h = mean(y1), col = "red")
abline(h = median(y1), col = "blue")

## NMAR via gencalib

set.seed(20230602)
N <- 100000
n_a <- 500
n_b <- 10000
n_b1 <- 0.7*n_b
n_b2 <- 0.3*n_b
x <- rchisq(N, 1)#rnorm(N, 2, 1)
x2 <- rbinom(N, 1, 0.7)
e <- rnorm(N)
y1 <- 1 + 2*x2 + 2*x + e
y2 <- 3 + x + 2*e
y3 <- 2.5 + 0.5*x^2 + e
#strata <- x < median(x) ## MAR
strata <- y1 < median(y1) ## NMAR
pop <- data.frame(x, x2, y1, y2, y3, strata)
sample_a <- pop[sample(1:N, n_a),]
sample_a$w_a <- N/n_a
pop1 <- subset(pop, strata == TRUE)
pop2 <- subset(pop, strata == FALSE)
sample_b <- rbind(pop1[sample(1:nrow(pop1), n_b1), ],
                  pop2[sample(1:nrow(pop2), n_b2), ])
sample_b$w_b <- N/n_b
p_quantiles <- seq(0.1, 0.9, 0.1)
tot_qualtiles <- quantile(x, p_quantiles)


B <- 500
ujemne_g3 <- ujemne_g2 <- 0

results_sim <- matrix(0, nrow = B, ncol = 6)
results_ci <- matrix(0, nrow=B,ncol=2)
for (b in 1:B) {
  set.seed(b)
  
  sample_b <- rbind(pop1[sample(1:nrow(pop1), n_b1), ],
                    pop2[sample(1:nrow(pop2), n_b2), ])
  sample_b$w_b <- N/n_b
  #sample_b <- sample_b[order(sample_b$x),]
  
  w_res <- calib_quantiles(X_q = sample_b$x, 
                           X = cbind(sample_b$x, sample_b$x2),
                           d = sample_b$w_b, 
                           N=N, 
                           totals = c(sum(x), sum(x2)),
                           totals_q = tot_qualtiles,
                           method = "linear",
                           backend = "sampling")
  
  sample_b$w_cal <- w_res$w
  przed <- svydesign(ids= ~1, weights = ~w_b, data = sample_b)
  po <- svydesign(ids= ~1, weights = ~w_cal, data = sample_b)
  
  Zs_x <- calib_quantiles_create_matrix(sample_b$x, 
                                        N, 
                                        quantile(sample_b$x, p_quantiles))
  
  Xs <- calib_quantiles_create_matrix(sample_b$x, N, tot_qualtiles)
  
  Xs_cal <- cbind(1, Xs, sample_b$x, sample_b$x2) ## all x
  sample_b$y1_med <- sample_b$y1 < median(sample_b$y1)
  sample_b$y1_perc <- sample_b$y1 < quantile(sample_b$y1, 0.75)
  
  Zs0 <- cbind(1, Xs, sample_b$y1, sample_b$x) ### instrumental variable y1
  Zs1 <- cbind(1, Xs, sample_b$y1_med, sample_b$y1) ### instrumental variable < med(y1) y1
  Zs2 <- cbind(1, Xs, sample_b$y1_perc, sample_b$y1) ### instrumental variable < perc(y1), y1
  Zs3 <- cbind(1, Xs, sample_b$y1_perc, sample_b$x) ### instrumental variable < perc(y1)
  
  totals <- c(N, p_quantiles, sum(x), sum(x2))
  d <- sample_b$w_b
  w_gcal0 <- gencalib(Xs = Xs_cal, Zs = Zs0, d = d, total = totals, method = "linear")   
  w_gcal1 <- gencalib(Xs = Xs_cal, Zs = Zs1, d = d, total = totals, method = "linear")
  w_gcal2 <- gencalib(Xs = Xs_cal, Zs = Zs2, d = d, total = totals, method = "linear")
  w_gcal3 <- gencalib(Xs = Xs_cal, Zs = Zs3, d = d, total = totals, method = "linear")
  
  if (any(w_gcal0 < 0)) print("ujemne wagi w w_gcal0")
  if (any(w_gcal1 < 0)) print("ujemne wagi w w_gcal1")
  if (any(w_gcal2 < 0)) {
    print("ujemne wagi w w_gcal2")
    ujemne_g2 <- ujemne_g2+1
  }
  if (any(w_gcal3 < 0)) {
    print("ujemne wagi w w_gcal3")
    ujemne_g3 <- ujemne_g3+1
  }
  
  sample_b$w_gcal0 <- w_gcal0*d
  sample_b$w_gcal1 <- w_gcal1*d
  sample_b$w_gcal2 <- w_gcal2*d
  sample_b$w_gcal3 <- w_gcal3*d
  
  po_gcal0 <- svydesign(ids= ~1, weights = ~w_gcal0, data = sample_b)
  po_gcal1 <- svydesign(ids= ~1, weights = ~w_gcal1, data = sample_b)
  po_gcal2 <- svydesign(ids= ~1, weights = ~w_gcal2, data = sample_b)
  po_gcal3 <- svydesign(ids= ~1, weights = ~w_gcal3, data = sample_b)
  
  results_sim[b,1] <- svyquantile(~y1, przed, quantiles = 0.5)$y1[1]
  results_sim[b,2] <- svyquantile(~y1, po, quantiles = 0.5)$y1[1]
  results_sim[b,3] <- svyquantile(~y1, po_gcal0, quantiles = 0.5)$y1[1]
  results_sim[b,4] <- svyquantile(~y1, po_gcal1, quantiles = 0.5)$y1[1]
  results_sim[b,5] <- svyquantile(~y1, po_gcal2, quantiles = 0.5)$y1[1]
  results_sim[b,6] <- svyquantile(~y1, po_gcal3, quantiles = 0.5)$y1[1]
  
  results_ci[b,] <- svyquantile(~y1, po_gcal3, quantiles = 0.5)$y1[c(2,3)]
}


boxplot(results_sim, ylim = c(3,4))
abline(h =  median(y1), col = "blue")

## rmse dla mediany
sqrt((apply(results_sim, 2, mean)-median(y1))^2 + apply(results_sim, 2, var))


