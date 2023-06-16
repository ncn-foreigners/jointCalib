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
strata <- x <= median(x)
pop <- data.frame(x, x2, y1, y2, y3, strata)
sample_a <- pop[sample(1:N, n_a),]
sample_a$w_a <- N/n_a
svy_a <- svydesign(ids= ~1, weights = ~ w_a, data = sample_a)
svy_a_cal <- calibrate(svy_a, formula= ~x, population=c(`(Intercept)`=N, x = sum(x)))
pop1 <- subset(pop, strata == TRUE)
pop2 <- subset(pop, strata == FALSE)

sample_b <- rbind(pop1[sample(1:nrow(pop1), n_b1), ],
                  pop2[sample(1:nrow(pop2), n_b2), ])
sample_b$w_b <- N/n_b

p_quantiles <- seq(0.1, 0.9, 0.1)
tot_qualtiles <- quantile(x, p_quantiles)

## calibrate weights
w_res <- calib_quantiles(y = sample_b$y1, 
                         x = sample_b$x,  ## calibration for quantiles of x
                         X_cat = cbind(sample_b$x, sample_b$x2), ## for totals
                         d = sample_b$w_b,  ## design weight
                         N=N,  ## population size
                         totals = c(sum(x), sum(x2)), ## totals for X_cat
                         quantiles = p_quantiles, ## probabilities (alpha / tau)
                         pop_quantiles = tot_qualtiles, ## corresponding population quantile
                         method = "linear", 
                         backend = "sampling")


sample_b$w_cal <- w_res$w
## before calibration
przed <- svydesign(ids= ~1, weights = ~w_b, data = sample_b)
## after calibration with quantiles
po <- svydesign(ids= ~1, weights = ~w_cal, data = sample_b)
## after calibration to totals
po_kal <- calibrate(design = przed, formula = ~x, population = c(`(Intercept)`= N, x = sum(x)), calfun = cal.raking)

## results
results <- data.frame(
before = svyquantile(~x, design = przed, quantiles = p_quantiles)$x[,1],
calib = svyquantile(~x, design = po, quantiles = p_quantiles)$x[,1],
known = tot_qualtiles)

## simulation study

B <- 500
results_sim <- matrix(0, nrow = B, ncol = 6)
for (b in 1:B) {
  set.seed(b)
  
  sample_b <- rbind(pop1[sample(1:nrow(pop1), n_b1), ],
                    pop2[sample(1:nrow(pop2), n_b2), ])
  sample_b$w_b <- N/n_b
  #sample_b <- sample_b[order(sample_b$x),]
  
  w_res <- calib_quantiles(y = sample_b$y1, 
                           x = sample_b$x, 
                           X_cat = cbind(sample_b$x, sample_b$x2),
                           d = sample_b$w_b, 
                           N=N, 
                           totals = c(sum(x), sum(x2)),
                           quantiles = p_quantiles,
                           pop_quantiles = tot_qualtiles,
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

## rmse dla mean
sqrt((apply(results_sim[, 1:3], 2, mean)-mean(y1))^2 + apply(results_sim[, 1:3], 2, var))

## rmse dla median
sqrt((apply(results_sim[, 4:6], 2, mean)-median(y1))^2 + apply(results_sim[, 4:6], 2, var))

colnames(results_sim) <- c("mean\nnaive", "mean\nkal-quant", "mean\nkal-tot", 
                       "median\nnaive", "median\nkal-quant", "median\nkal-tot")
boxplot(results_sim)
abline(h = mean(y1), col = "red")
abline(h = median(y1), col = "blue")



