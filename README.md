
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Overview

## Details

A small package for joint calibration of totals and quantiles. The
package combines the following approaches:

- Deville, J. C., and Särndal, C. E. (1992). [Calibration estimators in
  survey
  sampling](https://www.tandfonline.com/doi/abs/10.1080/01621459.1992.10475217).
  Journal of the American statistical Association, 87(418), 376-382.
- Harms, T. and Duchesne, P. (2006). [On calibration estimation for
  quantiles](https://www150.statcan.gc.ca/n1/pub/12-001-x/2006001/article/9255-eng.pdf).
  Survey Methodology, 32(1), 37.

which allows to calibrate weights to known (or estimated) totals and
quantiles jointly. As an backend for calibration
[sampling](https://cran.r-project.org/web/packages/sampling)
(`sampling::calib`),
[laeken](https://cran.r-project.org/web/packages/laeken)
(`laeken::calibWeights`) or
[survey](https://cran.r-project.org/web/packages/survey/index.html)
(`survey::grake`) package can be used.

Currently supports:

- calibration of totals,
- calibration of quantiles.

Further plans:

- generalized calibration via `sampling::gencalib`,
- empirical likelihood,
- calibrated / covariate balancing propensity score,
- calibration for Gini and other metrices,
- observational studies / causal inference
- …

## Funding

Work on this package is supported by the the National Science Center,
OPUS 22 grant no. 2020/39/B/HS4/00941.

## Installation

You can install the development version of `jointCalib` from GitHub
with:

``` r
# install.packages("remotes")
remotes::install_github("ncn-foreigners/jointCalib")
```

## Examples

### Example 1

Based on Haziza, D., and Lesage, É. (2016). A discussion of weighting
procedures for unit nonresponse. Journal of Official Statistics, 32(1),
129-145.

``` r
library(jointCalib)
library(survey)
library(laeken)
```

``` r
N <- 1000
x <- runif(N, 0, 80)
y <- exp(-0.1 + 0.1*x) + rnorm(N, 0, 300)
p <- rbinom(N, 1, prob = exp(-0.2 - 0.014*x))
probs <- seq(0.1, 0.9, 0.1)
quants_known <- list(x=quantile(x, probs))
totals_known <- c(x=sum(x))
df <- data.frame(x, y, p)
df_resp <- df[df$p == 1, ]
df_resp$d <- N/nrow(df_resp)
y_quant_true <- quantile(y, probs)
head(df_resp)
#>           x         y p        d
#> 3  64.21569 1273.1933 1 1.976285
#> 4  10.56583  625.9272 1 1.976285
#> 5  67.25068  764.2867 1 1.976285
#> 8  23.18348  212.4432 1 1.976285
#> 10 53.45651  504.6195 1 1.976285
#> 11 57.75862  499.3678 1 1.976285
```

### Using `jointCalib` package

example 1a: calibrate only quantiles (deciles)

``` r
result1 <- joint_calib(formula_quantiles = ~x,
                      data=df_resp,
                      dweights=df_resp$d,
                      N = N,
                      pop_quantiles = quants_known,
                      method = "linear",
                      backend = "sampling")
summary(result1$g)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.6747  0.7230  0.8576  1.0000  1.2341  1.9462
```

``` r
data.frame(true = quants_known$x, est = weightedQuantile(df_resp$x, result1$g*df_resp$d, probs))
#>          true       est
#> 10%  7.946884  8.019299
#> 20% 16.229667 16.234197
#> 30% 23.120271 23.183478
#> 40% 31.319196 31.330988
#> 50% 39.977623 40.110044
#> 60% 47.636452 47.740290
#> 70% 56.053328 56.038232
#> 80% 64.266002 64.263914
#> 90% 72.268093 72.446851
```

example 1b: calibrate only quantiles (deciles)

``` r
result2 <- joint_calib(formula_totals = ~x,
                       formula_quantiles = ~x,
                       data = df_resp,
                       dweights = df_resp$d,
                       N = N,
                       pop_quantiles = quants_known,
                       pop_totals = totals_known,
                       method = "linear",
                       backend = "sampling")
summary(result2$g)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.6139  0.7203  0.8796  1.0000  1.2113  2.0022
```

``` r
data.frame(true = quants_known$x, est = weightedQuantile(df_resp$x, result2$g*df_resp$d, probs))
#>          true       est
#> 10%  7.946884  7.809372
#> 20% 16.229667 16.234197
#> 30% 23.120271 23.043652
#> 40% 31.319196 31.330988
#> 50% 39.977623 39.643152
#> 60% 47.636452 47.428222
#> 70% 56.053328 56.038232
#> 80% 64.266002 64.263914
#> 90% 72.268093 71.642378
```

We can restrict weights to specific range using `logit`.

``` r
result3 <- joint_calib(formula_totals = ~x,
                       formula_quantiles = ~x,
                       data = df_resp,
                       dweights = df_resp$d,
                       N = N,
                       pop_quantiles = quants_known,
                       pop_totals = totals_known,
                       method = "logit",
                       backend = "sampling", 
                       maxit = 500,
                       bounds = c(0, 2))

summary(result3$g)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.6142  0.7199  0.8818  1.0000  1.2099  1.9525
```

``` r
data.frame(true = quants_known$x, est = weightedQuantile(df_resp$x, result3$g*df_resp$d, probs))
#>          true       est
#> 10%  7.946884  7.809372
#> 20% 16.229667 16.234197
#> 30% 23.120271 23.043652
#> 40% 31.319196 31.301508
#> 50% 39.977623 39.643152
#> 60% 47.636452 47.428222
#> 70% 56.053328 56.038232
#> 80% 64.266002 64.263914
#> 90% 72.268093 71.642378
```

### Using `survey` package

example 1a: calibrate only quantiles (deciles)

``` r
A <- joint_calib_create_matrix(X_q = model.matrix(~x+0, df_resp), N = N, pop_quantiles = quants_known)
colnames(A) <- paste0("quant_", gsub("\\D", "", names(quants_known$x)))
A <- as.data.frame(A)
df_resp <- cbind(df_resp, A)
head(df_resp)
#>           x         y p        d quant_10 quant_20    quant_30 quant_40
#> 3  64.21569 1273.1933 1 1.976285        0    0.000 0.00000e+00    0.000
#> 4  10.56583  625.9272 1 1.976285        0    0.001 1.00000e-03    0.001
#> 5  67.25068  764.2867 1 1.976285        0    0.000 0.00000e+00    0.000
#> 8  23.18348  212.4432 1 1.976285        0    0.000 3.54316e-31    0.001
#> 10 53.45651  504.6195 1 1.976285        0    0.000 0.00000e+00    0.000
#> 11 57.75862  499.3678 1 1.976285        0    0.000 0.00000e+00    0.000
#>    quant_50 quant_60 quant_70 quant_80 quant_90
#> 3     0.000    0.000    0.000    0.001    0.001
#> 4     0.001    0.001    0.001    0.001    0.001
#> 5     0.000    0.000    0.000    0.000    0.001
#> 8     0.001    0.001    0.001    0.001    0.001
#> 10    0.000    0.000    0.001    0.001    0.001
#> 11    0.000    0.000    0.000    0.001    0.001
m1 <- svydesign(ids = ~1, data = df_resp, weights = ~d)
quants_formula <- as.formula(paste("~", paste(colnames(A), collapse = "+")))
quants_formula
#> ~quant_10 + quant_20 + quant_30 + quant_40 + quant_50 + quant_60 + 
#>     quant_70 + quant_80 + quant_90
svytotal(quants_formula, m1)
#>            total     SE
#> quant_10 0.14822 0.0158
#> quant_20 0.28461 0.0201
#> quant_30 0.42292 0.0220
#> quant_40 0.53953 0.0222
#> quant_50 0.63636 0.0214
#> quant_60 0.71739 0.0200
#> quant_70 0.79051 0.0181
#> quant_80 0.88121 0.0144
#> quant_90 0.94862 0.0098
```

Calibration using `calibrate`

``` r
pop_totals <- c(N, probs)
names(pop_totals) <- c("(Intercept)", colnames(A))
m1_cal <- calibrate(m1, quants_formula, pop_totals)
svytotal(quants_formula, m1_cal)
#>          total SE
#> quant_10   0.1  0
#> quant_20   0.2  0
#> quant_30   0.3  0
#> quant_40   0.4  0
#> quant_50   0.5  0
#> quant_60   0.6  0
#> quant_70   0.7  0
#> quant_80   0.8  0
#> quant_90   0.9  0
```

Calibration using `grake` (low level)

``` r
g1 <- grake(mm = as.matrix(cbind(1, A)),
            ww = df_resp$d,
            calfun = cal.linear,
            population = pop_totals,
            bounds = list(lower = -Inf, upper = Inf),
            epsilon = 1e-7,
            verbose = FALSE,
            maxit = 50,
            variance = NULL)
summary(as.numeric(g1))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.6747  0.7230  0.8576  1.0000  1.2341  1.9462
```
