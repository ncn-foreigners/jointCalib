
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![R-CMD-check](https://github.com/ncn-foreigners/jointCalib/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ncn-foreigners/jointCalib/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

# Overview

## Details

A small package for joint calibration of totals and quantiles (see
[Beręsewicz and Szymkowiak (2023)](https://arxiv.org/abs/2308.13281)
working paper for details). The package combines the following
approaches:

- Deville, J. C., and Särndal, C. E. (1992). [Calibration estimators in
  survey
  sampling](https://www.tandfonline.com/doi/abs/10.1080/01621459.1992.10475217).
  Journal of the American statistical Association, 87(418), 376-382.
- Harms, T. and Duchesne, P. (2006). [On calibration estimation for
  quantiles](https://www150.statcan.gc.ca/n1/pub/12-001-x/2006001/article/9255-eng.pdf).
  Survey Methodology, 32(1), 37.
- Wu, C. (2005) [Algorithms and R codes for the pseudo empirical
  likelihood method in survey
  sampling](https://www150.statcan.gc.ca/n1/pub/12-001-x/2005002/article/9051-eng.pdf),
  Survey Methodology, 31(2), 239.
- Zhang, S., Han, P., and Wu, C. (2023) [Calibration Techniques
  Encompassing Survey Sampling, Missing Data Analysis and Causal
  Inference](https://onlinelibrary.wiley.com/doi/10.1111/insr.12518),
  International Statistical Review 91, 165–192.

which allows to calibrate weights to known (or estimated) totals and
quantiles jointly. As an backend for calibration
[sampling](https://CRAN.R-project.org/package=sampling)
(`sampling::calib`), [laeken](https://CRAN.R-project.org/package=laeken)
(`laeken::calibWeights`),
[survey](https://CRAN.R-project.org/package=survey) (`survey::grake`) or
[ebal](https://CRAN.R-project.org/package=ebal) (`ebal::eb`) package can
be used. One can also apply empirical likelihood using codes from Wu
(2005) with support of `stats::constrOptim` as used in Zhang, Han and Wu
(2022).

| backend    | method                                        | function called                 |
|------------|-----------------------------------------------|---------------------------------|
| `sampling` | `c("raking", "linear", "logit", "truncated")` | `sampling::calib`               |
| `laeken`   | `c("raking", "linear", "logit")`              | `laeken::calibWeights`          |
| `survey`   | `c("raking", "linear", "logit", "sinh")`      | `survey::grake`                 |
| `ebal`     | `eb`                                          | `ebal::eb`                      |
| `base`     | `el`                                          | R code and `stats::constrOptim` |

Currently supports:

- calibration of quantiles,
- calibration of quantiles and totals,
- calibration using standard calibration, empirical likelihood and
  entropy balancing method.

Further plans:

- generalized calibration via `sampling::gencalib`,
- calibrated / covariate balancing propensity score,
- calibration for Gini and other metrics,
- observational studies / causal inference
- …

## Funding

Work on this package is supported by the the National Science Centre,
OPUS 22 grant no. 2020/39/B/HS4/00941.

## Installation

You can install the development version of `jointCalib` from GitHub
with:

``` r
# install.packages("remotes")
remotes::install_github("ncn-foreigners/jointCalib")
```

## Examples

Load packages

``` r
library(jointCalib)
library(survey)
library(laeken)
library(ebal)
```

### Example 1 – census case

Based on Haziza, D., and Lesage, É. (2016). A discussion of weighting
procedures for unit nonresponse. Journal of Official Statistics, 32(1),
129-145.

``` r
set.seed(20230817)
N <- 1000
x <- runif(N, 0, 80)
#y <- 1000+10*x+rnorm(N, 0, 300)
#y <- exp(-0.1 + 0.1*x) + rnorm(N, 0, 300)
#y <- rbinom(N, 1, prob = 1/(exp(-0.5*(x-55))+1))
y <- 1300-(x-40)^2 + rnorm(N, 0, 300)
#p <- rbinom(N, 1, prob = 0.2+0.6*(1 + exp(-5 + x/8))^-1)
p <- rbinom(N, 1, prob = 0.07+0.45*(x/40-1)^2+0.0025*x)
#p <- rbinom(N, 1, prob = (1.2+0.024*x)^-1)
#p <- rbinom(N, 1, prob = exp(-0.2 - 0.014*x))
probs <- seq(0.1, 0.9, 0.1)
quants_known <- list(x=quantile(x, probs))
totals_known <- c(x=sum(x))
df <- data.frame(x, y, p)
df_resp <- df[df$p == 1, ]
df_resp$d <- N/nrow(df_resp)
y_quant_true <- quantile(y, probs)
head(df_resp)
#>           x        y p        d
#> 6  12.35687 444.9053 1 3.134796
#> 7  61.90319 403.9473 1 3.134796
#> 13 60.96079 923.4975 1 3.134796
#> 14 76.85300 124.4110 1 3.134796
#> 18 71.52828 422.0934 1 3.134796
#> 19 65.32808 740.4801 1 3.134796
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
result1
#> Weights calibrated using: linear withsamplingbackend.
#> Summary statistics for g-weights:
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.4562  0.6773  0.8180  1.0000  1.4500  2.4538 
#> Totals and precision (abs diff: 4.574665e-10)
#>       totals    precision
#>  [1,]  1e+03 4.557705e-10
#>  [2,]  1e-01 5.020984e-14
#>  [3,]  2e-01 1.055545e-13
#>  [4,]  3e-01 1.324496e-13
#>  [5,]  4e-01 1.580958e-13
#>  [6,]  5e-01 1.765255e-13
#>  [7,]  6e-01 1.991740e-13
#>  [8,]  7e-01 2.304823e-13
#>  [9,]  8e-01 2.882139e-13
#> [10,]  9e-01 3.552714e-13
```

``` r
data.frame(true = quants_known$x, est = weightedQuantile(df_resp$x, result1$g*df_resp$d, probs))
#>          true       est
#> 10%  7.078067  7.085003
#> 20% 14.831221 14.824424
#> 30% 23.146180 23.287657
#> 40% 31.641911 31.802986
#> 50% 39.033812 39.154276
#> 60% 47.527168 48.252065
#> 70% 54.984229 55.311953
#> 80% 64.073167 64.062629
#> 90% 71.565441 71.567274
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
#>  0.3563  0.6208  0.8199  1.0000  1.4368  2.5384
```

``` r
data.frame(true = quants_known$x, est = weightedQuantile(df_resp$x, result2$g*df_resp$d, probs))
#>          true       est
#> 10%  7.078067  7.085003
#> 20% 14.831221 14.824424
#> 30% 23.146180 23.287657
#> 40% 31.641911 31.802986
#> 50% 39.033812 39.154276
#> 60% 47.527168 48.252065
#> 70% 54.984229 55.311953
#> 80% 64.073167 64.062629
#> 90% 71.565441 71.567274
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
                       bounds = c(0, 3))

summary(result3$g)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.3911  0.6254  0.8140  1.0000  1.4325  2.5186
```

``` r
data.frame(true = quants_known$x, est = weightedQuantile(df_resp$x, result3$g*df_resp$d, probs))
#>          true       est
#> 10%  7.078067  7.085003
#> 20% 14.831221 14.824424
#> 30% 23.146180 22.872078
#> 40% 31.641911 31.297922
#> 50% 39.033812 39.154276
#> 60% 47.527168 48.252065
#> 70% 54.984229 55.311953
#> 80% 64.073167 64.529683
#> 90% 71.565441 71.567274
```

Empirical likelihood method can be applied using the following code

``` r
result4a <- joint_calib(formula_quantiles = ~ x,
                        data = df_resp,
                        dweights = df_resp$d,
                        N = N,
                        pop_quantiles = quants_known,
                        method = "el",
                        backend = "base")
summary(result4a$g)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.4563  0.6773  0.8180  1.0000  1.4500  2.4538

result4b <- joint_calib(formula_totals = ~ x,
                        formula_quantiles = ~ x,
                        data = df_resp,
                        dweights = df_resp$d,
                        N = N,
                        pop_quantiles = quants_known,
                        pop_totals = totals_known,
                        method = "el",
                        backend = "base")

summary(result4b$g)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.4414  0.6584  0.8109  1.0000  1.4256  2.8513

result4c <- calib_el(X = result4b$Xs[, c(1, 11)],
                     d = df_resp$d, 
                     totals = c(N,totals_known),
                     maxit = 50,
                     tol = 1e-8)

summary(result4c)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.7438  0.7910  0.8848  1.0000  1.2455  1.4923
```

``` r
data.frame(true = quants_known$x, 
           est = weightedQuantile(df_resp$x, result2$g*df_resp$d, probs),
           est_el0 = weightedQuantile(df_resp$x, result4c*df_resp$d, probs),
           est_el1 = weightedQuantile(df_resp$x, result4a$g*df_resp$d, probs),
           est_el2 = weightedQuantile(df_resp$x, result4b$g*df_resp$d, probs))
#>          true       est   est_el0   est_el1   est_el2
#> 10%  7.078067  7.085003  3.640246  7.085003  7.085003
#> 20% 14.831221 14.824424  8.024948 14.824424 14.824424
#> 30% 23.146180 23.287657 14.499018 23.287657 23.287657
#> 40% 31.641911 31.802986 24.010501 31.802986 31.802986
#> 50% 39.033812 39.154276 39.154276 38.892342 39.154276
#> 60% 47.527168 48.252065 54.685438 48.252065 48.252065
#> 70% 54.984229 55.311953 63.084405 54.954054 54.954054
#> 80% 64.073167 64.062629 70.050593 64.062629 64.062629
#> 90% 71.565441 71.567274 75.663078 71.567274 71.567274
```

Entropy balancing method can be applied using the following code

``` r
result5 <- joint_calib(formula_quantiles = ~ x,
                        data = df_resp,
                        dweights = df_resp$d,
                        N = N,
                        pop_quantiles = quants_known,
                        method = "eb")
summary(result5$g)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.4576  0.6775  0.8180  1.0003  1.4500  2.4538
```

Finally, compare all method with true Y distribution where

- `est_y1` – weights calibrated to quantiles only,
- `est_y2` – weights calibrated to quantiles and totals,
- `est_y3` – weights calibrated to quantiles and totals using `logit`
  distance and range limitations,
- `est_y4` – weights calibrated to means only using `el`.
- `est_y5` – weights calibrated to quantiles only using `el`.
- `est_y6` – weights calibrated to quantiles and totals using `el`.
- `est_y7` – weights calibrated to quantiles and totals using `eb`.

``` r
results_y <- data.frame(true_y = y_quant_true,
                        est_y1 = weightedQuantile(df_resp$y, result1$g * df_resp$d, probs),
                        est_y2 = weightedQuantile(df_resp$y, result2$g * df_resp$d, probs),
                        est_y3 = weightedQuantile(df_resp$y, result3$g * df_resp$d, probs),
                        est_y4 = weightedQuantile(df_resp$y, result4c * df_resp$d, probs),
                        est_y5 = weightedQuantile(df_resp$y, result4a$g * df_resp$d, probs),
                        est_y6 = weightedQuantile(df_resp$y, result4b$g * df_resp$d, probs),
                        est_y7 = weightedQuantile(df_resp$y, result5$g * df_resp$d, probs))

results_y
#>         true_y     est_y1     est_y2     est_y3     est_y4     est_y5
#> 10%  -56.97982  -36.18473  -36.18473  -36.18473 -149.00996  -36.18473
#> 20%  210.58301  228.75429  228.75429  228.75429   17.27107  228.75429
#> 30%  444.81652  413.22517  417.60155  413.22517  203.58302  413.22517
#> 40%  646.06099  622.82192  622.82192  622.82192  381.99560  622.82192
#> 50%  803.50842  789.89936  789.89936  789.89936  503.35849  789.89936
#> 60%  975.31803  925.84730  925.84730  925.84730  679.04530  925.84730
#> 70% 1123.44643 1023.75230 1023.75230 1023.75230  811.84389 1023.75230
#> 80% 1245.62645 1162.06504 1162.06504 1162.06504 1009.46703 1162.06504
#> 90% 1449.23819 1471.33301 1471.33301 1471.33301 1242.60357 1471.33301
#>         est_y6     est_y7
#> 10%  -36.18473  -36.18473
#> 20%  228.75429  228.75429
#> 30%  411.47088  413.22517
#> 40%  622.82192  622.82192
#> 50%  789.89936  789.89936
#> 60%  925.84730  925.84730
#> 70% 1023.75230 1023.75230
#> 80% 1162.06504 1162.06504
#> 90% 1471.33301 1471.33301
```

Here is the sum of squares and calibration of totals and means seems to
be the best.

``` r
apply(results_y[, -1], 2, FUN = function(x) sum((x-results_y[,1])^2))
#>    est_y1    est_y2    est_y3    est_y4    est_y5    est_y6    est_y7 
#>  22342.87  22085.51  22342.87 547195.99  22342.87  22456.79  22342.87
```

### Using `survey` package

example 1a: calibrate only quantiles (deciles)

``` r
A <- joint_calib_create_matrix(X_q = model.matrix(~x+0, df_resp), N = N, pop_quantiles = quants_known)
colnames(A) <- paste0("quant_", gsub("\\D", "", names(quants_known$x)))
A <- as.data.frame(A)
df_resp <- cbind(df_resp, A)
head(df_resp)
#>           x        y p        d quant_10 quant_20 quant_30 quant_40 quant_50
#> 6  12.35687 444.9053 1 3.134796        0    0.001    0.001    0.001    0.001
#> 7  61.90319 403.9473 1 3.134796        0    0.000    0.000    0.000    0.000
#> 13 60.96079 923.4975 1 3.134796        0    0.000    0.000    0.000    0.000
#> 14 76.85300 124.4110 1 3.134796        0    0.000    0.000    0.000    0.000
#> 18 71.52828 422.0934 1 3.134796        0    0.000    0.000    0.000    0.000
#> 19 65.32808 740.4801 1 3.134796        0    0.000    0.000    0.000    0.000
#>    quant_60 quant_70 quant_80 quant_90
#> 6     0.001    0.001    0.001    0.001
#> 7     0.000    0.000    0.001    0.001
#> 13    0.000    0.000    0.001    0.001
#> 14    0.000    0.000    0.000    0.000
#> 18    0.000    0.000    0.000    0.001
#> 19    0.000    0.000    0.000    0.001
m1 <- svydesign(ids = ~1, data = df_resp, weights = ~d)
quants_formula <- as.formula(paste("~", paste(colnames(A), collapse = "+")))
quants_formula
#> ~quant_10 + quant_20 + quant_30 + quant_40 + quant_50 + quant_60 + 
#>     quant_70 + quant_80 + quant_90
svytotal(quants_formula, m1)
#>            total     SE
#> quant_10 0.10972 0.0175
#> quant_20 0.23197 0.0237
#> quant_30 0.29154 0.0255
#> quant_40 0.34796 0.0267
#> quant_50 0.38871 0.0273
#> quant_60 0.43887 0.0278
#> quant_70 0.50784 0.0280
#> quant_80 0.63323 0.0270
#> quant_90 0.78100 0.0232
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
#>  0.4562  0.6773  0.8180  1.0000  1.4500  2.4538
```

## Example 2 – non-probability samples

## Example 3 – observational studies / causal inference
