
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
(`sampling::calib`) or
[laeken](https://cran.r-project.org/web/packages/laeken)
(`laeken::calibWeights`) package can be used.

Curently supports:

- calibration of totals
- calibration of quantiles

Further plans:

- calibration for Gini and other metrices
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
#>          x          y p        d
#> 1 17.06176 -694.23371 1 1.976285
#> 2 31.46439 -624.73315 1 1.976285
#> 4 17.83496   58.41928 1 1.976285
#> 5 76.45788 2328.42838 1 1.976285
#> 6 35.78455   78.80412 1 1.976285
#> 7 27.65472 -155.18430 1 1.976285
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
summary(result1$w)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   1.206   1.518   1.670   1.976   2.612   4.213
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
summary(result2$w)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   1.134   1.539   1.668   1.976   2.538   4.292
```

### Using `survey` package

example 1a: calibrate only quantiles (deciles)

``` r
A <- joint_calib_create_matrix(X_q = model.matrix(~x+0, df_resp), N = N, pop_quantiles = quants_known)
colnames(A) <- paste0("quant_", gsub("\\D", "", names(quants_known$x)))
A <- as.data.frame(A)
df_resp <- cbind(df_resp, A)
head(df_resp)
#>          x          y p        d quant_10 quant_20 quant_30 quant_40 quant_50
#> 1 17.06176 -694.23371 1 1.976285        0        0    0.001    0.001    0.001
#> 2 31.46439 -624.73315 1 1.976285        0        0    0.000    0.000    0.001
#> 4 17.83496   58.41928 1 1.976285        0        0    0.001    0.001    0.001
#> 5 76.45788 2328.42838 1 1.976285        0        0    0.000    0.000    0.000
#> 6 35.78455   78.80412 1 1.976285        0        0    0.000    0.000    0.001
#> 7 27.65472 -155.18430 1 1.976285        0        0    0.000    0.001    0.001
#>   quant_60 quant_70 quant_80 quant_90
#> 1    0.001    0.001    0.001    0.001
#> 2    0.001    0.001    0.001    0.001
#> 4    0.001    0.001    0.001    0.001
#> 5    0.000    0.000    0.000    0.000
#> 6    0.001    0.001    0.001    0.001
#> 7    0.001    0.001    0.001    0.001
m1 <- svydesign(ids = ~1, data = df_resp, weights = ~d)
quants_formula <- as.formula(paste("~", paste(colnames(A), collapse = "+")))
svytotal(quants_formula, m1)
#>            total     SE
#> quant_10 0.16383 0.0165
#> quant_20 0.29407 0.0203
#> quant_30 0.41245 0.0219
#> quant_40 0.53478 0.0222
#> quant_50 0.65125 0.0212
#> quant_60 0.73658 0.0196
#> quant_70 0.81229 0.0174
#> quant_80 0.88765 0.0140
#> quant_90 0.93468 0.0110
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
