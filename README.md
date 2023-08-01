
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
#> 2  74.92869 1821.4218 1 2.070393
#> 5  10.27340  193.4621 1 2.070393
#> 7  69.67968  832.8965 1 2.070393
#> 9  72.48948 1826.4118 1 2.070393
#> 10 13.31971  381.5234 1 2.070393
#> 11 27.76119  306.9979 1 2.070393
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
#>   1.216   1.611   1.931   2.070   2.629   3.499
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
#>   1.081   1.637   1.960   2.070   2.575   3.602
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
                       bounds = c(0.5, 2))

summary(result3$w/df_resp$d)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.5695  0.7870  0.9437  1.0000  1.2332  1.7322
```

### Using `survey` package

example 1a: calibrate only quantiles (deciles)

``` r
A <- joint_calib_create_matrix(X_q = model.matrix(~x+0, df_resp), N = N, pop_quantiles = quants_known)
colnames(A) <- paste0("quant_", gsub("\\D", "", names(quants_known$x)))
A <- as.data.frame(A)
df_resp <- cbind(df_resp, A)
head(df_resp)
#>           x         y p        d quant_10 quant_20 quant_30 quant_40 quant_50
#> 2  74.92869 1821.4218 1 2.070393        0    0.000    0.000    0.000    0.000
#> 5  10.27340  193.4621 1 2.070393        0    0.001    0.001    0.001    0.001
#> 7  69.67968  832.8965 1 2.070393        0    0.000    0.000    0.000    0.000
#> 9  72.48948 1826.4118 1 2.070393        0    0.000    0.000    0.000    0.000
#> 10 13.31971  381.5234 1 2.070393        0    0.001    0.001    0.001    0.001
#> 11 27.76119  306.9979 1 2.070393        0    0.000    0.000    0.001    0.001
#>    quant_60 quant_70 quant_80 quant_90
#> 2     0.000    0.000    0.000    0.000
#> 5     0.001    0.001    0.001    0.001
#> 7     0.000    0.000    0.000    0.001
#> 9     0.000    0.000    0.000    0.000
#> 10    0.001    0.001    0.001    0.001
#> 11    0.001    0.001    0.001    0.001
m1 <- svydesign(ids = ~1, data = df_resp, weights = ~d)
quants_formula <- as.formula(paste("~", paste(colnames(A), collapse = "+")))
quants_formula
#> ~quant_10 + quant_20 + quant_30 + quant_40 + quant_50 + quant_60 + 
#>     quant_70 + quant_80 + quant_90
svytotal(quants_formula, m1)
#>            total     SE
#> quant_10 0.17013 0.0171
#> quant_20 0.29870 0.0208
#> quant_30 0.41346 0.0224
#> quant_40 0.51577 0.0228
#> quant_50 0.62290 0.0221
#> quant_60 0.70969 0.0207
#> quant_70 0.78847 0.0186
#> quant_80 0.86377 0.0156
#> quant_90 0.92307 0.0121
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
