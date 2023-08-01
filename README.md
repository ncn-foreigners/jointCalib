
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
#>           x         y p       d
#> 1  30.43356  290.8317 1 2.04918
#> 3  21.19021 -181.8232 1 2.04918
#> 9  67.57426 1020.7038 1 2.04918
#> 11 55.74922  155.6848 1 2.04918
#> 12 76.09029 2305.1251 1 2.04918
#> 13 22.25656  272.1227 1 2.04918
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
#>   1.300   1.494   1.790   2.049   2.391   3.173
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
#>  0.8057  1.5897  1.9497  2.0492  2.5069  3.6588
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
#>  0.5486  0.7419  0.9370  1.0000  1.2265  1.7946
```

### Using `survey` package

example 1a: calibrate only quantiles (deciles)

``` r
A <- joint_calib_create_matrix(X_q = model.matrix(~x+0, df_resp), N = N, pop_quantiles = quants_known)
colnames(A) <- paste0("quant_", gsub("\\D", "", names(quants_known$x)))
A <- as.data.frame(A)
df_resp <- cbind(df_resp, A)
head(df_resp)
#>           x         y p       d quant_10 quant_20 quant_30 quant_40 quant_50
#> 1  30.43356  290.8317 1 2.04918        0        0    0.000    0.000    0.001
#> 3  21.19021 -181.8232 1 2.04918        0        0    0.001    0.001    0.001
#> 9  67.57426 1020.7038 1 2.04918        0        0    0.000    0.000    0.000
#> 11 55.74922  155.6848 1 2.04918        0        0    0.000    0.000    0.000
#> 12 76.09029 2305.1251 1 2.04918        0        0    0.000    0.000    0.000
#> 13 22.25656  272.1227 1 2.04918        0        0    0.001    0.001    0.001
#>    quant_60 quant_70 quant_80 quant_90
#> 1     0.001    0.001    0.001    0.001
#> 3     0.001    0.001    0.001    0.001
#> 9     0.000    0.000    0.000    0.001
#> 11    0.000    0.000    0.001    0.001
#> 12    0.000    0.000    0.000    0.000
#> 13    0.001    0.001    0.001    0.001
m1 <- svydesign(ids = ~1, data = df_resp, weights = ~d)
quants_formula <- as.formula(paste("~", paste(colnames(A), collapse = "+")))
quants_formula
#> ~quant_10 + quant_20 + quant_30 + quant_40 + quant_50 + quant_60 + 
#>     quant_70 + quant_80 + quant_90
svytotal(quants_formula, m1)
#>            total     SE
#> quant_10 0.15758 0.0165
#> quant_20 0.29467 0.0206
#> quant_30 0.40922 0.0223
#> quant_40 0.52374 0.0226
#> quant_50 0.61322 0.0221
#> quant_60 0.69896 0.0208
#> quant_70 0.79488 0.0183
#> quant_80 0.85959 0.0157
#> quant_90 0.93444 0.0112
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
summary(g1*df_resp$d)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   1.300   1.494   1.790   2.049   2.391   3.173
```
