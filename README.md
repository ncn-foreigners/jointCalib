
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
#>           x          y p        d
#> 1 38.859020 -201.76413 1 1.912046
#> 2 30.084190  614.18087 1 1.912046
#> 3 48.576182  157.16099 1 1.912046
#> 5 57.292801  456.68291 1 1.912046
#> 6  2.989587  268.84749 1 1.912046
#> 8 45.070184   22.00903 1 1.912046
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
#>  0.6705  0.6792  0.8574  1.0000  1.2163  1.8034
```

``` r
data.frame(true = quants_known$x, est = weightedQuantile(df_resp$x, result1$g*df_resp$d, probs))
#>         true       est
#> 10%  6.58857  6.750345
#> 20% 13.19709 13.220486
#> 30% 21.20933 21.687725
#> 40% 29.75575 29.949604
#> 50% 37.28935 37.359518
#> 60% 44.85285 44.845606
#> 70% 53.02708 53.015611
#> 80% 61.06943 61.574482
#> 90% 70.64060 70.689237
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
#>  0.6452  0.6919  0.8666  1.0000  1.2155  1.8404
```

``` r
data.frame(true = quants_known$x, est = weightedQuantile(df_resp$x, result2$g*df_resp$d, probs))
#>         true       est
#> 10%  6.58857  6.553297
#> 20% 13.19709 13.220486
#> 30% 21.20933 20.885609
#> 40% 29.75575 29.515912
#> 50% 37.28935 37.212618
#> 60% 44.85285 44.845606
#> 70% 53.02708 53.015611
#> 80% 61.06943 60.856741
#> 90% 70.64060 70.610593
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
#>  0.6452  0.6921  0.8675  1.0000  1.2166  1.8179
```

``` r
data.frame(true = quants_known$x, est = weightedQuantile(df_resp$x, result3$g*df_resp$d, probs))
#>         true       est
#> 10%  6.58857  6.553297
#> 20% 13.19709 13.103485
#> 30% 21.20933 20.885609
#> 40% 29.75575 29.515912
#> 50% 37.28935 37.212618
#> 60% 44.85285 44.845606
#> 70% 53.02708 53.015611
#> 80% 61.06943 60.856741
#> 90% 70.64060 70.610593
```

### Using `survey` package

example 1a: calibrate only quantiles (deciles)

``` r
A <- joint_calib_create_matrix(X_q = model.matrix(~x+0, df_resp), N = N, pop_quantiles = quants_known)
colnames(A) <- paste0("quant_", gsub("\\D", "", names(quants_known$x)))
A <- as.data.frame(A)
df_resp <- cbind(df_resp, A)
head(df_resp)
#>           x          y p        d quant_10 quant_20 quant_30      quant_40
#> 1 38.859020 -201.76413 1 1.912046    0.000    0.000    0.000  0.000000e+00
#> 2 30.084190  614.18087 1 1.912046    0.000    0.000    0.000 2.281127e-146
#> 3 48.576182  157.16099 1 1.912046    0.000    0.000    0.000  0.000000e+00
#> 5 57.292801  456.68291 1 1.912046    0.000    0.000    0.000  0.000000e+00
#> 6  2.989587  268.84749 1 1.912046    0.001    0.001    0.001  1.000000e-03
#> 8 45.070184   22.00903 1 1.912046    0.000    0.000    0.000  0.000000e+00
#>   quant_50    quant_60 quant_70 quant_80 quant_90
#> 1    0.000 1.00000e-03    0.001    0.001    0.001
#> 2    0.001 1.00000e-03    0.001    0.001    0.001
#> 3    0.000 0.00000e+00    0.001    0.001    0.001
#> 5    0.000 0.00000e+00    0.000    0.001    0.001
#> 6    0.001 1.00000e-03    0.001    0.001    0.001
#> 8    0.000 4.11706e-98    0.001    0.001    0.001
m1 <- svydesign(ids = ~1, data = df_resp, weights = ~d)
quants_formula <- as.formula(paste("~", paste(colnames(A), collapse = "+")))
quants_formula
#> ~quant_10 + quant_20 + quant_30 + quant_40 + quant_50 + quant_60 + 
#>     quant_70 + quant_80 + quant_90
svytotal(quants_formula, m1)
#>            total     SE
#> quant_10 0.14914 0.0156
#> quant_20 0.29637 0.0200
#> quant_30 0.42065 0.0216
#> quant_40 0.51625 0.0219
#> quant_50 0.59847 0.0215
#> quant_60 0.71510 0.0198
#> quant_70 0.79924 0.0175
#> quant_80 0.87572 0.0144
#> quant_90 0.94455 0.0100
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
#>  0.6705  0.6792  0.8574  1.0000  1.2163  1.8034
```
