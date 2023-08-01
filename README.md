
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
#>            x           y p        d
#> 2  18.131074 -183.410128 1 2.242152
#> 3  19.903402  -22.601732 1 2.242152
#> 4  37.050626 -280.201668 1 2.242152
#> 6  46.785617   -8.721672 1 2.242152
#> 10 14.447789  238.925141 1 2.242152
#> 11  8.625689  206.695376 1 2.242152
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
#>  0.6110  0.6969  0.8576  1.0000  1.1436  1.7840
```

``` r
data.frame(true = quants_known$x, est = weightedQuantile(df_resp$x, result1$g*df_resp$d, probs))
#>          true      est
#> 10%  9.995078 10.01146
#> 20% 17.674664 17.69873
#> 30% 25.888919 25.88402
#> 40% 34.145194 34.16041
#> 50% 42.696666 42.96362
#> 60% 51.238375 51.21526
#> 70% 58.200756 58.23054
#> 80% 65.416836 65.95623
#> 90% 72.673909 72.74871
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
#>  0.5702  0.7011  0.8438  1.0000  1.1703  1.8168
```

``` r
data.frame(true = quants_known$x, est = weightedQuantile(df_resp$x, result2$g*df_resp$d, probs))
#>          true      est
#> 10%  9.995078 10.01146
#> 20% 17.674664 17.69873
#> 30% 25.888919 25.88402
#> 40% 34.145194 34.16041
#> 50% 42.696666 42.96362
#> 60% 51.238375 51.21526
#> 70% 58.200756 58.23054
#> 80% 65.416836 65.95623
#> 90% 72.673909 72.74871
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
#>  0.5718  0.7011  0.8420  1.0000  1.1735  1.7981
```

``` r
data.frame(true = quants_known$x, est = weightedQuantile(df_resp$x, result3$g*df_resp$d, probs))
#>          true       est
#> 10%  9.995078  9.772978
#> 20% 17.674664 17.434639
#> 30% 25.888919 25.884016
#> 40% 34.145194 34.122365
#> 50% 42.696666 42.477097
#> 60% 51.238375 51.215263
#> 70% 58.200756 57.754228
#> 80% 65.416836 64.904944
#> 90% 72.673909 72.385134
```

### Using `survey` package

example 1a: calibrate only quantiles (deciles)

``` r
A <- joint_calib_create_matrix(X_q = model.matrix(~x+0, df_resp), N = N, pop_quantiles = quants_known)
colnames(A) <- paste0("quant_", gsub("\\D", "", names(quants_known$x)))
A <- as.data.frame(A)
df_resp <- cbind(df_resp, A)
head(df_resp)
#>            x           y p        d quant_10      quant_20 quant_30 quant_40
#> 2  18.131074 -183.410128 1 2.242152    0.000 6.073619e-202    0.001    0.001
#> 3  19.903402  -22.601732 1 2.242152    0.000  0.000000e+00    0.001    0.001
#> 4  37.050626 -280.201668 1 2.242152    0.000  0.000000e+00    0.000    0.000
#> 6  46.785617   -8.721672 1 2.242152    0.000  0.000000e+00    0.000    0.000
#> 10 14.447789  238.925141 1 2.242152    0.000  1.000000e-03    0.001    0.001
#> 11  8.625689  206.695376 1 2.242152    0.001  1.000000e-03    0.001    0.001
#>    quant_50 quant_60 quant_70 quant_80 quant_90
#> 2     0.001    0.001    0.001    0.001    0.001
#> 3     0.001    0.001    0.001    0.001    0.001
#> 4     0.001    0.001    0.001    0.001    0.001
#> 6     0.000    0.001    0.001    0.001    0.001
#> 10    0.001    0.001    0.001    0.001    0.001
#> 11    0.001    0.001    0.001    0.001    0.001
m1 <- svydesign(ids = ~1, data = df_resp, weights = ~d)
quants_formula <- as.formula(paste("~", paste(colnames(A), collapse = "+")))
quants_formula
#> ~quant_10 + quant_20 + quant_30 + quant_40 + quant_50 + quant_60 + 
#>     quant_70 + quant_80 + quant_90
svytotal(quants_formula, m1)
#>            total     SE
#> quant_10 0.16368 0.0175
#> quant_20 0.30717 0.0219
#> quant_30 0.42375 0.0234
#> quant_40 0.54036 0.0236
#> quant_50 0.67040 0.0223
#> quant_60 0.75785 0.0203
#> quant_70 0.81390 0.0184
#> quant_80 0.88341 0.0152
#> quant_90 0.94395 0.0109
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
#>  0.6110  0.6969  0.8576  1.0000  1.1436  1.7840
```
