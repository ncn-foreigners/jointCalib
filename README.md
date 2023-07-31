# Overview

A small package for joint calibration of totals and quantiles. The package combines the following approaches:

+ Deville, J. C., and Särndal, C. E. (1992). [Calibration estimators in survey sampling](https://www.tandfonline.com/doi/abs/10.1080/01621459.1992.10475217). Journal of the American statistical Association, 87(418), 376-382.
+ Harms, T. and Duchesne, P. (2006). [On calibration estimation for quantiles](https://www150.statcan.gc.ca/n1/pub/12-001-x/2006001/article/9255-eng.pdf). Survey Methodology, 32(1), 37.

which allows to calibrate weights to known (or estimated) totals and quantiles jointly. As an backend for calibration [sampling](https://cran.r-project.org/web/packages/sampling) (`sampling::calib`) or [laeken](https://cran.r-project.org/web/packages/laeken) (`laeken::calibWeights`) package can be used.

## Installation

You can install the development version of `jointCalib` from GitHub with:

```r
# install.packages("remotes")
remotes::install_github("ncn-foreigners/jointCalib")
```

## Examples 

### Example 1 

Based on Haziza, D., and Lesage, É. (2016). A discussion of weighting procedures for unit nonresponse. Journal of Official Statistics, 32(1), 129-145.

```r
library(jointCalib)
```

```r
N <- 1000
x <- runif(N, 0, 80)
y <- exp(-0.1 + 0.1*x) + rnorm(N, 0, 300)
p <- rbinom(N, 1, prob = exp(-0.2 - 0.014*x))
probs <- seq(0.1, 0.9, 0.1)
quants_known <- quantile(x, probs)
totals_known <- sum(x)
df <- data.frame(x, y, p)
df_resp <- df[df$p == 1, ]
df_resp$d <- N/nrow(df_resp)
y_quant_true <- quantile(y, probs)
```

example 1a: calibrate only quantiles (deciles)

```r
result1 <- joint_calib(X_q = as.matrix(df_resp$x),
                      d = df_resp$d,
                      N = N,
                      pop_quantiles = list(quants_known),
                      method = "linear",
                      backend = "sampling")
```

example 1b: calibrate with quantiles (deciles) and totals

```r
result2 <- joint_calib(X_q = as.matrix(df_resp$x),
                       X = as.matrix(df_resp$x),
                       d = df_resp$d,
                       N = N,
                       pop_quantiles = list(quants_known),
                       pop_totals = totals_known,
                       method = "linear",
                       backend = "sampling")
```


## Funding

Work on this package is supported by the the National Science Center, OPUS 22 grant no. 2020/39/B/HS4/00941.







