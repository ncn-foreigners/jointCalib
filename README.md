
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
- Zhang, S., Han, P., and Wu, C. (2023) [Calibration Techniques
  Encompassing Survey Sampling](https://doi.org/10.1111/insr.12518.),
  Missing Data Analysis and Causal Inference. International Statistical
  Review, 91: 165–192.

which allows to calibrate weights to known (or estimated) totals and
quantiles jointly. As an backend for calibration
[sampling](https://cran.r-project.org/web/packages/sampling)
(`sampling::calib`),
[laeken](https://cran.r-project.org/web/packages/laeken)
(`laeken::calibWeights`) or
[survey](https://cran.r-project.org/web/packages/survey/index.html)
(`survey::grake`) package can be used. One can also apply empirical
likelihood using `optim::constrOptim` function.

Currently supports:

- calibration of totals,
- calibration of quantiles,
- calibration using empirical likelihood method.

Further plans:

- generalized calibration via `sampling::gencalib`,
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
set.seed(20230817)
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
#>           x           y p        d
#> 1 33.916804 -487.869536 1 1.984127
#> 2 30.798712 -250.224643 1 1.984127
#> 3 28.778660  203.242132 1 1.984127
#> 4 27.633994  185.899405 1 1.984127
#> 7 61.903192   25.259454 1 1.984127
#> 8  6.673617    3.999977 1 1.984127
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
#>  0.6146  0.7013  0.8690  1.0000  0.9509  2.2909
```

``` r
data.frame(true = quants_known$x, est = weightedQuantile(df_resp$x, result1$g*df_resp$d, probs))
#>          true       est
#> 10%  7.078067  7.103456
#> 20% 14.831221 14.832920
#> 30% 23.146180 23.172026
#> 40% 31.641911 31.658152
#> 50% 39.033812 39.154276
#> 60% 47.527168 47.593709
#> 70% 54.984229 55.311953
#> 80% 64.073167 64.062629
#> 90% 71.565441 72.141764
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
#>  0.5676  0.7188  0.8872  1.0000  1.0043  2.3416
```

``` r
data.frame(true = quants_known$x, est = weightedQuantile(df_resp$x, result2$g*df_resp$d, probs))
#>          true       est
#> 10%  7.078067  7.103456
#> 20% 14.831221 14.832920
#> 30% 23.146180 23.085873
#> 40% 31.641911 31.658152
#> 50% 39.033812 39.154276
#> 60% 47.527168 47.593709
#> 70% 54.984229 54.954054
#> 80% 64.073167 64.062629
#> 90% 71.565441 71.293602
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
#>  0.5775  0.7165  0.8871  1.0000  1.0085  2.3353
```

``` r
data.frame(true = quants_known$x, est = weightedQuantile(df_resp$x, result3$g*df_resp$d, probs))
#>          true       est
#> 10%  7.078067  7.103456
#> 20% 14.831221 14.832920
#> 30% 23.146180 23.172026
#> 40% 31.641911 31.658152
#> 50% 39.033812 39.154276
#> 60% 47.527168 47.593709
#> 70% 54.984229 55.311953
#> 80% 64.073167 64.698127
#> 90% 71.565441 71.293602
```

Empirical likelihood method can be applied using the following code

``` r
result4a <- joint_calib(formula_quantiles = ~ x,
                        data = df_resp,
                        dweights = df_resp$d,
                        N = N,
                        pop_quantiles = quants_known,
                        method = "el",
                        backend = "optim")
summary(result4a$g)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.6176  0.7085  0.8652  1.0000  0.9547  2.2884

result4b <- joint_calib(formula_totals = ~ x,
                        formula_quantiles = ~ x,
                        data = df_resp,
                        dweights = df_resp$d,
                        N = N,
                        pop_quantiles = quants_known,
                        pop_totals = totals_known,
                        method = "el",
                        backend = "optim")

summary(result4b$g)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.6065  0.7247  0.8733  1.0000  1.0005  2.5797
```

``` r
data.frame(true = quants_known$x, 
           est = weightedQuantile(df_resp$x, result2$g*df_resp$d, probs),
           est_el1 = weightedQuantile(df_resp$x, result4a$g*df_resp$d, probs),
           est_el2 = weightedQuantile(df_resp$x, result4b$g*df_resp$d, probs))
#>          true       est   est_el1  est_el2
#> 10%  7.078067  7.103456  7.015648  6.94050
#> 20% 14.831221 14.832920 14.832920 14.83292
#> 30% 23.146180 23.085873 23.085873 23.07136
#> 40% 31.641911 31.658152 31.658152 31.40628
#> 50% 39.033812 39.154276 39.154276 38.92874
#> 60% 47.527168 47.593709 47.136099 47.13610
#> 70% 54.984229 54.954054 55.311953 55.31195
#> 80% 64.073167 64.062629 64.062629 64.69813
#> 90% 71.565441 71.293602 71.293602 72.14176
```

Finally, compare all method with true Y distribution where

- `est_y1` – weights calibrated to quantiles only,
- `est_y2` – weights calibrated to quantiles and totals,
- `est_y3` – weights calibrated to quantiles and totals using `logit`
  distance and range limitations,
- `est_y4` – weights calibrated to quantiles only using `el`.
- `est_y4` – weights calibrated to quantiles and totals using `el`.

``` r
results_y <- data.frame(true_y = y_quant_true,
                        est_y1 = weightedQuantile(df_resp$y, result1$g * df_resp$d, probs),
                        est_y2 = weightedQuantile(df_resp$y, result2$g * df_resp$d, probs),
                        est_y3 = weightedQuantile(df_resp$y, result3$g * df_resp$d, probs),
                        est_y4 = weightedQuantile(df_resp$y, result4a$g * df_resp$d, probs),
                        est_y5 = weightedQuantile(df_resp$y, result4b$g * df_resp$d, probs))

results_y
#>         true_y     est_y1     est_y2     est_y3     est_y4     est_y5
#> 10% -289.97527 -291.79309 -292.34399 -292.34399 -291.79309 -292.34399
#> 20% -146.15516 -154.36951 -146.90524 -146.90524 -154.36951 -146.90524
#> 30%  -37.31173  -25.12616  -25.02487  -25.02487  -25.12616  -25.02487
#> 40%   71.38290   95.73374   99.44012   99.44012   95.73374   99.44012
#> 50%  172.57996  194.32959  197.68361  197.68361  194.32959  197.68361
#> 60%  289.64478  323.74182  326.60335  326.60335  323.74182  328.19661
#> 70%  437.99144  473.52287  477.95125  477.95125  473.52287  477.95125
#> 80%  644.99052  653.93520  653.93520  653.93520  653.93520  654.96939
#> 90% 1223.23217 1111.48582 1111.48582 1111.48582 1111.48582 1136.33628
```

Here is the sum of squares and `el` seems to perform the best but it is
mainly due to smaller difference for 90% quantile.

``` r
apply(results_y[, -1], 2, FUN = function(x) sum((x-results_y[,1])^2))
#>   est_y1   est_y2   est_y3   est_y4   est_y5 
#> 16277.62 17104.52 17104.52 16277.62 12308.04
```

### Using `survey` package

example 1a: calibrate only quantiles (deciles)

``` r
A <- joint_calib_create_matrix(X_q = model.matrix(~x+0, df_resp), N = N, pop_quantiles = quants_known)
colnames(A) <- paste0("quant_", gsub("\\D", "", names(quants_known$x)))
A <- as.data.frame(A)
df_resp <- cbind(df_resp, A)
head(df_resp)
#>           x           y p        d quant_10 quant_20 quant_30 quant_40 quant_50
#> 1 33.916804 -487.869536 1 1.984127    0.000    0.000    0.000    0.000    0.001
#> 2 30.798712 -250.224643 1 1.984127    0.000    0.000    0.000    0.001    0.001
#> 3 28.778660  203.242132 1 1.984127    0.000    0.000    0.000    0.001    0.001
#> 4 27.633994  185.899405 1 1.984127    0.000    0.000    0.000    0.001    0.001
#> 7 61.903192   25.259454 1 1.984127    0.000    0.000    0.000    0.000    0.000
#> 8  6.673617    3.999977 1 1.984127    0.001    0.001    0.001    0.001    0.001
#>   quant_60 quant_70 quant_80 quant_90
#> 1    0.001    0.001    0.001    0.001
#> 2    0.001    0.001    0.001    0.001
#> 3    0.001    0.001    0.001    0.001
#> 4    0.001    0.001    0.001    0.001
#> 7    0.000    0.000    0.001    0.001
#> 8    0.001    0.001    0.001    0.001
m1 <- svydesign(ids = ~1, data = df_resp, weights = ~d)
quants_formula <- as.formula(paste("~", paste(colnames(A), collapse = "+")))
quants_formula
#> ~quant_10 + quant_20 + quant_30 + quant_40 + quant_50 + quant_60 + 
#>     quant_70 + quant_80 + quant_90
svytotal(quants_formula, m1)
#>            total     SE
#> quant_10 0.16270 0.0165
#> quant_20 0.28602 0.0201
#> quant_30 0.42857 0.0221
#> quant_40 0.54365 0.0222
#> quant_50 0.64881 0.0213
#> quant_60 0.75397 0.0192
#> quant_70 0.82738 0.0169
#> quant_80 0.90873 0.0128
#> quant_90 0.95238 0.0095
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
#>  0.6146  0.7013  0.8690  1.0000  0.9509  2.2909
```
