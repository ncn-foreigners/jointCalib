
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![R-CMD-check](https://github.com/ncn-foreigners/jointCalib/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ncn-foreigners/jointCalib/actions/workflows/R-CMD-check.yaml)
[![CRAN-s](https://www.r-pkg.org/badges/version/jointCalib)](https://CRAN.R-project.org/package=jointCalib)
[![CRAN-d](http://cranlogs.r-pkg.org/badges/grand-total/jointCalib?color=blue)](https://cran.r-project.org/package=jointCalib)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8355993.svg)](https://doi.org/10.5281/zenodo.8355993)
[![Dependencies](https://tinyverse.netlify.com/badge/jointCalib)](https://cran.r-project.org/package=jointCalib)
[![Mentioned in Awesome Official
Statistics](https://awesome.re/mentioned-badge.svg)](http://www.awesomeofficialstatistics.org)

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
  entropy balancing method,
- covariate distribution entropy balancing for ATT and QTT
  (distributional entropy balancing; DEB),
- covariate distribution balancing propensity score for ATE and QTE
  (distributional propensity score; DPS).

Further plans:

- generalized calibration via `sampling::gencalib`,
- calibration for Gini and other metrics,
- …

## Funding

Work on this package is supported by the the National Science Centre,
OPUS 22 grant no. 2020/39/B/HS4/00941.

## Installation

You can install CRAN version of the package using

``` r
install.packages("jointCalib")
```

You can install the development version of `jointCalib` from GitHub
with:

``` r
# install.packages("remotes")
remotes::install_github("ncn-foreigners/jointCalib")
```
