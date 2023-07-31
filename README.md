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

## Examples -- TBA

## Funding

Work on this package is supported by the the National Science Center, OPUS 22 grant no. 2020/39/B/HS4/00941.







