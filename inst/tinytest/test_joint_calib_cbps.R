# generate data
set.seed(123)
N <- 1500
X1 <- rnorm(N)
X2 <- rnorm(N)
X3 <- rbinom(N, size = 1, prob = .5)
X1X3 <- X1*X3
D_star <- 0.5*X1 + 0.3*X2 + 0.2*X1*X2 - 0.5*X1*X3 - 1
D <- ifelse(D_star > rnorm(N), 1, 0)
y <- 0.5*D + X1 + X2 + X2*X3 + rnorm(N)
dat <- data.frame(D = D, X1 = X1, X2 = X2, X3 = X3, X1X3=X1X3, Y = y)

## general tests

expect_silent(
  result <- joint_calib_cbps(formula_means = ~ X1 + X2 + X3,
                             formula_quantiles = ~ X1 + X2,
                             treatment = ~ D,
                             data = dat)
)

expect_equal(
  as.numeric(result$coefficients[,1]),
  c(-3.33937327206783, 0.961941047142501, 1.16056056568591, -0.0403507231216043,
    151.984307944511, 24.8020507938985, 129.064924158176, 111.06168185942,
    78.7900124714788, 161.701095459623)
)

## expect class CBPS
expect_inherits(
  joint_calib_cbps(formula_means = ~ X1 + X2 + X3,
                   formula_quantiles = ~ X1 + X2,
                   treatment = ~ D,
                   data = dat),
  "CBPS"
)

## expect error for X1:X3

expect_error(
  joint_calib_cbps(formula_means = ~ X1 + X2 + X3,
                   formula_quantiles = ~ X1 + X2 + X1:X3,
                   treatment = ~ D,
                   data = dat)
)

## check one x for quantiles
expect_silent(
  result <- joint_calib_cbps(formula_means = ~ X1,
                             formula_quantiles = ~ X1,
                             treatment = ~ D,
                             data = dat)
)


## test variable selection

expect_error(
  joint_calib_cbps(formula_means = ~ X1 + X2,
                   formula_quantiles = ~ X1 + X2,
                   treatment = ~ D,
                   data = dat,
                   variable_selection = TRUE)
)

### this takes close to 10 sec
expect_silent(
  res <- joint_calib_cbps(formula_means = ~ X1 + X2,
                          formula_quantiles = ~ X1 + X2,
                          treatment = ~ D,
                          data = dat,
                          variable_selection = TRUE,
                          target = ~ Y)
)




