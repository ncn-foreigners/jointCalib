# generate data

set.seed(123)
N <- 1500
X1 <- rnorm(N)
X2 <- rnorm(N)
X3 <- rbinom(N, size = 1, prob = .5)
X1X3 <- X1*X3
D_star <- 0.5*X1 + 0.3*X2 + 0.2*X1*X2 - 0.5*X1*X3 - 1
D <- ifelse(D_star > rnorm(N), 1, 0) # Treatment indicator
y <- 0.5*D + X1 + X2 + X2*X3 + rnorm(N) # Outcome
dat <- data.frame(D = D, X1 = X1, X2 = X2, X3 = X3, X1X3=X1X3, Y = y)

## general tests

expect_silent(
  joint_calib_att(
    formula_means = ~ X1 + X2 + X3,
    formula_quantiles = ~ X1 + X2,
    treatment = ~ D,
    data = dat,
    method = "raking"
  )
)

expect_silent(
  joint_calib_att(
    formula_means = ~ X1 + X2 + X3,
    formula_quantiles = ~ X1 + X2,
    treatment = ~ D,
    data = dat,
    method = "el"
  )
)

## test list of probs

expect_silent(
  joint_calib_att(
    formula_means = ~ X1 + X2 + X3,
    formula_quantiles = ~ X1 + X2,
    treatment = ~ D,
    data = dat,
    probs = list(X1=0.5, X2 = c(0.25, 0.5)),
    method = "raking"
  )
)

## test error if repeated

expect_error(
  joint_calib_att(
    formula_means = ~ X1 + X2 + X3,
    formula_quantiles = ~ X1 + X2,
    treatment = ~ D,
    data = dat,
    probs = list(X1=0.5, X2 = c(0.5, 0.5)),
    method = "raking"
  )
)

## test for empty list

expect_error(
  joint_calib_att(
    formula_means = ~ X1 + X2 + X3,
    formula_quantiles = ~ X1 + X2,
    treatment = ~ D,
    data = dat,
    probs = list(X1=numeric()),
    method = "raking"
  )
)

