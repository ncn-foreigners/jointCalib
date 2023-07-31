## create data
set.seed(123)
N <- 1000
x <- runif(N, 0, 80)
y <- 1000+10*x+rnorm(N, 0, 300)
p <- rbinom(N, 1, prob = (1.2+0.024*x)^-1)
probs <- seq(0.1, 0.9, 0.1)
quants_known <- quantile(x, probs)
totals_known <- sum(x)
df <- data.frame(x, y, p)
df_resp <- df[df$p == 1, ]
df_resp$d <- N/nrow(df_resp)


# example 1 -- only quantiles ---------------------------------------------

expect_silent(
  result1 <- joint_calib(X_q = as.matrix(df_resp$x),
                         d = df_resp$d,
                         N = N,
                         pop_quantiles = list(quants_known),
                         method = "linear",
                         backend = "sampling")

)

## check if sum up to N
expect_equal(
  sum(result1$w),
  N
)

## check the difference with quantiles

expect_equal(
  colSums(result1$Xs*result1$w)[-1],
  probs
)

## check with quantiles
expect_true(
  {
    quants_estimated <- laeken::weightedQuantile(x = df_resp$x, probs = probs, weights = result1$w)
    sum(( quants_estimated - quants_known)^2)/length(quants_known) < 0.05
  }
)

## check for nonnegative weights
expect_true(
  all(result$w > 0)
)


# example 2 - quantiles and totals ----------------------------------------

expect_silent(
  result2 <- joint_calib(X_q = as.matrix(df_resp$x),
                         X = as.matrix(df_resp$x),
                         d = df_resp$d,
                         N = N,
                         pop_quantiles = list(quants_known),
                         pop_totals = totals_known,
                         method = "linear",
                         backend = "sampling")

)

## check if sum up to N
expect_equal(
  sum(result2$w),
  N
)

## check the difference with quantiles

expect_equal(
  colSums(result2$Xs*result2$w)[2:10],
  probs
)

## check with quantiles
expect_true(
  {
    quants_estimated <- laeken::weightedQuantile(x = df_resp$x, probs = probs, weights = result2$w)
    sum(( quants_estimated - quants_known)^2)/length(quants_known) < 0.05
  }
)

## check for nonnegative weights
expect_true(
  all(result$w > 0)
)
