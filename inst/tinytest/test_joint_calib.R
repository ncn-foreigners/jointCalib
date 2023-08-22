## create data
library(survey)
set.seed(123)
N <- 1000
x <- runif(N, 0, 80)
y <- 1000+10*x+rnorm(N, 0, 300)
p <- rbinom(N, 1, prob = (1.2+0.024*x)^-1)
probs <- seq(0.1, 0.9, 0.1)
quants_known <- list(x=quantile(x, probs))
totals_known <- c(x=sum(x))
df <- data.frame(x, y, p)
df_resp <- df[df$p == 1, ]
df_resp$d <- N/nrow(df_resp)
svydes <- svydesign(ids=~1, data=df, weights = ~1)
quants_known_svyq <- svyquantile(~ x, svydes, probs)

# example 1 -- only quantiles ---------------------------------------------

expect_silent(
  result1 <- joint_calib(formula_quantiles = ~x,
                         data = df_resp,
                         dweights = df_resp$d,
                         N = N,
                         pop_quantiles = quants_known,
                         method = "linear",
                         backend = "sampling")

)

## check if sum up to N
expect_equal(
  sum(result1$g*df_resp$d),
  N
)

## check the difference with quantiles

expect_equal(
  colSums(result1$Xs*result1$g*df_resp$d)[-1],
  probs
)

## check with quantiles
expect_true(
  {
    quants_estimated <- laeken::weightedQuantile(x = df_resp$x, probs = probs, weights = result1$g*df_resp$d)
    sum(( quants_estimated - quants_known$x)^2)/length(quants_known$x) < 0.05
  }
)

## check for nonnegative weights
expect_true(
  all(result1$g > 0)
)


# example 2 - quantiles and totals ----------------------------------------

expect_silent(
  result2 <- joint_calib(formula_quantiles = ~x,
                         formula_totals = ~x,
                         data = df_resp,
                         dweights =  df_resp$d,
                         N = N,
                         pop_quantiles = quants_known,
                         pop_totals = totals_known,
                         method = "linear",
                         backend = "sampling")

)

## check if sum up to N
expect_equal(
  sum(result2$g*df_resp$d),
  N
)

## check the difference with quantiles

expect_equal(
  unname(colSums(result2$Xs*result2$g*df_resp$d)[2:10]),
  probs
)

## check with quantiles
expect_true(
  {
    quants_estimated <- laeken::weightedQuantile(x = df_resp$x, probs = probs, weights = result2$g*df_resp$d)
    sum(( quants_estimated - quants_known$x)^2)/length(quants_known$x) < 0.05
  }
)

## check for nonnegative weights
expect_true(
  all(result2$g > 0)
)

# example 3 - empirical likelihood ----------------------------------------

expect_silent(
  result3 <- joint_calib(formula_quantiles = ~x,
                         formula_totals = ~x,
                         data = df_resp,
                         dweights =  df_resp$d,
                         N = N,
                         pop_quantiles = quants_known,
                         pop_totals = totals_known,
                         method = "el")

)

## check if sum up to N
expect_equal(
  sum(result3$g*df_resp$d),
  N
)

## check the difference with quantiles

expect_equal(
  unname(colSums(result3$Xs*result3$g*df_resp$d)[2:10]),
  probs,
  tolerance = 1e-3
)

## check with quantiles
expect_true(
  {
    quants_estimated <- laeken::weightedQuantile(x = df_resp$x, probs = probs, weights = result3$g*df_resp$d)
    sum(( quants_estimated - quants_known$x)^2)/length(quants_known$x) < 0.05
  }
)

## check for nonnegative weights
expect_true(
  all(result3$g > 0)
)

# example 4 - ebal ----------------------------------------

expect_silent(
  result4 <- joint_calib(formula_quantiles = ~x,
                         formula_totals = ~x,
                         data = df_resp,
                         dweights =  df_resp$d,
                         N = N,
                         pop_quantiles = quants_known,
                         pop_totals = totals_known,
                         method = "eb")

)

## check if sum up to N
expect_equal(
  sum(result4$g*df_resp$d),
  N,
  tolerance = 0.01
)

## check the difference with quantiles

expect_equal(
  unname(colSums(result4$Xs*result4$g*df_resp$d)[2:10]),
  probs,
  tolerance = 1e-3
)

## check with quantiles
expect_true(
  {
    quants_estimated <- laeken::weightedQuantile(x = df_resp$x, probs = probs, weights = result4$g*df_resp$d)
    sum(( quants_estimated - quants_known$x)^2)/length(quants_known$x) < 0.05
  }
)

## check for nonnegative weights
expect_true(
  all(result4$g > 0)
)


# test svyquantile --------------------------------------------------------

expect_silent(
  result5 <- joint_calib(formula_quantiles = ~x,
                         data = df_resp,
                         dweights = df_resp$d,
                         N = N,
                         pop_quantiles = quants_known_svyq,
                         method = "linear",
                         backend = "sampling")

)
