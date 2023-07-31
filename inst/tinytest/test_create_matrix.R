## first example
set.seed(123)
N <- 1000
x <- as.matrix(rnorm(N))
quants <- list(quantile(x, c(0.25,0.5,0.75)))
A <- joint_calib_create_matrix(x, N, quants)

expect_equal(
  dim(A),
  c(1000, 3)
)

expect_equal(
  colSums(A),
  c(0.25, 0.50, 0.75),
  tolerance = 0.01
)

## second example
set.seed(123)
x1 <- rnorm(N)
x2 <- rchisq(N, 1)
x <- cbind(x1, x2)
quants <- list(quantile(x1, 0.5), quantile(x2, c(0.1, 0.75, 0.9)))
B <- joint_calib_create_matrix(x, N, quants)

expect_equal(
  dim(B),
  c(1000, 4)
)

expect_equal(
  colSums(B),
  c(0.50, 0.1, 0.75, 0.9),
  tolerance = 0.01
)

