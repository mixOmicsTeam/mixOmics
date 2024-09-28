context("tplsda")

test_that(
  "tplsda agrees with mixOmics pls",
  cpde = {
    n <- 4
    p <- 5
    q <- 7
    t <- 1
    k <- min(n, p, q)
    ncomp_input <- 2

    set.seed(1)
    test_x <- array(rnorm(n * p * t, mean = 0, sd = 5), dim = c(n, p, t))
    test_y <- array(rnorm(n * q * t, mean = 0, sd = 3), dim = c(n, q, t))
  }
)