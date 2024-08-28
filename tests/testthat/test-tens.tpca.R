context("tpca")

test_that(
  "basic tpca sense checks",
  code = {
    n <- 2
    p <- 4
    t <- 3
    ncomp_input <- 1
    test_tensor <- array(1:24, dim = c(n, p, t))
    tpca_obj <- tpca(test_tensor, ncomp = ncomp_input)
    expect_equal(length(tpca_obj$explained_variance), ncomp_input)
    expect_equal(dim(tpca_obj$variates), c(n, ncomp_input))
    expect_equal(dim(tpca_obj$loadings), c(p, ncomp_input))
  }
)
