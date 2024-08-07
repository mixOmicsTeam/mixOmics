context("tpca")

test_that(
  "basic tpca sense checks",
  code = {
    test_tensor <- array(1:24, dim = c(2, 4, 3))
    tpca_obj <- tpca(test_tensor)
    expect_equal(length(tpca_obj$explained_variance), tpca_obj$ncomp)
  }
)
