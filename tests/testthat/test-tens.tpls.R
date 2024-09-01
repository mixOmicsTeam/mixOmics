context("tpls")

test_that(
  "tpls: tsvdm-tpls mode agrees with tpca",
  code = {
    # it actually takes a lot of coercing to make these two outputs comparable
    # because their default configurations are aimed at completely different
    # problems. that being said, this test here is probably still useful as a
    # sense check.
    #
    # problem 1: we perform the tsvdm decomposition in the hat-space ('m' space)
    # so the tsvdm of ft(m(x)) %fp% m(x) "XtX" really has nothing to do with the
    # tsvdm decomposition of m(X) "X" directly anymore.
    #
    # problem 2: even if it did (by using m == identity), squaring the singular
    # values means that the k_t_flatten_sort compression is different and picks
    # out a different compressed matrix.
    #
    # problem 3: loadings (and therefore projections) may differ in sign, think
    # this is dependent on the svd solver implementation?
    #
    # problem 4: if you use centering (on for tpca and tpls by default), the
    # faces of each tensor loses one rank, so only rank - 1 columns of each face
    # in the loadings tensor ("v") will match.
    #
    # problem 5: tpls computes the svd on XtX, which is rank-deficient if X is
    # non-square. this means that the last few columns of "v" are absolutely
    # crap and should not be compared. these columns may not even exist in
    # tpca depending on the input dimensions.
    n <- 4
    p <- 5
    t <- 3
    k <- min(n, p)
    ncomp_input <- 2

    # this test fails if `test_tensor` is rank deficient, as is the case if we
    # use array(1:24, dim = c(3, 4, 2))
    set.seed(1)
    test_tensor <- array(rnorm(n * p * t, mean = 0, sd = 5), dim = c(n, p, t))

    # just use the identity as m transform to get over problem 1
    transforms <- matrix_to_m_transforms(diag(t))
    m <- transforms$m
    minv <- transforms$minv

    tpca_obj <- tpca(
      test_tensor,
      ncomp = ncomp_input,
      m = m,
      minv = minv,
      # turn off centering to get over problem 4
      center = FALSE,
      # return full tensors to get over problem 2
      matrix_output = FALSE
    )

    tpls_obj <- tpls(
      test_tensor, test_tensor,
      ncomp = ncomp_input,
      m = m,
      minv = minv,
      mode = "tsvdm",
      # turn off centering to get over problem 4
      center = FALSE,
      # return full tensors to get over problem 2
      matrix_output = FALSE
    )

    # compare elementwise absolute values to get over problem 3
    expect_equal(
      abs(tpca_obj$loadings),
      # only compare k columns to get over problem 5
      abs(tpls_obj$y_loadings[, 1:k, ])
    )

    expect_equal(
      abs(tpca_obj$variates),
      abs(tpls_obj$y_projected)
    )
  }
)
