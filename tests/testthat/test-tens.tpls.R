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
    # out a different compressed matrix
    n <- 4
    p <- 5
    t <- 3
    ncomp_input <- 2

    # this test fails if `test_tensor` is rank deficient, as is the case if we
    # use array(1:24, dim = c(3, 4, 2))
    set.seed(1)
    test_tensor <- array(rnorm(n * p * t, mean = 0, sd = 5), dim = c(n, p, t))

    # just use the identity as m transform to get over problem 1
    transforms <- matrix_to_m_transforms(diag(t))

    tpca_obj <- tpca(
      test_tensor,
      ncomp = ncomp_input,
      m = transforms$m,
      minv = transforms$minv,
      # return full tensors to get over problem 2
      matrix_output = FALSE
    )

    tpls_obj <- tpls(
      test_tensor, test_tensor,
      ncomp = ncomp_input,
      mode = "tsvdm",
      m = transforms$m,
      minv = transforms$minv,
      # return full tensors to get over problem 2
      matrix_output = FALSE
    )

    expect_equal(tpca_obj$loadings, tpls_obj$y_loadings)
    expect_equal(tpca_obj$variates, tpls_obj$y_projected)
  }
)
