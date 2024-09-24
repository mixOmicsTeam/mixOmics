context("tpls")

#' @description Unnames and ensures the top row in a matrix-type output is all
#' positive. This allows for comparison between pls results that are the same
#' but just differ by a negative sign due to svd solver or other implementation
#' detail.
.make_signs_consistent <- function(mat) {
  mat <- unname(mat)
  ncols <- dim(mat)[2]
  for (i in seq_len(ncols)) {
    if (mat[1, i] < 0) {
      mat[, i] <- -mat[, i]
    }
  }
  return(mat)
}

test_that(
  "tpls: tsvdm-tpls mode agrees with tpca",
  code = {
    # it actually takes a lot of coercing to make these two outputs comparable
    # because their default configurations are aimed at completely different
    # problems. that being said, this test here is probably still useful as a
    # sense check.
    #
    # problem 1: squaring the singular values means that the k_t_flatten_sort
    # compression is different and picks out a different compressed matrix. Sort
    # order of the squared values is not necessarily the same.
    #
    # problem 2: loadings (and therefore projections) may differ in sign, think
    # this is dependent on the svd solver implementation?
    #
    # problem 3: if you use centering (on for tpca and tpls by default), the
    # faces of each tensor loses one rank, so only rank - 1 columns of each face
    # in the loadings tensor ("v") will match.
    #
    # problem 4: tpls computes the svd on XtX, which is rank-deficient if X is
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

    transforms <- dctii_m_transforms(t)
    m <- transforms$m
    minv <- transforms$minv

    tpca_obj <- tpca(
      test_tensor,
      ncomp = ncomp_input,
      m = m,
      minv = minv,
      # turn off centering to get over problem 3
      center = FALSE,
      # return full tensors to get over problem 1
      matrix_output = FALSE
    )

    tpls_obj <- tpls(
      test_tensor, test_tensor,
      ncomp = ncomp_input,
      m = m,
      minv = minv,
      mode = "tsvdm",
      # turn off centering to get over problem 3
      center = FALSE,
      # return full tensors to get over problem 1
      matrix_output = FALSE
    )

    # compare elementwise absolute values to get over problem 2
    expect_equal(
      abs(tpca_obj$loadings),
      # only compare k columns to get over problem 4
      abs(tpls_obj$y_loadings[, 1:k, ])
    )

    expect_equal(
      abs(tpca_obj$variates),
      abs(tpls_obj$y_projected)
    )
  }
)

test_that(
  "tpls agrees with MixOmics pls",
  code = {
    n <- 4
    p <- 5
    q <- 7
    t <- 1
    k <- min(n, p, q)
    ncomp_input <- 2

    set.seed(1)
    test_x <- array(rnorm(n * p * t, mean = 0, sd = 5), dim = c(n, p, t))
    test_y <- array(rnorm(n * q * t, mean = 0, sd = 3), dim = c(n, q, t))

    mixomics_pls <- pls(
      test_x[, , 1],
      test_y[, , 1],
      ncomp = ncomp_input,
      scale = FALSE,
      mode = "canonical"
    )

    transforms <- matrix_to_m_transforms(diag(1))
    tensor_pls <- tpls(
      test_x,
      test_y,
      ncomp = ncomp_input,
      m = transforms$m,
      minv = transforms$minv,
      mode = "canonical"
    )

    expect_equal(
      .make_signs_consistent(mixomics_pls$loadings$X),
      .make_signs_consistent(tensor_pls$x_loadings)
    )

    expect_equal(
      .make_signs_consistent(mixomics_pls$loadings$Y),
      .make_signs_consistent(tensor_pls$y_loadings)
    )

    expect_equal(
      .make_signs_consistent(mixomics_pls$variates$X),
      .make_signs_consistent(tensor_pls$x_projected)
    )

    expect_equal(
      .make_signs_consistent(mixomics_pls$variates$Y),
      .make_signs_consistent(tensor_pls$y_projected)
    )
  }
)
