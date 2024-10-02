context("tplsda")

#' Unnames and ensures the top row in a matrix-type output is all positive. This
#' allows for comparison between pls results that are the same but just differ
#' by a negative sign due to svd solver or other implementation detail.
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
  "tplsda agrees with mixOmics plsda",
  code = {
    n <- 6
    p <- 7
    q <- 9
    t <- 1
    k <- min(n, p, q)
    ncomp_input <- 2

    set.seed(1)
    test_x <- array(rnorm(n * p * t, mean = 0, sd = 5), dim = c(n, p, t))
    test_y <- letters[1:n]

    mixomics_plsda <- plsda(
      test_x[, , 1],
      test_y,
      ncomp = ncomp_input,
      scale = FALSE
    )

    # bltodo: we do not even need to specify identity transforms here right?
    transforms <- matrix_to_m_transforms(diag(1))
    tensor_plsda <- tplsda(
      test_x,
      test_y,
      ncomp = ncomp_input,
      m = transforms$m,
      minv = transforms$minv
    )

    expect_equal(
      .make_signs_consistent(mixomics_plsda$loadings$X),
      .make_signs_consistent(tensor_plsda$x_loadings)
    )

    expect_equal(
      .make_signs_consistent(mixomics_plsda$loadings$Y),
      .make_signs_consistent(tensor_plsda$y_loadings)
    )

    expect_equal(
      .make_signs_consistent(mixomics_plsda$variates$X),
      .make_signs_consistent(tensor_plsda$x_projected)
    )

    expect_equal(
      .make_signs_consistent(mixomics_plsda$variates$Y),
      .make_signs_consistent(tensor_plsda$y_projected)
    )
  }
)
