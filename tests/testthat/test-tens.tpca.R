context("tpca")

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
  "basic tpca sense checks",
  code = {
    n <- 4
    p <- 5
    t <- 3
    ncomp_input <- 2
    test_tensor <- array(1:(n * p * t), dim = c(n, p, t))
    tpca_obj <- tpca(test_tensor, ncomp = ncomp_input)
    expect_equal(length(tpca_obj$explained_variance), ncomp_input)
    expect_equal(dim(tpca_obj$variates), c(n, ncomp_input))
    expect_equal(dim(tpca_obj$loadings), c(p, ncomp_input))
  }
)

test_that(
  "tpca agrees with mixOmics pca",
  code = {
    n <- 3
    p <- 10
    t <- 1
    k <- min(n, p)
    ncomp_input <- 2

    set.seed(1)
    test_x <- array(rnorm(n * p * t, mean = 0, sd = 5), dim = c(n, p, t))

    mixomics_pca <- pca(
      test_x[, , 1],
      ncomp = ncomp_input
    )

    # bltodo: we do not even need to specify identity transforms here right?
    transforms <- matrix_to_m_transforms(diag(1))
    tensor_pca <- tpca(
      test_x,
      ncomp = ncomp_input,
      m = transforms$m,
      minv = transforms$minv
    )

    expect_equal(
      .make_signs_consistent(mixomics_pca$loadings$X),
      .make_signs_consistent(tensor_pca$loadings)
    )
  }
)
