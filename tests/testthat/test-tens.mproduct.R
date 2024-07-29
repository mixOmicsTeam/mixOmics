context("tensor m-transforms")

test_vector <- c(19, 3, 6, 11)
test_tensor1 <- array(1:24, dim = c(2, 4, 3))
test_tesor2 <- array(1:36, dim = c(4, 3, 3))

#' @description
#' ========================
#' Use for internal testing
#' ========================
#' Performs a DCT-II transform using the stats::fft algorithm. Produces
#' identical output to the scipy.fft.dct implementation in Python.
#' @param vec A numeric real vector to transform.
#' @param ortho Logical, if TRUE the output is orthonormal.
#' @return A numeric vector representing the DCT-II of the input.
#' @seealso docs.scipy.org/doc/scipy/tutorial/fft.html
dctii_scipy_equivalent <- function(vec, ortho = TRUE) {
  # https://scipy.github.io/devdocs/reference/generated/scipy.fftpack.dct.html
  # https://stackoverflow.com/questions/11215162
  # NOTE: this is dtt output scaled by 2
  n <- length(vec)
  p <- exp(complex(imaginary = pi / 2 / n) * (seq(2 * n) - 1))
  unscaled_res <- Re(stats::fft(c(vec, rev(vec))) / p)[1:n]

  if (ortho) {
    return(
      c(
        sqrt(1 / (4 * n)) * unscaled_res[1],
        sqrt(1 / (2 * n)) * tail(unscaled_res, n - 1)
      )
    )
  } else {
    return(unscaled_res)
  }
}

#' @description
#' ========================
#' Use for internal testing
#' ========================
naive_facewise_product <- function(a, b) {
  na <- dim(a)[1]
  pa <- dim(a)[2]
  ta <- dim(a)[3]

  nb <- dim(b)[1]
  pb <- dim(b)[2]
  tb <- dim(b)[3]

  stopifnot(ta == tb)
  t <- ta

  stopifnot(pa == nb)
  fp_ab <- array(0, dim = c(na, pb, t))

  for (i in 1:t) {
    fp_ab[, , i] <- a[, , i] %*% b[, , i]
  }

  return(fp_ab)
}

test_that(
  "forward and inverse dctii transforms invert each other",
  code = {
    transforms <- dctii_m_transforms(dim(test_tensor1)[3])
    m <- transforms$m
    minv <- transforms$minv
    expect_equal(test_tensor1, minv(m(test_tensor1)))
  }
)

test_that(
  "facewise product works",
  code = {
    expect_equal(
      naive_facewise_product(test_tensor1, test_tesor2),
      .binary_facewise(test_tensor1, test_tensor2)
    )
  }
)
