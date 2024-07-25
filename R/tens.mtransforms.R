# ==============================================================================
# m-transforms for the tensor-tensor m-product algebra
# ==============================================================================

#' @description
#' Performs a DCT-II transform using the stats::fft algorithm. Produces 
#' identical output to the scipy.fft.dct implementation in Python.
#' @param vec A numeric real vector to transform.
#' @param ortho Logical, if TRUE the output is orthonormal.
#' @return A numeric vector representing the DCT-II of the input.
#' @author Brendan Lu
#' @seealso docs.scipy.org/doc/scipy/tutorial/fft.html
#' @export
dctii <- function(vec, ortho = TRUE) {
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

# library(microbenchmark)
# test_vec <- runif(1000, min = 0, max = 100)
# microbenchmark(
#   dttpackage = dtt::dct(test_vec) * 2,
#   fftapproach = dctii(test_vec, ortho = FALSE)
# )

#' @description
#' Performs inverse DCT-II transform using the DDT package. Produces identical
#' output to the scipy.fft.dct implementation in Python.
#' @param vec A numeric real vector to transform.
#' @param ortho Logical, if TRUE the output is orthonormal.
#' @return A numeric vector representing the DCT-II of the input.
#' @author Brendan Lu
#' @seealso docs.scipy.org/doc/scipy/tutorial/fft.html
#' @export
idctii <- function(vec, ortho = TRUE) {
  # TODO: find someone smart who can implement this using stats::fft
  # https://dsp.stackexchange.com/questions/51311
  # modify outputs from dtt R library to be consistent with the definitions
  # found in scipt.fft.idct
  n <- length(vec)
  unscaled_res <- dtt::dct(vec, inverted = TRUE) / 2

  if (ortho) {
    return(
      (unscaled_res * 2 * n - vec[1]) * (1 / sqrt(2 * n)) + vec[1] / sqrt(n)
    )
  } else {
    return(unscaled_res)
  }
}
