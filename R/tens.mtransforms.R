# ==============================================================================
# m-transforms for the tensor-tensor m-product algebra
# ==============================================================================

#' @examples
#' \dontrun{
#' library(microbenchmark)
#' test <- runif(1000, min = 0, max = 100)
#' microbenchmark(
#'   dttpackage = dtt::dct(test) * 2,
#'   fftapproach = dctii(test)
#' )
#' }

dctii <- function(vec, ortho = TRUE) {
  # perform dct transform using fft algorithm, faster than dtt implementation
  # produces same results as scipy.fft.dct in Python
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

idctii <- function(vec, ortho = TRUE) {
  # TODO: find someone smart who can implement this using stats::fft
  # modify outputs from dtt R library to be consistent with the definitions
  # found in scipt.fft.idct
  # https://scipy.github.io/devdocs/reference/generated/scipy.fftpack.dct.html
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
