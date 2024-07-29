context("tensor m-transforms")

test_vec <- c(19, 3, 6, 11)

#' @description
#' =============================================================================
#' Use for internal testing
#' =============================================================================
#' Performs a DCT-II transform using the stats::fft algorithm. Produces
#' identical output to the scipy.fft.dct implementation in Python.
#' @param vec A numeric real vector to transform.
#' @param ortho Logical, if TRUE the output is orthonormal.
#' @return A numeric vector representing the DCT-II of the input.
#' @seealso docs.scipy.org/doc/scipy/tutorial/fft.html
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

test_that(
  "forward and inverse dctii transforms invert each other",
  code = {
    ortho_configurations <- c(TRUE, FALSE)
    for (ortho_config in ortho_configurations) {
      expect_equal(
        test_vec,
        dctii(idctii(test_vec, ortho = ortho_config), ortho = ortho_config)
      )
    }
  }
)