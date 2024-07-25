context("tensor m-transforms")

test_vec <- c(19, 3, 6, 11)

test_that(
  "forward and inverse dctii transforms invert each other",
  code = {
    ortho_configurations <- c(TRUE, FALSE)
    for (ortho_config in ortho_configurations) {
      expect_equal(
        test_vec,
        dctii(idctii(test_vec, ortho = TRUE), ortho = TRUE)
      )
    }
  }
)

test_that(
  "fast fft-based dctii transform reproducible with dtt library implementation",
  code = {
    expect_equal(
      dtt::dct(test_vec) * 2,
      dctii(test_vec, ortho = FALSE)
    )
  }
)
