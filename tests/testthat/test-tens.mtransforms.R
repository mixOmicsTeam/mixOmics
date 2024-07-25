context("tensor m-transforms")

test_vec <- c(19, 3, 6, 11)

test_that(
  "forward and inverse dctii transforms invert each other",
  code = {
    ortho_configurations <- c(TRUE, FALSE)
    for (ortho_config in ortho_configurations) {
      expect_equal(test_vec, dctii(idct(test_vec, ortho = TRUE), ortho = TRUE))
    }
  }
)
