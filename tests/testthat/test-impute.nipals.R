context("impute.nipals")

test_that("impute.nipals returns X unchanged when there are no missing values", {
  X <- matrix(runif(20), nrow = 5, ncol = 4)
  expect_output(imputed_X <- impute.nipals(X, ncomp = 2), "no missing values in 'X' to impute")
  expect_equal(imputed_X, X)
})

test_that("impute.nipals warns when ncomp is too low", {
  X <- matrix(runif(20), nrow = 5, ncol = 4)
  X[2, 3] <- NA
  expect_message(impute.nipals(X, ncomp = 2), "consider high 'ncomp' for more accurate imputation")
})

test_that("impute.nipals performs imputation correctly", {
  X <- matrix(runif(20), nrow = 5, ncol = 4)
  X[2, 3] <- NA
  ncomp_valid <- min(dim(X))  # Ensure ncomp is valid (in this case, ncomp = 4)
  
  imputed_X <- impute.nipals(X, ncomp = ncomp_valid)  # Use valid ncomp
  
  expect_false(any(is.na(imputed_X)))  # No missing values after imputation
  expect_equal(dim(imputed_X), dim(X))  # Dimensions should remain the same
})