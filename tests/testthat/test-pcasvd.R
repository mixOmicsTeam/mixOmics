context("pcasvd")

# Test 1: Check if the function runs without errors
test_that("pcasvd runs without error", {
  X = matrix(rnorm(100), nrow = 10, ncol = 10)
  res = mixOmics:::pcasvd(X, ncomp = 3)
  expect_s3_class(res, "pca")
})

# Test 2: Check if the number of components is handled correctly
test_that("pcasvd handles ncomp correctly", {
  X = matrix(rnorm(100), nrow = 10, ncol = 10)
  res = mixOmics:::pcasvd(X, ncomp = 2)
  expect_equal(ncol(res$rotation), 2)
})

# Test 3: Check if the scaling and centering options work
test_that("pcasvd scales and centers data when specified", {
  X = matrix(rnorm(100), nrow = 10, ncol = 10)
  res = mixOmics:::pcasvd(X, scale = TRUE, center = TRUE)
  expect_true(!is.null(res$scale))
  expect_true(!is.null(res$center))
})

# Test 4: Check if 'retx = FALSE' works as expected (returns only sdev and rotation)
test_that("pcasvd returns sdev and rotation when retx = FALSE", {
  X = matrix(rnorm(100), nrow = 10, ncol = 10)
  res = mixOmics:::pcasvd(X, retx = FALSE)
  expect_true("sdev" %in% names(res))
  expect_true("rotation" %in% names(res))
  expect_false("X" %in% names(res))
})
