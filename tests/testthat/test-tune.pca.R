context("tune.pca")

test_that("tune.pca and tune(method='pca') are equivalent", {
  data(srbct)
  X <- srbct$gene[1:20, 1:200]
  object1 <- tune.pca(X, ncomp = 2, center = TRUE, scale = TRUE)
  object2 <- tune(method = "pca", X, ncomp = 2, center = TRUE, scale = TRUE)
  # expect results the same
  expect_equal(object1$ncomp, object2$ncomp)
  expect_equal(object1$sdev[1], object2$sdev[1])
  # check can plot outputs without errors
  expect_silent(plot(object1))
  expect_silent(plot(object2))
})
