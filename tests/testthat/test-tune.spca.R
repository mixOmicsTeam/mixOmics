context("tune.spca")
library(BiocParallel)

test_that("tune.spca works in serial and parallel", {
  data(srbct)
  X <- srbct$gene[1:20, 1:200]
  grid.keepX <- seq(5, 35, 10)
  object_serial <- tune.spca(X,ncomp = 2, 
                      folds = 5, 
                      test.keepX = grid.keepX, nrepeat = 3,
                      BPPARAM = SerialParam(RNGseed = NULL), seed = 5212)
  expect_equal(object_serial$choice.keepX[[1]], 35)
  expect_equal(object_serial$choice.keepX[[2]], 5)
  .expect_numerically_close(object_serial$cor.comp$comp1[1,2], 0.3994544)
  object_parallel <- tune.spca(X,ncomp = 2, 
                      folds = 5, 
                      test.keepX = grid.keepX, nrepeat = 3,
                      BPPARAM = MulticoreParam(), seed = 5212)
  expect_equal(object_parallel$choice.keepX[[1]], 35)
  expect_equal(object_parallel$choice.keepX[[2]], 5)
  .expect_numerically_close(object_parallel$cor.comp$comp1[1,2], 0.3994544)
  # test can plot outputs
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot(object_serial))
  expect_silent(plot(object_parallel))
})

test_that("tune.spca same result in serial and parallel", {
  data(srbct)
  X <- srbct$gene[1:20, 1:200]
  grid.keepX <- seq(5, 35, 10)
  serial_time <- system.time(
    object_serial <- tune.spca(X,ncomp = 2, 
                               folds = 5, 
                               test.keepX = grid.keepX, nrepeat = 20,
                               BPPARAM = SerialParam(RNGseed = 10000), seed = 5212) # RNGseed in BPPARAM ignored
  )
  parallel_time <- system.time(
    object_parallel <- tune.spca(X,ncomp = 2, 
                               folds = 5, 
                               test.keepX = grid.keepX, nrepeat = 20,
                               BPPARAM = MulticoreParam(), seed = 5212)
  )
  # expect parallel faster - doesn't actually work in this case so commented out!
  # expect_true(serial_time[3] > parallel_time[3])
  # expect results the same
  expect_equal(object_parallel$choice.keepX[[1]], object_serial$choice.keepX[[1]])
  expect_equal(object_parallel$choice.keepX[[2]], object_serial$choice.keepX[[2]])
  .expect_numerically_close(object_parallel$cor.comp$comp1[3,3], object_serial$cor.comp$comp1[3,3])
})

test_that("tune.spca works with NA input", {
  data(srbct)
  X <- srbct$gene[1:20, 1:200]
  set.seed(5212) # set here although this actually doesnt affect tune.spca
  na.feats <- sample(1:ncol(X), 20)
  for (c in na.feats) {
    na.samples <- sample(1:nrow(X),1) #  sample.int(3, 1)
    X[na.samples, c] <- NA
  }
  grid.keepX <- seq(5, 35, 10)
  
  expect_warning({object <- tune.spca(X,ncomp = 2, 
                                      folds = 5, 
                                      test.keepX = grid.keepX, nrepeat = 3,
                                      BPPARAM = SerialParam(), seed = 5212)},
                 "NAs present")
  expect_equal(object$choice.keepX[[1]], 15)
  expect_equal(object$choice.keepX[[2]], 5)
})

test_that("tune.spca and tune(method='spca') are equivalent", {
  data(srbct)
  X <- srbct$gene[1:20, 1:200]
  grid.keepX <- seq(5, 35, 10)
  object1 <- tune.spca(X, ncomp = 2, folds = 5, test.keepX = grid.keepX, nrepeat = 3,
                             BPPARAM = SerialParam(), seed = 5212)
  object2 <- tune(method = "spca",X, ncomp = 2, folds = 5, test.keepX = grid.keepX, nrepeat = 3,
                  BPPARAM = SerialParam(), seed = 5212)
  # expect results the same
  expect_equal(object1$choice.keepX[[1]], object2$choice.keepX[[1]])
  expect_equal(object1$choice.keepX[[2]], object2$choice.keepX[[2]])
  .expect_numerically_close(object1$cor.comp$comp1[3,3], object2$cor.comp$comp1[3,3])
  # test can plot outputs
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot(object1))
  expect_silent(plot(object2))
})
  