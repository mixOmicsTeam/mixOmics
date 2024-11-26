context("tune.rcc")
library(BiocParallel)

test_that("tune.rcc works with Mfold method", code = {
  # set up data
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  
  # run
  tune.rcc.res <- tune.rcc(X, Y, validation = "Mfold", seed = 20)
  
  # check outputs
  expect_equal(class(tune.rcc.res), "tune.rcc")
  expect_equal(tune.rcc.res$opt.lambda1, 0.5005)
  expect_equal(tune.rcc.res$grid1, c(0.00100, 0.25075, 0.50050, 0.75025, 1.00000))
})

test_that("tune.rcc works with loo method", code = {
  # set up data
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  
  # run
  tune.rcc.res <- tune.rcc(X, Y, validation = "loo", seed = 20)
  
  # check outputs
  expect_equal(class(tune.rcc.res), "tune.rcc")
  expect_equal(tune.rcc.res$opt.lambda1, 0.25075)
  expect_equal(tune.rcc.res$grid1, c(0.00100, 0.25075, 0.50050, 0.75025, 1.00000))
})
  
test_that("tune.rcc works in parallel same as in series", code = {
  # set up data
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  
  # run in series
  tune.rcc.res <- tune.rcc(X, Y, validation = "Mfold",
                           BPPARAM = SerialParam(RNGseed = NULL), seed = 12)
  # run in parallel
  tune.rcc.res.parallel <- tune.rcc(X, Y, validation = "Mfold",
                           BPPARAM = SnowParam(workers = 2, RNGseed = NULL), seed = 12)
  
  # check outputs
  expect_equal(class(tune.rcc.res), "tune.rcc")
  expect_equal(tune.rcc.res$opt.lambda2, tune.rcc.res.parallel$opt.lambda2)
  expect_equal(tune.rcc.res$grid2, tune.rcc.res.parallel$grid2)
})

test_that("tune.rcc and tune(method='rcc') are equivalent", {
  
  # set up data
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  
  # run independently
  tune.rcc.res.1 <- tune.rcc(X, Y, validation = "Mfold",
                           BPPARAM = SerialParam(RNGseed = NULL), seed = 12,
                           grid1 = c(0.001, 0.2, 0.6, 1),
                           grid2 = c(0.001, 0.2, 0.6, 1))
  
  # run in tune wrapper
  tune.rcc.res.2 <- tune(X, Y, validation = "Mfold", 
                           BPPARAM = SerialParam(), seed = 12,
                           grid1 = c(0.001, 0.2, 0.6, 1),
                           grid2 = c(0.001, 0.2, 0.6, 1),
                         method = "rcc")
  
  # check outputs
  expect_equal(class(tune.rcc.res.1), "tune.rcc")
  expect_equal(tune.rcc.res.1$opt.lambda2, tune.rcc.res.2$opt.lambda2)
  expect_equal(tune.rcc.res.1$grid2, tune.rcc.res.2$grid2)
  
})

