context("tune.splsda")
library(BiocParallel)

test_that("tune.spls works and is the same in parallel and when run in tune wrapper", code = {
  
  # set up data
  data(breast.tumors)
  X = breast.tumors$gene.exp
  Y = as.factor(breast.tumors$sample$treatment)
  
  # run in serial
  tune.splsda.res.1 = tune.splsda(X, Y, ncomp = 2, nrepeat = 1, logratio = "none",
                                test.keepX = c(5, 10, 15), folds = 2, dist = "max.dist",
                                BPPARAM = SerialParam(), seed = 42)
  
  # run in parallel
  tune.splsda.res.2 = tune.splsda(X, Y, ncomp = 2, nrepeat = 1, logratio = "none",
                                test.keepX = c(5, 10, 15), folds = 2, dist = "max.dist",
                                BPPARAM = SnowParam(workers = 2), seed = 42)
  
  # in tune wrapper in serial
  tune.splsda.res.3 = tune(X, Y, ncomp = 2, nrepeat = 1, logratio = "none",
                                test.keepX = c(5, 10, 15), folds = 2, dist = "max.dist",
                                BPPARAM = SerialParam(), seed = 42,
                         method = "splsda")
  
  # in tune wrapper in parallel
  tune.splsda.res.4 = tune(X, Y, ncomp = 2, nrepeat = 1, logratio = "none",
                         test.keepX = c(5, 10, 15), folds = 2, dist = "max.dist",
                         BPPARAM = SnowParam(workers = 2), seed = 42,
                         method = "splsda")
  
  
  # check outputs
  expect_equal(class(tune.splsda.res.1), "tune.splsda")
  expect_equal(class(tune.splsda.res.2), "tune.splsda")
  expect_equal(class(tune.splsda.res.3), "tune.splsda")
  expect_equal(class(tune.splsda.res.4), "tune.splsda")
  expect_equal(unname(tune.splsda.res.1$choice.keepX), c(10,15))
  expect_equal(unname(tune.splsda.res.2$choice.keepX), c(10,15))
  expect_equal(unname(tune.splsda.res.3$choice.keepX), c(10,15))
  expect_equal(unname(tune.splsda.res.4$choice.keepX), c(10,15))
  .expect_numerically_close(tune.splsda.res.1$error.rate[1,1], 0.3111111)
  .expect_numerically_close(tune.splsda.res.2$error.rate[1,1], 0.3111111)
  .expect_numerically_close(tune.splsda.res.3$error.rate[1,1], 0.3111111)
  .expect_numerically_close(tune.splsda.res.4$error.rate[1,1], 0.3111111)
  
})

test_that("tune.splsda works when test.keepX = NULL and gives same result as perf()", code = {
  
  # set up data
  data(breast.tumors)
  X <- breast.tumors$gene.exp
  Y <- as.factor(breast.tumors$sample$treatment)
  
  # tune on components only
  tune_res <- suppressWarnings(
    tune.splsda(X, Y, ncomp = 2, logratio = "none",
                     nrepeat = 1, folds = 3,
                     test.keepX = NULL, seed = 20)
  )
  
  # run perf
  splsda_res <- splsda(X, Y, ncomp = 2)
  perf_res <- suppressWarnings(
    perf(splsda_res, ncomp = 2, nrepeat = 1, folds = 3, seed = 20)
  )
  
  # check results match
  expect_equal(tune_res$error.rate$overall[1,2], perf_res$error.rate$overall[1,2])
})
