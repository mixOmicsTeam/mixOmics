context("tune.plsda")
library(BiocParallel)

test_that("tune.plsda works and is the same perf alone and in tune wrapper", code = {
  
  # set up data
  data(breast.tumors)
  X = breast.tumors$gene.exp
  Y = as.factor(breast.tumors$sample$treatment)
  
  # run alone
  tune.plsda.res.1 = suppressWarnings(
    tune.splsda(X, Y, ncomp = 2, nrepeat = 1, logratio = "none",
                                  folds = 2,
                                  BPPARAM = SerialParam(), seed = 42)
  )
  
  # run in wrapper
  tune.plsda.res.2 = suppressWarnings(
    tune(X, Y, ncomp = 2, nrepeat = 1, logratio = "none",
                                  folds = 2,
                                  BPPARAM = SnowParam(workers = 2), seed = 42,
                          method = "plsda")
  )
  
  # run perf
  model <- plsda(X, Y, ncomp = 2, logratio = "none")
  tune.plsda.res.3 = suppressWarnings(
    perf(model, ncomp = 2, nrepeat = 1,
                           folds = 2,
                           BPPARAM = SerialParam(), seed = 42,
                           method = "splsda")
  )
  
  
  # check outputs format
  expect_equal(class(tune.plsda.res.1)[1], "perf")
  expect_equal(class(tune.plsda.res.2)[1], "perf")
  expect_equal(class(tune.plsda.res.3)[1], "perf")
  # check outputs values
  .expect_numerically_close(tune.plsda.res.1$error.rate$overall[1,1], 0.5106383)
  .expect_numerically_close(tune.plsda.res.2$error.rate$overall[1,1], 0.5106383)
  .expect_numerically_close(tune.plsda.res.3$error.rate$overall[1,1], 0.5106383)
  
})