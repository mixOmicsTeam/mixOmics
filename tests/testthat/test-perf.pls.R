context("perf.pls")
library(BiocParallel)

## ------------------------------------------------------------------------ ##
## Test perf.mixo_pls()

test_that("perf() works on pls object", code = {
  
  # set up data
  data("liver.toxicity")
  X <- liver.toxicity$gene[1:500]
  Y <- liver.toxicity$clinic
  pls.obg <- pls(Y, X, ncomp = 4)
  class(pls.obg) # "mixo_pls"
  
  set.seed(12)
  pls.perf.obj <- perf(pls.obg, validation = "Mfold", folds = 4, 
                       progressBar = FALSE, nrepeat = 3,
                       BPPARAM = SerialParam())
  trueVals <- c(-0.017, -0.294, -0.431, -0.622)
  testVals <- round(pls.perf.obj$measures$Q2.total$summary[, "mean"], 3)
  expect_equal(trueVals, testVals)
  expect_equal(names(pls.perf.obj), c("call", "measures", "features"))
})

test_that("perf() works on pls object in serial and parallel", code = {
  
  # set up data
  data("liver.toxicity")
  X <- liver.toxicity$gene[1:500]
  Y <- liver.toxicity$clinic
  pls.obg <- pls(Y, X, ncomp = 4)
  class(pls.obg) # "mixo_pls"
  
  # run in serial with seed
  set.seed(12)
  pls.perf.obj <- perf(pls.obg, validation = "Mfold", folds = 4, 
                       progressBar = FALSE, nrepeat = 3,
                       BPPARAM = SerialParam(seed = 12))
  trueVals <- c(-0.017, -0.294, -0.431, -0.622)
  testVals <- round(pls.perf.obj$measures$Q2.total$summary[, "mean"], 3)
  expect_equal(trueVals, testVals)
  
  # run in parallel with seed
  set.seed(12)
  pls.perf.obj.p <- perf(pls.obg, validation = "Mfold", folds = 4, 
                       progressBar = FALSE, nrepeat = 3,
                       BPPARAM = SnowParam(seed = 12, workers = 2))
  trueVals <- c(-0.017, -0.294, -0.431, -0.622)
  testVals <- round(pls.perf.obj.p$measures$Q2.total$summary[, "mean"], 3)
  expect_equal(trueVals, testVals)
  # check the same
  expect_equal(pls.perf.obj$measures$RSS.upred$summary, pls.perf.obj.p$measures$RSS.upred$summary)
})

test_that("perf() works on pls object in serial and parallel - progress bar", code = {
  
  # set up data
  data("liver.toxicity")
  X <- liver.toxicity$gene[1:500]
  Y <- liver.toxicity$clinic
  pls.obg <- pls(Y, X, ncomp = 4)
  class(pls.obg) # "mixo_pls"
  
  # run in serial with seed
  set.seed(12)
  pls.perf.obj <- perf(pls.obg, validation = "Mfold", folds = 4, 
                       progressBar = TRUE, nrepeat = 3,
                       BPPARAM = SerialParam(seed = 12))
  trueVals <- c(-0.017, -0.294, -0.431, -0.622)
  testVals <- round(pls.perf.obj$measures$Q2.total$summary[, "mean"], 3)
  expect_equal(trueVals, testVals)
  
  # run in parallel with seed
  set.seed(12)
  pls.perf.obj.p <- perf(pls.obg, validation = "Mfold", folds = 4, 
                         progressBar = TRUE, nrepeat = 3,
                         BPPARAM = SnowParam(seed = 12, workers = 2))
  trueVals <- c(-0.017, -0.294, -0.431, -0.622)
  testVals <- round(pls.perf.obj.p$measures$Q2.total$summary[, "mean"], 3)
  expect_equal(trueVals, testVals)
  # check the same
  expect_equal(pls.perf.obj$measures$RSS.upred$summary, pls.perf.obj.p$measures$RSS.upred$summary)
})

test_that("perf() works on pls with nzv features (all modes)", code = {
  
  # set up data
  data("liver.toxicity")
  X <- liver.toxicity$gene[1:1000]
  Y <- liver.toxicity$clinic
  class(pls.obg) # "mixo_pls"
  
  # to reproduce error, we need to induce some features to have near zero variance
  X[, c(1, 23, 62, 234, 789)] <- 0
  
  modes <- c("classic", "regression", "canonical")
  trueVals <- list(c(0.031,  0.007,  0.003, -0.006),
                   c(0.006, -0.222, -0.379, -0.553),
                   c(0.089, -0.473, -1.238, -2.228))
  
  for (m in 1:3) {
    set.seed(12)
    suppressWarnings(pls.obg <- pls(Y, X, ncomp = 4, mode = modes[m]))
    suppressWarnings(pls.perf.obj <- perf(pls.obg, validation = "Mfold", folds = 4, 
                         progressBar = F, 
                         nrepeat = 3))
    
    testVals <- round(pls.perf.obj$measures$Q2.total$summary[, "mean"], 3)
    expect_equal(trueVals[[m]], testVals)
  }
  
})

## ------------------------------------------------------------------------ ##
## Test perf.mixo_spls()

test_that("perf() works on spls object in serial and parallel", code = {
  
  # set up data
  data("liver.toxicity")
  X <- liver.toxicity$gene[1:500]
  Y <- liver.toxicity$clinic
  model.spls = spls(X, Y, ncomp = 7, mode = 'regression',
                    keepX = c(rep(10, 7)), keepY = c(rep(4, 7)))
  class(pls.obg) # "mixo_pls"
  
  set.seed(12)
  model.spls.val <- perf(model.spls, validation = "Mfold", folds = 10,
                         BPPARAM = SerialParam(RNGseed = 12))
  set.seed(12)
  model.spls.val.p <- perf(model.spls, validation = "Mfold", folds = 10,
                         BPPARAM = SnowParam(workers = 2, RNGseed = 12))
  
  # check values are expected
  trueVals <- c(0.171, -0.126, -0.251, -0.533, -0.725, -1.050, -1.345)
  testVals <- round(model.spls.val$measures$Q2.total$summary[, "mean"], 3)
  expect_equal(trueVals, testVals)
  expect_equal(names(model.spls.val), c("call", "measures", "features"))
  # check the same
  expect_equal(model.spls.val$measures$RSS.upred$summary, model.spls.val.p$measures$RSS.upred$summary)
})

