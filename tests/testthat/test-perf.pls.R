context("perf.pls")
library(BiocParallel)

## ------------------------------------------------------------------------ ##
## Test perf.mixo_pls()

test_that("perf() works on pls object in serial and parallel", code = {
  
  # set up data
  data("liver.toxicity")
  X <- liver.toxicity$gene[1:500]
  Y <- liver.toxicity$clinic
  pls.obg <- pls(Y, X, ncomp = 2)
  class(pls.obg) # "mixo_pls"
  
  # run in serial with seed
  pls.perf.obj <- perf(pls.obg, validation = "Mfold", folds = 2, 
                       progressBar = FALSE, nrepeat = 1,
                       BPPARAM = SerialParam(RNGseed = 1000),
                       seed = 12)
  trueVals <- c(0.041, -0.177)
  testVals <- round(pls.perf.obj$measures$Q2.total$summary[, "mean"], 3)
  expect_equal(trueVals, testVals)
  
  # run in parallel with seed - even if set RNGseed as different that is overwritten and gives reproducible results
  pls.perf.obj.p <- perf(pls.obg, validation = "Mfold", folds = 2, 
                       progressBar = FALSE, nrepeat = 1,
                       BPPARAM = SnowParam(workers = 2, RNGseed = 600),
                       seed = 12)
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
  
  # to reproduce error, we need to induce some features to have near zero variance
  X[, c(1, 23, 62, 234, 789)] <- 0
  
  modes <- c("classic", "regression", "canonical")
  trueVals <- list(c(0.021,  -0.014),
                   c(-0.071, -0.260),
                   c(0.089, -0.415))
  
  for (m in 1:3) {
    suppressWarnings(pls.obg <- pls(Y, X, ncomp = 2, mode = modes[m]))
    suppressWarnings(pls.perf.obj <- perf(pls.obg, validation = "Mfold", folds = 2, 
                         progressBar = F, 
                         nrepeat = 1, BPPARAM = SerialParam(),
                         seed = 2))
    
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
  model.spls = spls(X, Y, ncomp = 2, mode = 'regression',
                    keepX = c(10, 10), keepY = c(4, 4))
  class(model.spls) # "mixo_pls"
  
  # run perf in series
  model.spls.val <- perf(model.spls, validation = "Mfold", folds = 2,
                         BPPARAM = SerialParam(),
                         seed = 12)
  # run perf in parallel
  model.spls.val.p <- perf(model.spls, validation = "Mfold", folds = 2,
                         BPPARAM = SnowParam(workers = 2),
                         seed = 12)
  
  # check values are expected
  trueVals <- c(0.096, -0.190)
  testVals <- round(model.spls.val$measures$Q2.total$summary[, "mean"], 3)
  expect_equal(trueVals, testVals)
  expect_equal(names(model.spls.val), c("call", "measures", "features"))
  # check the same
  expect_equal(model.spls.val$measures$RSS.upred$summary, model.spls.val.p$measures$RSS.upred$summary)
})

