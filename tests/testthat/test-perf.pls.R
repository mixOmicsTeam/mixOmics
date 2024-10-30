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
  # run perf
  pls.perf.obj <- perf(pls.obg, validation = "Mfold", folds = 4, 
                       progressBar = FALSE, nrepeat = 3,
                       BPPARAM = SerialParam(),
                       seed = 12)
  trueVals <- c(0.009, -0.222, -0.332, -0.471)
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
  pls.perf.obj <- perf(pls.obg, validation = "Mfold", folds = 4, 
                       progressBar = FALSE, nrepeat = 3,
                       BPPARAM = SerialParam(),
                       seed = 12)
  trueVals <- c(0.009, -0.222, -0.332, -0.471)
  testVals <- round(pls.perf.obj$measures$Q2.total$summary[, "mean"], 3)
  expect_equal(trueVals, testVals)
  
  # run in parallel with seed
  pls.perf.obj.p <- perf(pls.obg, validation = "Mfold", folds = 4, 
                       progressBar = FALSE, nrepeat = 3,
                       BPPARAM = SnowParam(workers = 2),
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
  trueVals <- list(c(0.032,  0.008,  0.003, -0.006),
                   c(0.022, -0.175, -0.312, -0.437),
                   c(0.088, -0.475, -1.238, -2.218))
  
  for (m in 1:3) {
    suppressWarnings(pls.obg <- pls(Y, X, ncomp = 4, mode = modes[m]))
    suppressWarnings(pls.perf.obj <- perf(pls.obg, validation = "Mfold", folds = 4, 
                         progressBar = F, 
                         nrepeat = 3, BPPARAM = SerialParam(),
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
  model.spls = spls(X, Y, ncomp = 7, mode = 'regression',
                    keepX = c(rep(10, 7)), keepY = c(rep(4, 7)))
  class(model.spls) # "mixo_pls"
  
  # run perf in series
  model.spls.val <- perf(model.spls, validation = "Mfold", folds = 10,
                         BPPARAM = SerialParam(),
                         seed = 12)
  # run perf in parallel
  model.spls.val.p <- perf(model.spls, validation = "Mfold", folds = 10,
                         BPPARAM = SnowParam(workers = 2),
                         seed = 12)
  
  # check values are expected
  trueVals <- c(0.160, -0.114, -0.210, -0.531, -0.707, -0.962, -1.322)
  testVals <- round(model.spls.val$measures$Q2.total$summary[, "mean"], 3)
  expect_equal(trueVals, testVals)
  expect_equal(names(model.spls.val), c("call", "measures", "features"))
  # check the same
  expect_equal(model.spls.val$measures$RSS.upred$summary, model.spls.val.p$measures$RSS.upred$summary)
})

