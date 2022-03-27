context("perf.pls")

test_that("perf() works on pls object", code = {
  library(mixOmics)
  
  data("liver.toxicity")
  
  # reducing number of features to reduce run time
  X <- liver.toxicity$gene[1:500]
  Y <- liver.toxicity$clinic
  
  set.seed(12)
  pls.obg <- pls(Y, X, ncomp = 4)
  pls.perf.obj <- perf(pls.obg, validation = "Mfold", folds = 4, 
                       progressBar = F, 
                       nrepeat = 3)
  
  trueVals <- c(-0.017, -0.294, -0.431, -0.622)
  testVals <- round(pls.perf.obj$measures$Q2.total$summary[, "mean"], 3)
  
  expect_equal(trueVals, testVals)
})

test_that("perf() works on pls with nzv features (all modes)", code = {
  library(mixOmics)
  
  data("liver.toxicity")
  
  # reducing number of features to reduce run time
  X <- liver.toxicity$gene[, 1:1000]
  Y <- liver.toxicity$clinic
  
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