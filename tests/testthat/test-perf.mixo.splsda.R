context("perf.mixo_splsda")
library(BiocParallel)

# Helper function to suppress progress bar output
run_perf_silent <- function(model, validation, folds, nrepeat = 1, BPPARAM, progressBar = TRUE, dist = "all", auc = FALSE) {
  sink(tempfile())  # Suppress progress bar output
  result <- suppressWarnings(
    perf(model, validation = validation, folds = folds, nrepeat = nrepeat, BPPARAM = BPPARAM, progressBar = progressBar, dist = dist, auc = auc)
  )
  sink()  # Restore output
  return(result)
}

## ------------------------------------------------------------------------ ##
## Test perf.mixo_plsda()

test_that("perf.mixo_plsda", code = {
  data(liver.toxicity)
  X <- liver.toxicity$gene
  Y <- liver.toxicity$treatment$Dose.Group
  
  set.seed(12)
  res <- plsda(X, Y, ncomp = 2)
  class(res) # "mint.splsda" "mixo_splsda" "mixo_spls"   "DA"        
  out <- perf(res, validation = "Mfold", folds = 3, nrepeat = 3)
  
  ground.ncomp <- matrix(c(2,1,2,2,1,2), ncol = 3, byrow=T,
                         dimnames = list(c("overall", "BER"),
                                         c("max.dist", "centroids.dist", "mahalanobis.dist")))
  
  expect_equal(out$choice.ncomp, ground.ncomp)
})

test_that("perf.mixo_plsda in serial and parallel", code = {
  data(liver.toxicity)
  X <- liver.toxicity$gene
  Y <- liver.toxicity$treatment$Dose.Group
  res <- plsda(X, Y, ncomp = 2)
  
  # Serial execution
  set.seed(12)
  out <- run_perf_silent(res, validation = "Mfold", folds = 3, nrepeat = 1, BPPARAM = SerialParam(RNGseed = 12))
  
  # Parallel execution
  set.seed(12)
  out.parallel <- run_perf_silent(res, validation = "Mfold", folds = 3, nrepeat = 1, BPPARAM = SnowParam(workers = 2, RNGseed = 12))
  
  # Expect the same result
  expect_equal(out$study.specific.error, out.parallel$study.specific.error)
})

## ------------------------------------------------------------------------ ##
## Test perf.mixo_splsda()

test_that("perf.mixo_splsda", code = {
  data(srbct)
  X <- srbct$gene
  Y <- srbct$class
  
  set.seed(12)
  srbct.splsda <- splsda(X, Y, ncomp = 2, keepX = rep(10, 2))
  class(srbct.splsda)  # Check model class
  
  # Use fewer folds and repeats to speed up the test
  out <- run_perf_silent(srbct.splsda, validation = "Mfold", folds = 5, nrepeat = 1, BPPARAM = SerialParam(RNGseed = 12))
  
  ground.ncomp <- matrix(c(2,2,2,2,2,2), ncol = 3, byrow = TRUE,
                         dimnames = list(c("overall", "BER"),
                                         c("max.dist", "centroids.dist", "mahalanobis.dist")))
  expect_equal(out$choice.ncomp, ground.ncomp)
})

test_that("perf.mixo_splsda in serial and parallel with fewer folds", code = {
  data(srbct)
  X <- srbct$gene
  Y <- srbct$class
  
  set.seed(12)
  srbct.splsda <- splsda(X, Y, ncomp = 2, keepX = rep(10, 2))
  
  # Serial execution
  set.seed(12)
  out <- run_perf_silent(srbct.splsda, validation = "Mfold", folds = 3, nrepeat = 1, BPPARAM = SerialParam(RNGseed = 12))
  
  # Parallel execution
  set.seed(12)
  out.parallel <- run_perf_silent(srbct.splsda, validation = "Mfold", folds = 3, nrepeat = 1, BPPARAM = SnowParam(workers = 2, RNGseed = 12))
  
  # Expect the same result
  expect_equal(out$study.specific.error, out.parallel$study.specific.error)
})
