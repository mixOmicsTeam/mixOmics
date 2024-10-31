context("perf.mixo_splsda")
library(BiocParallel)

## ------------------------------------------------------------------------ ##
## Test perf.mixo_plsda()

test_that("perf.mixo_plsda in serial and parallel", code = {
  data(liver.toxicity)
  X <- liver.toxicity$gene
  Y <- liver.toxicity$treatment$Dose.Group
  res <- plsda(X, Y, ncomp = 2)
  
  # Serial execution
  out <- suppressWarnings(
    perf(res, validation = "Mfold", folds = 2, nrepeat = 1, BPPARAM = SerialParam(), seed = 12)
  )
  
  # Parallel execution
  out.parallel <- suppressWarnings(
    perf(res, validation = "Mfold", folds = 2, nrepeat = 1, BPPARAM = SnowParam(workers = 2), seed = 12)
  )
  
  # Check result
  ground.ncomp <- matrix(c(2,2,2,2,2,2), ncol = 3, byrow=T,
                         dimnames = list(c("overall", "BER"),
                                         c("max.dist", "centroids.dist", "mahalanobis.dist")))
  expect_equal(out$choice.ncomp, ground.ncomp)
  
  # Expect the same result
  expect_equal(out$error.rate$overall, out.parallel$error.rate$overall)
})

## ------------------------------------------------------------------------ ##
## Test perf.mixo_splsda() (commented out as was slow)
# 
# test_that("perf.mixo_splsda in serial and parallel with fewer folds", code = {
#   
#   # set up data and model
#   data(srbct)
#   X <- srbct$gene
#   Y <- srbct$class
#   srbct.splsda <- splsda(X, Y, ncomp = 2, keepX = rep(10, 2))
#   
#   # Serial execution
#   out <- suppressWarnings(
#     perf(srbct.splsda, validation = "Mfold", folds = 2, nrepeat = 1, 
#               BPPARAM = SerialParam(RNGseed = 1000), seed = 12)
#   )
#   
#   # Parallel execution
#   out.parallel <- suppressWarnings(
#     perf(srbct.splsda, validation = "Mfold", folds = 2, nrepeat = 1, 
#                                   BPPARAM = SnowParam(workers = 2, RNGseed = 1050), seed = 12)
#   )
#   
#   # Check result
#   ground.ncomp <- matrix(c(2,2,2,2,2,2), ncol = 3, byrow=T,
#                          dimnames = list(c("overall", "BER"),
#                                          c("max.dist", "centroids.dist", "mahalanobis.dist")))
#   expect_equal(out$choice.ncomp, ground.ncomp)
#   
#   # Expect the same result
#   expect_equal(out$study.specific.error, out.parallel$study.specific.error)
#   
# })
