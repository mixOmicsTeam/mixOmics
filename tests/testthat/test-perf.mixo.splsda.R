context("perf.mixo_splsda")
library(BiocParallel)

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
  # in serial
  set.seed(12)
  out <- perf(res, validation = "Mfold", folds = 3, nrepeat = 3, BPPARAM = SerialParam(RNGseed = 12))
  # in parallel
  set.seed(12)
  out.parallel <- perf(res, validation = "Mfold", folds = 3, nrepeat = 3, BPPARAM = SnowParam(workers = 2, RNGseed = 12))
  # expect same result
  expect_equal(out$study.specific.error, out.parallel$study.specific.error)
  
})

test_that("perf.mixo_plsda in serial and parallel with progress bar", code = {
  data(liver.toxicity)
  X <- liver.toxicity$gene
  Y <- liver.toxicity$treatment$Dose.Group
  res <- plsda(X, Y, ncomp = 2)
  # in serial
  set.seed(44)
  out <- perf(res, validation = "Mfold", folds = 3, nrepeat = 10, 
              BPPARAM = SerialParam(RNGseed = 12), progressBar = TRUE)
  # in parallel
  set.seed(44)
  out.parallel <- perf(res, validation = "Mfold", folds = 3, nrepeat = 10, 
                       BPPARAM = SnowParam(workers = 2, RNGseed = 12), progressBar = TRUE)
  # expect same result
  expect_equal(out$study.specific.error, out.parallel$study.specific.error)
  
})

test_that("does not allow for class with 1 associated sample", code = {
    data(liver.toxicity)
    X <- liver.toxicity$gene
    Y <- liver.toxicity$treatment$Dose.Group
    # create a class with one sample only
    Y[c(1)] <- 'random.class'

    res <- plsda(X, Y, ncomp = 2)

    expect_error(perf(res, validation = "Mfold", folds = 3, nrepeat = 3),
                 "single assocaited")
})

## ------------------------------------------------------------------------ ##
## Test perf.mixo_splsda()

test_that("perf.mixo_splsda", code = {
  data(srbct)
  X <- srbct$gene
  Y <- srbct$class
  
  set.seed(12)
  srbct.splsda <- splsda(X, Y, ncomp = 2, keepX = rep(10, 2))
  class(res) #  "mixo_splsda" "mixo_spls"   "DA"   
  set.seed(45)
  out <- suppressWarnings(
    perf(srbct.splsda, validation = "Mfold", folds = 8,
                dist = "all", auc = TRUE)
  )
  
  ground.ncomp <- matrix(c(2,2,2,2,2,2), ncol = 3, byrow=T,
                         dimnames = list(c("overall", "BER"),
                                         c("max.dist", "centroids.dist", "mahalanobis.dist")))
  expect_equal(out$choice.ncomp, ground.ncomp)
})

test_that("perf.mixo_splsda in serial and parallel", code = {
  data(srbct)
  X <- srbct$gene
  Y <- srbct$class
  
  set.seed(12)
  srbct.splsda <- splsda(X, Y, ncomp = 2, keepX = rep(10, 2))
  class(res) #  "mixo_splsda" "mixo_spls"   "DA"   

  # in serial
  set.seed(12)
  out <- suppressWarnings(
    perf(srbct.splsda, validation = "Mfold", folds = 12, dist = "all", auc = FALSE,
         BPPARAM = SerialParam(RNGseed = 12))
  )
  # in parallel
  set.seed(12)
  out.parallel <- suppressWarnings(
    perf(srbct.splsda, validation = "Mfold", folds = 12, dist = "all", auc = FALSE,
         BPPARAM = SnowParam(workers = 2, RNGseed = 12))
  )
  # expect same result
  expect_equal(out$study.specific.error, out.parallel$study.specific.error)
  
})

test_that("perf.mixo_splsda in serial and parallel with progress bars", code = {
  data(srbct)
  X <- srbct$gene
  Y <- srbct$class
  
  set.seed(12)
  srbct.splsda <- splsda(X, Y, ncomp = 2, keepX = rep(10, 2))
  class(res) #  "mixo_splsda" "mixo_spls"   "DA"   
  
  # in serial
  set.seed(12)
  out <- suppressWarnings(
    perf(srbct.splsda, validation = "Mfold", folds = 5, dist = "all", auc = FALSE,
         BPPARAM = SerialParam(RNGseed = 12), progressBar = TRUE)
  )
  # in parallel
  set.seed(12)
  out.parallel <- suppressWarnings(
    perf(srbct.splsda, validation = "Mfold", folds = 5, dist = "all", auc = FALSE,
         BPPARAM = SnowParam(workers = 2, RNGseed = 12), progressBar = TRUE)
  )
  # expect same result
  expect_equal(out$study.specific.error, out.parallel$study.specific.error)
  
})
