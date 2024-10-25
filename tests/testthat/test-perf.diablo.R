context("perf.diablo")
library(BiocParallel)

## ------------------------------------------------------------------------ ##
## Test perf.sgccda()

test_that("perf.diablo works ", {
    data(nutrimouse)
    Y = nutrimouse$diet
    data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
    design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
    nutrimouse.sgccda <- block.splsda(X=data, Y = Y,design = design, keepX = list(gene=c(10,10), lipid=c(15,15)), ncomp = 2, scheme = "horst")
    perf = perf(nutrimouse.sgccda, folds = 3, nrepeat = 2)
    expect_is(perf, "perf.sgccda.mthd")
})

test_that("perf.diablo works with and without perf wrapper and with auroc", {
    data(nutrimouse)
    nrep <- 3
    folds <- 2
    Y = nutrimouse$diet
    data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
    design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
    
    nutrimouse.sgccda <- block.splsda(X=data,
                                      Y = Y,
                                      design = design,
                                      keepX = list(gene=c(10,10), lipid=c(15,15)),
                                      ncomp = 2,
                                      scheme = "horst")
    # run perf.sgccda() directly
    set.seed(100)
    perf.res12 = perf.sgccda(nutrimouse.sgccda, folds = folds, nrepeat = nrep, auc = TRUE, 
                             progressBar = FALSE, BPPARAM = SerialParam(RNGseed = 100))
    choices <- unname(perf.res12$choice.ncomp$AveragedPredict[,1])
    expect_equal(choices, c(1,1))
    aucs <- round(unname(perf.res12$auc$comp1[,1]), 2)
    expect_equal(aucs, c(0.96, 0.76, 0.64, 0.53, 0.70))
    # run perf.sgccda() through perf()
    perf.res42 = perf(nutrimouse.sgccda, folds = folds, nrepeat = nrep, auc = TRUE, 
                      BPPARAM = SerialParam(RNGseed = 100), seed = 100, progressBar = FALSE)
    choices <- unname(perf.res42$choice.ncomp$AveragedPredict[,1])
    expect_equal(choices, c(1,1))
    aucs <- round(unname(perf.res42$auc$comp1[,1]), 2)
    expect_equal(aucs, c(0.96, 0.76, 0.64, 0.53, 0.70))
})

test_that("perf.diablo works with and without parallelisation", {
  data(nutrimouse)
  Y = nutrimouse$diet
  data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
  
  nutrimouse.sgccda <- block.splsda(X=data,
                                    Y = Y,
                                    design = 'full',
                                    keepX = list(gene=c(5,10), lipid=c(10,15)),
                                    ncomp = 2,
                                    scheme = "horst")
  class(nutrimouse.sgccda) # "block.splsda" "block.spls"   "sgccda"       "sgcca"        "DA"
  
  # serial
  set.seed(50)
  perf <- perf(nutrimouse.sgccda, folds = 6, BPPARAM = SerialParam(RNGseed = 50))
  # parallel
  set.seed(50)
  perf.parallel <- perf(nutrimouse.sgccda, folds = 6, BPPARAM = SnowParam(workers = 2, RNGseed = 50))
  # check result same
  expect_equal(perf$weights, perf.parallel$weights)
})

test_that("perf.diablo works with and without parallelisation with progress bar", {
  data(nutrimouse)
  Y = nutrimouse$diet
  data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)

  nutrimouse.sgccda <- block.splsda(X=data,
                                    Y = Y,
                                    design = 'full',
                                    keepX = list(gene=c(5,10), lipid=c(10,15)),
                                    ncomp = 2,
                                    scheme = "horst")
  class(nutrimouse.sgccda) # "block.splsda" "block.spls"   "sgccda"       "sgcca"        "DA"

  # serial
  set.seed(100)
  sink(tempfile()) # added to hide progress bars during testing
  perf <- perf(nutrimouse.sgccda, folds = 6, BPPARAM = SerialParam(RNGseed = 100), progressBar = TRUE)
  sink()
  # parallel
  set.seed(100)
  sink(tempfile()) # added to hide progress bars during testing
  perf.parallel <- perf(nutrimouse.sgccda, folds = 6, BPPARAM = SnowParam(workers = 2, RNGseed = 100), progressBar = TRUE)
  sink()
  # check result same
  expect_equal(perf$weights, perf.parallel$weights)
})
