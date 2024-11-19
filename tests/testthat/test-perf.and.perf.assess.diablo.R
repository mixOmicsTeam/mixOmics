context("perf.diablo")
library(BiocParallel)

## ------------------------------------------------------------------------ ##
## Test perf.sgccda()

test_that("perf.diablo works with with auroc", {
  # set up data and model
    data(nutrimouse)
    data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
    design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
    nutrimouse.sgccda <- block.splsda(X = data,
                                      Y = nutrimouse$diet,
                                      design = design,
                                      keepX = list(gene=c(10,10), lipid=c(15,15)),
                                      ncomp = 2)
    
    # run perf model - set seed as 100, ignores RNGseed
    perf.res <- perf(nutrimouse.sgccda, folds = 2, nrepeat = 1, auc = TRUE, 
                      BPPARAM = SerialParam(RNGseed = 100000), seed = 100, progressBar = FALSE)
    true_aucs <- c(0.95, 0.62, 0.68, 0.54, 0.77)
    aucs <- round(unname(perf.res$auc$comp1[,1]), 2)
    expect_equal(aucs, true_aucs)
    
    # run in parallel
    perf.res.parallel = perf(nutrimouse.sgccda, folds = 2, nrepeat = 1, auc = TRUE, 
                      BPPARAM = SnowParam(workers = 2), seed = 100, progressBar = FALSE)
    aucs <- round(unname(perf.res.parallel$auc$comp1[,1]), 2)
    expect_equal(aucs, true_aucs)
    expect_equal(perf.res$weights, perf.res.parallel$weights)
    
    # run perf.assess
    perf.assess.res <- perf.assess(nutrimouse.sgccda, folds = 2, nrepeat = 1, auc = TRUE, 
                            BPPARAM = SerialParam(RNGseed = 100000), seed = 100, progressBar = FALSE)
    expect_equal(perf.res$weights$comp2, perf.assess.res$weights$comp2)
    
    # run perf.assess in parallel
    perf.assess.res.parallel <- perf.assess(nutrimouse.sgccda, folds = 2, nrepeat = 1, auc = TRUE, 
                                   BPPARAM = SnowParam(workers = 2), seed = 100, progressBar = FALSE)
    expect_equal(perf.res$weights$comp2, perf.assess.res.parallel$weights$comp2)
    
})
