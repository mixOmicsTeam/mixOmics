context("test-perf.diablo")

test_that("perf.diablo works with and without parallel processing and with and without auroc", {
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

    perf.res11 = perf(nutrimouse.sgccda, folds = folds, nrepeat = nrep)
    expect_is(perf.res11, "perf.sgccda.mthd")

    perf.res41 = perf(nutrimouse.sgccda, folds = folds, nrepeat = nrep, cpus = 2)
    expect_is(perf.res41, "perf.sgccda.mthd")
    
    perf.res12 = perf.sgccda(nutrimouse.sgccda, folds = folds, nrepeat = nrep, auc = TRUE)
    expect_true("auc" %in% names(perf.res12))
    expect_true("auc.study" %in% names(perf.res12))
    
    perf.res42 = perf.sgccda(nutrimouse.sgccda, folds = folds, nrepeat = nrep, auc = TRUE, cpus = 2)
    expect_true("auc" %in% names(perf.res42))
    expect_true("auc.study" %in% names(perf.res42))

})
