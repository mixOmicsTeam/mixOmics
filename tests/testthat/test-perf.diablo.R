context("test-perf.diablo")

test_that("perf.diablo works with and without parallel processing", {
    data(nutrimouse)
    nrep <- 6
    seed <- 100
    Y = nutrimouse$diet
    data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
    design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
    
    nutrimouse.sgccda <- block.splsda(X=data,
                                      Y = Y,
                                      design = design,
                                      keepX = list(gene=c(10,10), lipid=c(15,15)),
                                      ncomp = 2,
                                      scheme = "horst")
    set.seed(42)
    perf.res11 = perf(nutrimouse.sgccda, folds = 2, nrepeat = nrep)
    expect_is(perf.res11, "perf.sgccda.mthd")
    set.seed(42)
    perf.res41 = perf(nutrimouse.sgccda, folds = 2, nrepeat = nrep, cpus=2)
    expect_is(perf.res41, "perf.sgccda.mthd")
})
