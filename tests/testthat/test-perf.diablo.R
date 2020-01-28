context("test-perf.diablo")

test_that("perf.diablo works with and without parallel processing and with auroc", {
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
    
    RNGversion(.mixo_rng()) ## in case RNG changes!
    set.seed(100)
    perf.res12 = perf.sgccda(nutrimouse.sgccda, folds = folds, nrepeat = nrep, auc = TRUE, progressBar = TRUE)
    choices <- unname(perf.res12$choice.ncomp$AveragedPredict[,1])
    expect_equal(choices, c(2,2))
    aucs <- round(unname(perf.res12$auc$comp1[,1]), 2)
    expect_equal(aucs, c(0.92, 0.74, 0.66, 0.55, 0.69))

    # with cpus seed must be designated after cluster is created so I made it possible
    # by listening to ... in perf for seed. Results are different even with seeds but reproducible with same cpus
    # the hassle of making it fully reproducible is a bit too arduous
    
    perf.res42 = perf.sgccda(nutrimouse.sgccda, folds = folds, nrepeat = nrep, auc = TRUE, cpus = 2, seed = 100, progressBar = TRUE)
    choices <- unname(perf.res42$choice.ncomp$AveragedPredict[,1])
    expect_equal(choices, c(1,1))
    aucs <- round(unname(perf.res42$auc$comp1[,1]), 2)
    expect_equal(aucs,c(0.97, 0.66, 0.6, 0.59, 0.79))
})
