context("test-tune.splsda")

test_that("tune.splsda works", {
    data(breast.tumors)
    X = breast.tumors$gene.exp
    Y = as.factor(breast.tumors$sample$treatment)
    RNGversion("3.5.0") ## in case RNG changes!
    RNGkind("L'Ecuyer-CMRG") ## works for parallel as well
    
    set.seed(100)
    tune = tune.splsda(X, Y, ncomp = 2, nrepeat = 3, logratio = "none",
                       test.keepX = c(5, 10, 15), folds = 3, dist = "max.dist")
    
    expect_equal(tune$choice.ncomp$ncomp, 1L)
    expect_equal(tune$choice.keepX, c(comp1 = 10, comp2 = 15))
    

    set.seed(100)
    tune2 = tune.splsda(X, Y, ncomp = 2, nrepeat = 3, logratio = "none",
                       test.keepX = c(5, 10, 15), folds = 3, dist = "max.dist", cpus = 2)
    .almost_identical(tune, tune2)
    
})
