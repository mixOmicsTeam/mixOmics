context("test-tune.spls")

test_that("tune.spls works with and without parallel", code = {
    data("nutrimouse")
    X <- nutrimouse$gene
    Y <- nutrimouse$lipid
    test.keepX <- c(5,10)
    ncomp <- 2
    nrepeat <- 3
    
    RNGversion("3.5.0") ## in case RNG changes!
    RNGkind("L'Ecuyer-CMRG") ## works for parallel as well
    
    set.seed(100)
    tune.spls11 <- tune.spls(X, Y, ncomp = ncomp, test.keepX = test.keepX,
                              nrepeat = nrepeat)
    expect_is(tune.spls11, "tune.spls")
    
    set.seed(100)
    tune.spls31 <- tune.spls(X, Y, ncomp = ncomp, test.keepX = test.keepX,
                             nrepeat = nrepeat, cpus = 2)
    expect_identical2(tune.spls11, tune.spls31)
})
