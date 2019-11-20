context("test-tune.spls")

test_that("tune.spls works with and without parallel", code = {
    data("nutrimouse")
    X <- nutrimouse$gene
    Y <- nutrimouse$lipid
    test.keepX <- c(5,10)
    ncomp <- 2
    nrepeat <- 3
    
    tune.spls11 <- tune.spls(X, Y, ncomp = ncomp, test.keepX = test.keepX,
                              nrepeat = nrepeat)
    expect_is(tune.spls11, "tune.spls")

    
    tune.spls31 <- tune.spls(X, Y, ncomp = ncomp, test.keepX = test.keepX,
                             nrepeat = nrepeat, cpus = 3)
    expect_is(tune.spls31, "tune.spls")
})
