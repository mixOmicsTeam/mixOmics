context("tune.spls")

test_that("tune.spls works with and without parallel", code = {
    data("nutrimouse")
    X <- nutrimouse$gene
    Y <- nutrimouse$lipid
    test.keepX <- c(5,10)
    ncomp <- 2
    nrepeat <- 3
    folds <- 3
    
    RNGversion(.mixo_rng()) ## in case RNG changes!
    
    ## ----------- tune.spls works fine ----------- 
    set.seed(42)
    tune.spls11 <-
        tune.spls(
            X = X,
            Y = Y,
            ncomp = ncomp,
            measure = 'cor',
            test.keepX = test.keepX,
            test.keepY = test.keepX,
            folds = folds,
            nrepeat = nrepeat
        )
    expect_equal(c(comp1 = 5, comp2 = 5), tune.spls11$choice.keepX)
    expect_equal(c(comp1 = 5, comp2 = 5), tune.spls11$choice.keepY)
    expect_is(tune.spls11, "tune.spls")
    
    ## ----------- tune(method="spls") same as tune.spls ----------- 
    set.seed(42)
    tune.spls12 <-
        tune(
            method = "spls",
            X = X,
            Y = Y,
            ncomp = ncomp,
            measure = 'cor',
            test.keepX = test.keepX,
            test.keepY = test.keepX,
            folds = folds,
            nrepeat = nrepeat
        )
    
    .almost_identical(tune.spls11, tune.spls12)
    
    ## ----------- tune.spls parallel works fine ----------- 
    # set.seed(12)
    # tune.spls31 <- tune.spls(X, Y, ncomp = ncomp, test.keepX = test.keepX, folds = folds,
    #                          nrepeat = nrepeat, BPPARAM = SnowParam(workers = 2, RNGseed = 12))
                             # TODO tune.spls runs full model nrepeat times for each keepX and keepY!
                             # TODO tune.spls fails with SnowParam
    # expect_is(tune.spls31, "tune.spls")
    # set.seed(12)
    # tune.spls13 <- tune(method = "spls", X = X, Y = Y, ncomp = ncomp, test.keepX = test.keepX, folds = folds,
    #                     nrepeat = nrepeat, cpus = 2)
    # .almost_identical(tune.spls13, tune.spls31)
})
