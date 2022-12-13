test_that("tune.spls works ", code = {
    data("nutrimouse")
    X <- nutrimouse$gene
    Y <- nutrimouse$lipid
    
    set.seed(42)
    tune.spls.res = suppressWarnings(tune.spls(X, Y, ncomp = 3,
                              test.keepX = seq(5, 10, 5),
                              test.keepY = seq(3, 6, 3), measure = "cor",
                              folds = 5, nrepeat = 3, progressBar = F,
                              BPPARAM = SerialParam(RNGseed = 5212)))
    
    expect_equal(class(tune.spls.res), c("tune.pls", "tune.spls"))
    expect_equal(unname(tune.spls.res$choice.keepX), c(5,5,5))
    expect_equal(unname(tune.spls.res$choice.keepY), c(3,3,6))
})

test_that("tune.spls works in parallel", code = {
    
    library(BiocParallel)
    data("nutrimouse")
    X <- nutrimouse$gene
    Y <- nutrimouse$lipid
    
    set.seed(42)
    tune.spls.res = suppressWarnings(tune.spls(X, Y, ncomp = 3,
                                      test.keepX = seq(1, 5, 1),
                                      test.keepY = seq(3, 6, 3), measure = "cor",
                                      folds = 5, nrepeat = 3, progressBar = F,
                                      BPPARAM = SnowParam(RNGseed = 5212, workers = 4)))
    
    expect_equal(class(tune.spls.res), c("tune.pls", "tune.spls"))
    expect_equal(unname(tune.spls.res$choice.keepX), c(1,1,1))
    expect_equal(unname(tune.spls.res$choice.keepY), c(3,3,3))
})

test_that("tune.spls and tune(method='spls') are equivalent", {
    
    data("nutrimouse")
    X <- nutrimouse$gene
    Y <- nutrimouse$lipid
    
    set.seed(42)
    tune.spls.res.1 = suppressWarnings(tune.spls(X, Y, ncomp = 3,
                                      test.keepX = seq(1, 2, 1),
                                      test.keepY = seq(3, 6, 3),
                                      folds = 5, nrepeat = 3, progressBar = F,
                                      BPPARAM = SerialParam(RNGseed = 5212)))
    
    tune.spls.res.2 = suppressWarnings(tune(method = "spls", X, Y, ncomp = 3,
                                      test.keepX = seq(1, 2, 1),
                                      test.keepY = seq(3, 6, 3),
                                      folds = 5, nrepeat = 3, progressBar = F,
                                      BPPARAM = SerialParam(RNGseed = 5212)))
    
    expect_equal(tune.spls.res.1$measure.pred, tune.spls.res.2$measure.pred)
})
