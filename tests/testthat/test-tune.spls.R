test_that("tune.spls works ", code = {
    data("nutrimouse")
    X <- nutrimouse$gene
    Y <- nutrimouse$lipid
    
    set.seed(42)
    tune.spls.res = tune.spls(X, Y, ncomp = 3,
                              test.keepX = seq(5, 10, 5),
                              test.keepY = seq(3, 6, 3), measure = "cor",
                              folds = 5, nrepeat = 3, progressBar = F)
    
    expect_equal(class(tune.spls.res), c("tune.pls", "tune.spls"))
})

test_that("tune.spls works in parallel", code = {
    
    library(BiocParallel)
    data("nutrimouse")
    X <- nutrimouse$gene
    Y <- nutrimouse$lipid
    
    set.seed(42)
    tune.spls.res = suppressWarnings(tune.spls(X, Y, ncomp = 3,
                                      test.keepX = seq(1, 10, 1),
                                      test.keepY = seq(3, 6, 3), measure = "cor",
                                      folds = 5, nrepeat = 3, progressBar = F,
                                      BPPARAM = SnowParam(workers = 4)))
    
    expect_equal(class(tune.spls.res), c("tune.pls", "tune.spls"))
})