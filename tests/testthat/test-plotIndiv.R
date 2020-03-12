context("plotIndiv")

test_that("plotIndiv.rcc works without ind.names", code = {
    data(nutrimouse)
    X <- nutrimouse$lipid
    Y <- nutrimouse$gene
    nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
    
    plotIndiv.res <- tryCatch(expr = {plotIndiv(nutri.res, group = nutrimouse$genotype, ind.names = FALSE, legend = TRUE)},
                              error = function(e) e,
                              warning = function(w) w
    )
    
    expect_is(plotIndiv.res$graph, "ggplot")
})
