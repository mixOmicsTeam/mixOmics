context('pca')

test_that("pca works with and without reconts as expected", code = {
    data(multidrug)
    pca.res1 <- pca(multidrug$ABC.trans, ncomp = 4, scale = TRUE, reconst = TRUE)
    expect_is(pca.res1, 'pca')
    
    pca.res2 <- pca(multidrug$ABC.trans, ncomp = 4, scale = TRUE, reconst = FALSE)
    expect_is(pca.res2, 'pca')
    
    expect_condition(expect_identical(pca.res1$x, pca.res2$x))
})
