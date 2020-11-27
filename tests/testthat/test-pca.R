context('pca')

test_that("pca works as expected", code = {
    data(multidrug)
    pca.res1 <- pca(multidrug$ABC.trans, ncomp = 2, scale = TRUE)
    expect_is(pca.res1, 'pca')
})
