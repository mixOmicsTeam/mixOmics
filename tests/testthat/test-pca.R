context("pca") ## context helps devtools::test() to categorise tests
## the 'call'slot in all results must be ignored before comparison for identicality
## as it would naturally won't be identical

test_that('pca results for nutrimouse are the same when using either matrix or MAE class',{
  pca.res1 <- pca(X = nutrimouse$lipid)
  pca.res2 <- pca(data = nutrimouse.mae, X = "lipid")
  expect_identical(pca.res1[-1], pca.res2[-1])
})

test_that('assay can be input as a variables',{
  pca.res1 <- pca(X = nutrimouse$lipid)
  pca.assay <- 'lipid'
  pca.res3 <- pca(data = nutrimouse.mae, X = pca.assay)
  expect_identical(pca.res1[-1], pca.res3[-1])
})

test_that('pca returns appropriate error for invalid assay/data',{
  ## expect error
  expect_error(pca(data = nutrimouse.mae, X = "not-an-assay"), class = "inv_assay")
  expect_error(pca(data = nutrimouse, X = "lipid"), class = "inv_data")
  expect_error(pca(X = "lipid"), class = "inv_data")
})

test_that('pca entry checker works',{
  pca.res1 <- pca(X = nutrimouse$lipid)
  ## expect error
  expect_error(pca(data = nutrimouse.mae, X = "not-an-assay"), class = "inv_assay")
  expect_error(pca(data = nutrimouse, X = "lipid"), class = "inv_data")
  expect_error(pca(X = "lipid"), class = "inv_data")
})
