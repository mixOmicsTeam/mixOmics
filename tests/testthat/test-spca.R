context("spca") ## context helps devtools::test() to categorise tests
## the 'call'slot in all results must be ignored before comparison for identicality
## as it would naturally won't be identical

## ------- parameters used are defined in `testthat/helper-*.R` files
test_that('spca results for nutrimouse are the same when using either matrix or MAE class',{
  spca.res1 <- spca(X = nutrimouse$lipid, keepX = c(5,5))
  spca.res2 <- spca(data = nutrimouse.mae, X = "lipid", keepX = c(5,5))
  expect_identical(spca.res1[-1], spca.res2[-1])
})

test_that('assay name can be input as a variable',{
  spca.res1 <- spca(X = nutrimouse$lipid, keepX = c(5,5))
  spca.res3 <- spca(data = nutrimouse.mae, X = pca.assay, keepX = c(5,5))
  expect_identical(spca.res1[-1], spca.res3[-1])
})

test_that('spca returns appropriate error for invalid assay/data',{
  ## expect error
  expect_error(spca(data = nutrimouse.mae, X = "not-an-assay", keepX = c(5,5)), class = "inv_assay")
  expect_error(spca(data = nutrimouse, X = "lipid", keepX = c(5,5)), class = "inv_data")
  expect_error(spca(X = "lipid"), class = "inv_X")
})
