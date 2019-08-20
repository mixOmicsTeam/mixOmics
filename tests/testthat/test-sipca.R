context("sipca") ## context helps devtools::test() to categorise tests
## the 'call'slot in all results must be ignored before comparison for identicality
## as it would naturally won't be identical

test_that('sipca results for nutrimouse are the same when using either matrix or MAE class',{
  sipca.res1 <- sipca(X = liver.toxicity$gene, ncomp = 3, mode = "deflation", keepX = c(50,50,50))
  sipca.res2 <- sipca(X = "gene", data = liver.toxicity.mae, ncomp = 3, mode = "deflation", keepX = c(50,50,50))
  expect_identical(sipca.res1[-1], sipca.res2[-1])
})

test_that('assay name can be input as a variable',{
  sipca.res1 <- sipca(X = nutrimouse$lipid, ncomp = 2, keepX = c(5,5))
  sipca.res3 <- sipca(X = pca.assay, data = nutrimouse.mae, ncomp = 2, keepX = c(5,5))
  expect_identical(sipca.res1[-1], sipca.res3[-1])
})

test_that('sipca returns appropriate error for invalid assay/data',{
  ## expect error
  expect_error(sipca(data = nutrimouse.mae, X = "not-an-assay", ncomp = 2,  keepX = c(5,5)), class = "inv_assay")
  expect_error(sipca(data = nutrimouse, X = "lipid", ncomp = 2, keepX = c(5,5)), class = "inv_data")
  expect_error(sipca(X = "lipid"), class = "inv_X")
})
