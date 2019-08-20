context("ipca")

test_that('ipca results for nutrimouse are the same when using either list or MAE class (with and without quotes for Y)',{
  ipca.res1 <- ipca(X = nutrimouse$gene)
  ipca.res2 <- ipca(data = nutrimouse.mae, X = "gene")
  ## ignore call slot and expect identical
  expect_identical(ipca.res1[-1], ipca.res2[-1])
})

test_that('ipca returns appropriate error for invalid assay/data',{
  ## expect error
  expect_error(ipca(data = nutrimouse.mae, X = "not-an-assay"), class = "inv_assay")
  expect_error(ipca(data = nutrimouse, X = "lipid"), class = "inv_data")
  expect_error(ipca(X = "lipid"), class = "inv_matrix")
})
