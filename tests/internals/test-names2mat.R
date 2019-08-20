context("internals")

test_that(".get_xy works when Y is colData or assay name ", {
  expect_true(all((lapply(.get_xy(list(data = miniACC, X="RNASeq2GeneNorm", Y="gistict")), class))=="matrix"))
  expect_true(all((lapply(.get_xy(list(data = miniACC, X="RNASeq2GeneNorm", Y="years_to_birth")), class))=="matrix"))
})

test_that(".get_xy fails when Y is invalid", {
  expect_error(.get_xy(list(data = miniACC, X="RNASeq2GeneNorm", Y="invalid y")), class = "inv_xy")
  expect_error(.get_xy(list(data = miniACC, X="inv x", Y="years_to_birth")), class = "inv_xy")
})
