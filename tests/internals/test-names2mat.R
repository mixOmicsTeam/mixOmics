context("internals")

test_that(".names2mat works when Y is colData or assay name ", {
  expect_true(all((lapply(.names2mat(list(data = miniACC, X="RNASeq2GeneNorm", Y="gistict")), class))=="matrix"))
  expect_true(all((lapply(.names2mat(list(data = miniACC, X="RNASeq2GeneNorm", Y="years_to_birth")), class))=="matrix"))
})

test_that(".names2mat fails when Y is invalid", {
  expect_error(.names2mat(list(data = miniACC, X="RNASeq2GeneNorm", Y="invalid y")), class = "inv_xy")
  expect_error(.names2mat(list(data = miniACC, X="inv x", Y="years_to_birth")), class = "inv_xy")
})
