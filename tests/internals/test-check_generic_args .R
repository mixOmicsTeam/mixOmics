context("internals")
library(MultiAssayExperiment)

test_that("check_generic_args works when expectattions match arguments",{
 ## Expect= c("xy", "formula", "formula_mae", "xy_mae")
  X=assay(miniACC,1)
  Y=assay(miniACC,2)
  ##---- "xy"
  res <- try(check_generic_args(X =X, Y=Y,formula = NULL, data = NULL, Expect = "xy" ))
  expect_true(class(res)=="list")
  expect_true(is.numeric(res$Y) & is.matrix(res$Y))
  expect_true(is.matrix(res$X))
  ##---- "formula"
  res <- try(check_generic_args(X=NULL, Y=NULL,formula = Y~X, data = NULL, Expect = "formula" ))
  expect_true(class(res)=="list")
  expect_true(is.numeric(res$Y) & is.matrix(res$Y))
  expect_true(is.matrix(res$X))
  ##---- 'formula_mae'
  res <- try(check_generic_args(X=NULL, Y=NULL,formula = RNASeq2GeneNorm ~ gistict, data = miniACC, Expect = "formula_mae" ))
  expect_true(class(res)=="list")
  expect_true(is.numeric(res$Y) & is.matrix(res$Y))
  expect_true(is.matrix(res$X))
  ##---- 'xy_mae'
  ### quoted
  res <- try(check_generic_args(X="gistict", Y="RNASeq2GeneNorm",formula = NULL, data = miniACC, Expect = "xy_mae" ))
  expect_true(class(res)=="list")
  expect_true(is.numeric(res$Y) & is.matrix(res$Y))
  expect_true(is.matrix(res$X))
  ###  symbol
  res <- try(check_generic_args(X=gistict, Y=RNASeq2GeneNorm,formula = NULL, data = miniACC, Expect = "xy_mae" ))
  expect_true(class(res)=="list")
  expect_true(is.numeric(res$Y) & is.matrix(res$Y))
  expect_true(is.matrix(res$X))
})

test_that("check_generic_args fails when expectattions do not match arguments",{
  ## Expect= c("xy", "formula", "formula_mae", "xy_mae")
  X=assay(miniACC,1)
  Y=assay(miniACC,2)
  ##---- "xy"
  expect_condition(check_generic_args(X =X, Y=Y,formula = Y~X, data = NULL, Expect = "xy" ), class = "args_conflict")
  expect_condition(check_generic_args(X =X, Y=Y,formula = NULL, data = miniACC, Expect = "xy" ), class = "args_conflict")
  ##---- "formula"
  expect_condition(check_generic_args(X=NULL, Y=NULL,formula = Y~X, data = miniACC, Expect = "formula" ), class = "args_conflict")
  ##---- 'formula_mae'
  expect_condition(check_generic_args(X=NULL, Y=Y,formula = RNASeq2GeneNorm ~ gistict, data = miniACC, Expect = "formula_mae" ), class = "args_conflict")
  ##---- 'xy_mae'
  expect_condition(check_generic_args(X=X, Y=Y,formula = RNASeq2GeneNorm ~ gistict, data = miniACC, Expect = "xy_mae" ), class = "args_conflict")
  expect_condition(check_generic_args(X=gistictz, Y=RNASeq2GeneNorm,formula = NULL, data = miniACC, Expect = "xy_mae" ), class = "invalid_X")
  expect_condition(check_generic_args(X=gistict, Y=invalidY,formula = NULL, data = miniACC, Expect = "xy_mae" ), class = "invalid_Y")
})

