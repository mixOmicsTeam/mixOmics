context("plsda")
## the tests are parameterised in helper-pls.R (any helper-*.R file in ./testthatis run
## before tests) so we can try different datasets if need be, and also to cope with
## object inputs (e.g. formla stored in a variable)
##

test_that("plsda works for 'ANY' methods",{
  
  ## suppress NZV warnings
  suppressMessages({
    plsda.res.xy <-          plsda(X = Xm_Yc, Y = Ycn )
    plsda.res.xy.ret <-      pls(X = Xm_Yc, Y = Ycn, ret.call = TRUE )
    expect_true("mixo_plsda" %in% class( plsda.res.xy))
    expect_identical(plsda.res.xy.ret$call$X, Xm_Yc)
  })
})


test_that("plsda works for 'formula NOT data' methods",{

  ## suppress NZV warnings
  ## suppress NZV warnings
  suppressMessages({
    plsda.res.formula <-          plsda(formula = Ycn~Xm_Yc)
    expect_true("mixo_plsda" %in% class(plsda.res.formula))
  })
})

test_that("plsda works for 'formula AND data' methods",{

  ## suppress NZV warnings
  suppressMessages({
    plsda.res.xy <-          plsda(X = Xm_Yc, Y = Ycn )
    plsda.res.formula.mae <- plsda(formula = f_Yc, data = mae_data)
    expect_identical(pls.res.xy[-1] ,pls.res.formula.mae[-1])
  })
})



## ------ correct error with invalid signature combination for X,Y, formula, data
test_that("plsda fails with invalid signature and produces appropriate error",{

  expect_error(plsda(X=Xm_Ya, Y=Ycn,formula = Y~X ), class = "inv_signature")
  expect_condition(plsda(formula=Y~Z), regexp = 'not found')
  expect_condition(plsda(X=NULL, Y=Yam,formula = RNASeq2GeneNorm ~ gistict, data = mae_data ), class = "inv_signature")
})

## ------ correct error with invalid assays
test_that("plsda fails with invalid assay and produces appropriate error",{

  ##---- "xy"
  expect_condition(plsda(formula = Y~X, data = mae_data ), class = "inv_XY")
  expect_condition(plsda(X = "invalidX", Y="invalidY", data = mae_data ), class = "inv_XY")

  ##---- "formula"
  expect_condition(plsda(formula = Y~X, data = mae_data ), class = "inv_XY")
  ##---- 'formula_mae'
  expect_condition(plsda(formula = wrong_LHS ~ gistict, data = mae_data ), class = "inv_XY")
})

## ------ correct error with invalid formula format
test_that("plsda fails with invalid formula formats and produces expected errors",{
  expect_condition(plsda(formula = Y~X+Y), class = "inv_sformula")
  expect_condition(plsda(formula = Y+U~X), class = "inv_sformula")
})

## ------ correct error with invalid formula elements
test_that("plsda fails with invalid formula formats and produces expected errors",{
  expect_condition(plsda(formula = Y~U), class = "simpleError")
  expect_condition(plsda(data = mae_data, formula = foo~bar), class = "inv_XY")
})

## ------ correct error with non-numeric/factor Y coldata
test_that("plsda fails with invalid Y",{
  expect_condition(plsda(X=Xa , Y=Y_inv , data = mae_data), class = "inv_XY")
})

## ------ correct error with legacy code
test_that("plsda fails with invalid Y",{
  expect_condition(plsda(Xm_Ya ,Yam , ncomp = 2), class = "defunct")
})

