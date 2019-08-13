context("pls")
## the tests are parameterised in helper-pls.R (any helper-*.R file in ./testthatis run
## before tests) so we can try different datasets if need be, and also to cope with
## object inputs (e.g. formla stored in a variable)
##

## ------ pls works with numeric ~ matrix
test_that("pls produces identical 'mixo_pls' classes for designated valid signatures when Y is a column data",{

  ## suppress NZV warnings
  suppressMessages({
    pls.res.xy <-          pls(X =Xm_Yc, Y=Ycn )
    expect_true(class( pls.res.xy)=="mixo_pls")

    pls.res.formula <-     pls(formula =Ycn ~ Xm_Yc)
    expect_equal(pls.res.xy[-1] ,pls.res.formula[-1])

    pls.res.formula.mae <- pls(formula = f_Yc, data = mae_data)
    expect_equal(pls.res.xy[-1] ,pls.res.formula.mae[-1])

    pls.res.xy.mae <-      pls(X=Xa , Y=Yc , data = mae_data)
    expect_equal(pls.res.formula.mae[-1] , pls.res.xy.mae[-1] )
    })
})

## ------ pls works with matrix ~ matrix
test_that("pls produces identical 'mixo_pls' classes for designated valid signatures when Y is an assay",{
  ## suppress NZV warnings
  suppressMessages({
    pls.res.xy <-          pls(X =Xm_Ya, Y=Yam )
    expect_true(class( pls.res.xy)=="mixo_pls")

    pls.res.formula <-     pls(formula =Yam ~ Xm_Ya)
    expect_equal(pls.res.xy[-1] ,pls.res.formula[-1])

    pls.res.formula.mae <- pls(formula = f_Ya, data = mae_data)
    expect_equal(pls.res.xy[-1] ,pls.res.formula.mae[-1])

    pls.res.xy.mae <-      pls(X=Xa , Y=Ya , data = mae_data)
    expect_equal(pls.res.formula.mae[-1] , pls.res.xy.mae[-1] )
  })
})

## ------ correct error with invalid signature combination for X,Y, formula, data
test_that("pls fails with invalid signature and produces appropriate error",{

  expect_error(pls(X=Xm_Ya, Y=Ycn,formula = Y~X ))
  expect_condition(pls(X=Y~Z), regexp = 'must be a numeric matrix')
  expect_condition(pls(X=NULL, Y=Yam,formula = RNASeq2GeneNorm ~ gistict, data = mae_data ), class = "inv_signature")
})

## ------ correct error with invalid assays
test_that("pls fails with invalid assay and produces appropriate error",{

  ##---- "xy"
  expect_condition(pls(formula = Y~X, data = mae_data ), class = "inv_xy")
  expect_condition(pls(X = "invalidX", Y="invalidY", data = mae_data ), class = "inv_xy")

  expect_condition(pls(X=Xm_Ya, Y=Y_inv,data = mae_data),  class = "inv_xy")
  ##---- "formula"
  expect_condition(pls(formula = Y~X, data = mae_data ), class = "inv_xy")
  ##---- 'formula_mae'
  expect_condition(pls(formula = wrong_LHS ~ gistict, data = mae_data ), class = "inv_xy")
})

## ------ correct error with invalid formula format
test_that("pls fails with invalid formula formats and produces expected errors",{
  expect_condition(pls(formula = Y~X+Y), class = "inv_sformula")
  expect_condition(pls(formula = Y+U~X), class = "inv_sformula")
})

## ------ correct error with invalid formula elements
test_that("pls fails with invalid formula formats and produces expected errors",{
  expect_condition(pls(formula = Y~U), class = "simpleError")
})

## ------ correct error with non-numeric/factor Y coldata
test_that("pls fails with invalid Y",{
  expect_condition(pls(X=Xa , Y=Y_inv , data = mae_data), class = "inv_xy")
})
