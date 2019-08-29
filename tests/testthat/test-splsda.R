context("splsda")
## the tests are parameterised in helper-pls.R (any helper-*.R file in ./testthatis run
## before tests) so we can try different datasets if need be, and also to cope with
## object inputs (e.g. formla stored in a variable)
##

test_that("splsda works for 'ANY' methods",{
  
  ## suppress NZV warnings
  suppressMessages({
    splsda.res.xy <-          splsda(X = Xm_Yc, Y = Ycn )
    splsda.res.xy.ret <-      pls(X = Xm_Yc, Y = Ycn, ret.call = TRUE )
    expect_true("mixo_splsda" %in% class( splsda.res.xy))
    expect_identical(splsda.res.xy.ret$call$X, Xm_Yc)
  })
})


test_that("splsda works for 'formula NOT data' methods",{
  
  ## suppress NZV warnings
  ## suppress NZV warnings
  suppressMessages({
    splsda.res.formula <-          splsda(formula = Ycn~Xm_Yc)
    expect_true("mixo_splsda" %in% class(splsda.res.formula))
  })
})

test_that("splsda works for 'formula AND data' methods",{
  
  ## suppress NZV warnings
  suppressMessages({
    splsda.res.xy <-          splsda(X = Xm_Yc, Y = Ycn )
    splsda.res.formula.mae <- splsda(formula = f_Yc, data = mae_data)
    expect_identical(splsda.res.xy[-1] ,splsda.res.formula.mae[-1])
  })
})



## ------ correct error with invalid signature combination for X,Y, formula, data
test_that("splsda fails with invalid signature and produces appropriate error",{
  
  expect_error(splsda(X=Xm_Ya, Y=Ycn,formula = Y~X ), class = "inv_signature")
  expect_condition(splsda(formula=Y~Z), regexp = 'not found')
  expect_condition(splsda(X=NULL, Y=Yam,formula = RNASeq2GeneNorm ~ gistict, data = mae_data ), class = "ambiguousMethodSelection")
})

## ------ correct error with invalid assays
test_that("splsda fails with invalid assay and produces appropriate error",{
  
  ##---- "xy"
  expect_condition(splsda(formula = Y~X, data = mae_data ), class = "inv_XY")
  expect_condition(splsda(X = "invalidX", Y="invalidY", data = mae_data ), class = "inv_XY")
  
  ##---- "formula"
  expect_condition(splsda(formula = Y~X, data = mae_data ), class = "inv_XY")
  ##---- 'formula_mae'
  expect_condition(splsda(formula = wrong_LHS ~ gistict, data = mae_data ), class = "inv_XY")
})

## ------ correct error with invalid formula format
test_that("splsda fails with invalid formula formats and produces expected errors",{
  expect_condition(splsda(formula = Y~X+Y), class = "inv_formula")
  expect_condition(splsda(formula = Y+U~X), class = "inv_formula")
})

## ------ correct error with invalid formula elements
test_that("splsda fails with invalid formula formats and produces expected errors",{
  expect_condition(splsda(formula = Y~U), class = "simpleError")
  expect_condition(splsda(data = mae_data, formula = foo~bar), class = "inv_XY")
})

## ------ correct error with legacy code
test_that("splsda fails with invalid Y",{
  expect_condition(splsda(Xm_Ya ,Yam , ncomp = 2), class = "defunct")
})

