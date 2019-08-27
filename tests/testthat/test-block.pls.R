context("block.spls")

## the tests are parameterised in helper-block.pls.R (any helper-*.R file in ./testthatis run
## before tests) so we can try different datasets if need be, and also to cope with
## object inputs (e.g. formla stored in a variable)
##

test_that("block.pls works",{
  
  ## suppress NZV warnings
  suppressMessages({
    bpls.res.xy <-          block.pls(X = X_bpls, Y = Y_bpls )
    expect_true("block.pls" %in% class( bpls.res.xy))
    bpls.res.xy.ret <-      block.pls(X = X_bpls, Y = Y_bpls, ret.call = TRUE )
    expect_identical(unclass(bpls.res.xy), unclass(bpls.res.xy.ret)[-1])
    bpls.res.formula.xy <-          block.pls(formula = Y_bpls ~ X_bpls[[1]] + X_bpls[[2]])
    ## because the names will be different :/
    expect_equal(unname(bpls.res.xy$variates), unname(bpls.res.formula.xy$variates))
    bpls.res.formula.data <-          block.pls(data = data_bpls, formula = formula_bpls_assay)
    expect_identical(bpls.res.formula.data, bpls.res.xy)
  })
})


test_that("block.pls fails properly", {
   expect_error(block.pls(X = X_bpls, Y = Y_bpls, formula = foo~yu+zoo , data = data_bpls), class = "inv_signature")
  expect_error(block.pls(formula = foo~yu , data = data_bpls), class = "inv_formula")
 })
