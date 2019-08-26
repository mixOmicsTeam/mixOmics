context("block.spls")

## the tests are parameterised in helper-block.pls.R (any helper-*.R file in ./testthatis run
## before tests) so we can try different datasets if need be, and also to cope with
## object inputs (e.g. formla stored in a variable)
##

test_that("block.pls works for 'ANY' methods",{
  
  ## suppress NZV warnings
  suppressMessages({
    bpls.res.xy <-          block.pls(X = X_bpls, Y = Y_bpls )
    bpls.res.xy.ret <-      pls(X = Xm_Yc, Y = Ycn, ret.call = TRUE )
    expect_true(class( pls.res.xy) == "mixo_pls")
    expect_true(class( pls.res.xy.ret) == "mixo_pls")
    expect_identical(unclass(pls.res.xy), unclass(pls.res.xy.ret)[-1])
  })
})


test_that("block.pls works", {
   block.spls.res1 <- block.pls(X = X_bpls, Y = Y_bpls)
   block.spls.res2 <- block.pls(formula = formula_bpls_mat)
   block.spls.res3 <- block.pls(formula = formula_bpls_assay)
 })
