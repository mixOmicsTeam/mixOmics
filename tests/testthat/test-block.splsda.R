context("block.splsda")

## the tests are parameterised in helper-block.pls.R (any helper-*.R file in ./testthatis run
## before tests) so we can try different datasets if need be, and also to cope with
## object inputs (e.g. formla stored in a variable)
##

test_that("block.splsda works",{
  
  ## suppress NZV warnings
  suppressMessages({
    bplsda.res.xy <-          block.splsda(X = X_bpls, Y = Y_bpls_DA)
    expect_true("block.splsda" %in% class( bplsda.res.xy))
    bplsda.res.xy.ret <-      block.splsda(X = X_bpls, Y = Y_bpls_DA, ret.call = TRUE )
    expect_equal(unclass(bplsda.res.xy), unclass(bplsda.res.xy.ret)[-1])
    bplsda.res.data <-          block.splsda(data = data_bpls, X = Xa_bpls, Y = Yc_bpls)
    expect_identical(bplsda.res.data, bplsda.res.xy)
  })
})

test_that("block.splsda fails properly", {
  
  expect_error(block.splsda(X = X_bpls, Y = Y_bpls_DA, formula = foo~yu+zoo , data = data_bpls), class = "inv_signature")
  expect_error(block.splsda(formula = foo~yu , data = data_bpls), class = "inv_formula")
})
