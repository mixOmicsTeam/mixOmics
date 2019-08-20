# context("splsda")
# ## the tests are parameterised in helper-spls.R (any helper-*.R file in ./testthatis run
# ## before tests) so we can try different datasets if need be, and also to cope with
# ## object inputs (e.g. formla stored in a variable)
# ##
# 
# ## ------ splsda works with numeric ~ matrix
# test_that("spls produces identical 'mixo_spls' classes for designated valid signatures when Y is a column data",{
# 
#   ## suppress NZV warnings
#   suppressMessages({
#     spls.res.xy <-          splsda(keepX=c(10,10), X =Xm_Yc, Y=Ycn )
#     expect_true(class( spls.res.xy)=="mixo_spls")
# 
#     spls.res.formula <-     spls(keepX=c(10,10), formula =Ycn ~ Xm_Yc)
#     expect_equal(spls.res.xy[-1] ,spls.res.formula[-1])
# 
#     spls.res.formula.mae <- spls(keepX=c(10,10), formula = f_Yc, data = mae_data)
#     expect_equal(spls.res.xy[-1] ,spls.res.formula.mae[-1])
# 
#     spls.res.xy.mae <-      spls(keepX=c(10,10), X=Xa , Y=Yc , data = mae_data)
#     expect_equal(spls.res.formula.mae[-1] , spls.res.xy.mae[-1] )
#   })
# })
# 
# ## ------ spls works with matrix ~ matrix
# test_that("spls produces identical 'mixo_spls' classes for designated valid signatures when Y is an assay",{
#   ## suppress NZV warnings
#   suppressMessages({
#     spls.res.xy <-          spls(X =Xm_Ya, Y=Yam )
#     expect_true(class( spls.res.xy)=="mixo_spls")
# 
#     spls.res.formula <-     spls(formula =Yam ~ Xm_Ya)
#     expect_equal(spls.res.xy[-1] ,spls.res.formula[-1])
# 
#     spls.res.formula.mae <- spls(formula = f_Ya, data = mae_data)
#     expect_equal(spls.res.xy[-1] ,spls.res.formula.mae[-1])
# 
#     spls.res.xy.mae <-      spls(X=X , Y=Ya , data = mae_data)
#     expect_equal(spls.res.formula.mae[-1] , spls.res.xy.mae[-1] )
#   })
# })
# 
# ## ------ correct error with invalid signature combination for X,Y, formula, data
# test_that("spls fails with invalid signature and produces appropriate error",{
# 
#   expect_error(spls(X=Xm_Ya, Y=Ycn,formula = Y~X ))
#   expect_condition(spls(X=Y~Z ), regexp = 'must be a numeric matrix')
#   expect_condition(spls(X=NULL, Y=Yam,formula = RNASeq2GeneNorm ~ gistict, data = mae_data ), class = "inv_signature")
# })
# 
# ## ------ correct error with invalid assays
# test_that("spls fails with invalid assay and produces appropriate error",{
# 
#   ##---- "xy"
#   expect_condition(spls(formula = Y~X, data = mae_data ), class = "inv_xy")
#   expect_condition(spls(X = "invalidX", Y="invalidY", data = mae_data ), class = "inv_xy")
# 
#   expect_condition(spls(X=Xm_Ya, Y=Y_inv,data = mae_data),  class = "inv_xy")
#   ##---- "formula"
#   expect_condition(spls(formula = Y~X, data = mae_data ), class = "inv_xy")
#   ##---- 'formula_mae'
#   expect_condition(spls(formula = wrong_LHS ~ gistict, data = mae_data ), class = "inv_xy")
# })
# 
# ## ------ correct error with invalid formula format
# test_that("spls fails with invalid formula formats and produces expected errors",{
#   expect_condition(spls(formula = Y~X+Y), class = "inv_sformula")
#   expect_condition(spls(formula = Y+U~X), class = "inv_sformula")
# })
# 
# ## ------ correct error with invalid formula elements
# test_that("spls fails with invalid formula formats and produces expected errors",{
#   expect_condition(spls(formula = Y~U), class = "simpleError")
# })
# 
# ## ------ correct error with non-numeric/factor Y coldata
# test_that("spls fails with invalid Y",{
#   expect_condition(spls(X=Xa , Y=Y_inv , data = mae_data), class = "inv_xy")
# })
