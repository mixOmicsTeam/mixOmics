test_that(".getXY works for all combinations of X,Y,data, and formula", {
  ## formula=Y_matrix~X_matrix === Y=Y_matrix, X=X_matrix
  res1 <- .getXY(list(formula = Yam~Xm_Ya))
  res2 <- .getXY(list(Y=Yam, X=Xm_Ya))
  expect_identical(res1, res2)
  ## formula=Y_assay~X_assay, data=MAE  === Y=Y_assay, X=X_assay, data=MAE
  res1 <- .getXY(list( formula = f_Ya, data=mae_data))
  res2 <- .getXY(list(X=Xa, Y=Ya,  data=mae_data))
  expect_identical(res1, res2)
  ## formula=Y_colData~X_assay, data=MAE  === Y=colData(MAE)[,"Y_colData"], X=assay(MAE, Xa)
  res1 <- .getXY(list( formula = f_Yc, data=mae_data))
  res2 <- .getXY(list(X=Xm_Yc, Y=Ycn))
  expect_identical(res1, res2)
  ## Y=Y_colData, X=X_assay, data=MAE  === Y=colData(MAE)[,"Y_colData"], X=assay(MAE, Xa)
  res1 <- .getXY(list( X=Xa, Y=Yc, data=mae_data))
  res2 <- .getXY(list(X=Xm_Yc, Y=Ycn))
  expect_identical(res1, res2)
  ## Y=Y_assay, X=X_assay, data=MAE  === Y=assay(MAE, Y_assay), X=assay(MAE, Xa)
  res1 <- .getXY(list( X=Xa, Y=Ya, data=mae_data))
  res2 <- .getXY(list(X=Xm_Ya, Y=Yam))
  expect_identical(res1, res2)
})

