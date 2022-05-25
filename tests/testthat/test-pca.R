###############################################################################
### ============================ MISSING TESTS ============================ ###
###############################################################################



###############################################################################
### ============================ GROUND TRUTHS ============================ ###
###############################################################################

Test.Data <- readRDS(system.file("testdata", "testdata-pca.rda", package = "mixOmics"))
Testable.Components <- Test.Data$tc
Ground.Truths <- Test.Data$gt

###############################################################################
### ================================ BASIC ================================ ###
###############################################################################


test_that("(pca:basic): multidrug", {
  
  testable.components <- Testable.Components$basic.pca
  GT <- Ground.Truths$basic.pca
  
  data(multidrug)
  X <- multidrug$ABC.trans[1:10, 1:10]
  
  res.pca <- pca(X)
  
  invisible(capture.output(TT <- dput(res.pca[testable.components])))
    
  expect_equal(TT, GT)
  
})


###############################################################################
### ================================ DATA ================================= ###
###############################################################################


test_that("(pca:data): srbct", {

  testable.components <- Testable.Components$srbct.pca
  GT <- Ground.Truths$srbct.pca

  data(srbct)
  X <- srbct$gene[1:10, 1:10]

  res.pca <- pca(X)

  invisible(capture.output(TT <- dput(res.pca[testable.components])))

  expect_equal(TT, GT)

})


test_that("(pca:data): liver.toxicity", {

  testable.components <- Testable.Components$liver.toxicity.pca
  GT <- Ground.Truths$liver.toxicity.pca

  data(liver.toxicity)
  X <- liver.toxicity$gene[1:10, 1:10]

  res.pca <- pca(X)

  invisible(capture.output(TT <- dput(res.pca[testable.components])))

  expect_equal(TT, GT)

})


test_that("(pca:data): nutrimouse", {

  testable.components <- Testable.Components$nutrimouse.pca
  GT <- Ground.Truths$nutrimouse.pca

  data(nutrimouse)
  X <- nutrimouse$lipid[1:10, 1:10]

  res.pca <- pca(X)

  invisible(capture.output(TT <- dput(res.pca[testable.components])))

  expect_equal(TT, GT)

})


###############################################################################
### ============================== PARAMETER ============================== ###
###############################################################################


test_that("(pca:parameter): ncomp", {

  testable.components <- Testable.Components$ncomp.pca
  GT <- Ground.Truths$ncomp.pca

  data(nutrimouse)
  X <- nutrimouse$lipid[1:10, 1:10]

  res.pca <- pca(X, ncomp = 3)

  invisible(capture.output(TT <- dput(res.pca[testable.components])))

  expect_equal(TT, GT)
})


test_that("(pca:parameter): center", {

  testable.components <- Testable.Components$center.pca
  GT <- Ground.Truths$center.pca

  data(nutrimouse)
  X <- nutrimouse$lipid[1:10, 1:10]

  res.pca <- pca(X, center = F)

  invisible(capture.output(TT <- dput(res.pca[testable.components])))

  expect_equal(TT, GT)
})


test_that("(pca:parameter): scale", {

  testable.components <- Testable.Components$scale.pca
  GT <- Ground.Truths$scale.pca

  data(nutrimouse)
  X <- nutrimouse$lipid[1:10, 1:10]

  res.pca <- pca(X, scale = T)

  invisible(capture.output(TT <- dput(res.pca[testable.components])))

  expect_equal(TT, GT)
})


test_that("(pca:parameter): logratio", {

  testable.components <- Testable.Components$logratio.pca
  GT <- Ground.Truths$logratio.pca

  data(nutrimouse)
  X <- nutrimouse$lipid[1:10, 1:10]

  X <- X + 0.01
  res.pca <- pca(X, logratio = "CLR")

  invisible(capture.output(TT <- dput(res.pca[testable.components])))

  expect_equal(TT, GT)
})


test_that("(pca:parameter): ilr.offset", {

  testable.components <- Testable.Components$ilr.offset.pca
  GT <- Ground.Truths$ilr.offset.pca

  data(nutrimouse)
  X <- nutrimouse$lipid[1:10, 1:10]

  res.pca <- pca(X, ilr.offset = 0.5)

  invisible(capture.output(TT <- dput(res.pca[testable.components])))

  expect_equal(TT, GT)
})


test_that("(pca:parameter): multilevel", {

  testable.components <- Testable.Components$multilevel.pca
  GT <- Ground.Truths$multilevel.pca

  data(vac18)
  X <- vac18$genes[1:10, 1:10]
  ml <- c(rep(1,2), rep(2,2), rep(3,2), rep(4,2), rep(5,2))

  res.pca <- pca(X, multilevel = ml)

  invisible(capture.output(TT <- dput(res.pca[testable.components])))

  expect_equal(TT, GT)
})


###############################################################################
### ================================ ERROR ================================ ###
###############################################################################


test_that("(pca:error): catches all invalid values of ncomp", {
  
  data(multidrug)
  X <- multidrug$ABC.trans[1:10, 1:10]
  
  expect_error(pca(X, ncomp = "number"),
               "`ncomp` must be numeric",
               fixed=T)
  
  expect_error(pca(X, ncomp = -1),
               "invalid value for 'ncomp'.",
               fixed=T)
  
  expect_error(pca(X, ncomp = 20),
               "use smaller 'ncomp'",
               fixed=T)
})


test_that("(pca:error): catches all invalid values of center", {
  
  data(multidrug)
  X <- multidrug$ABC.trans[1:10, 1:10]
  
  expect_error(pca(X, center = "number"),
               "'center' should be either a logical value or a numeric vector of length equal to the number of columns of 'X'.",
               fixed=T)
  
  expect_error(pca(X, center = c(2,3,4)),
               "'center' should be either a logical value or a numeric vector of length equal to the number of columns of 'X'.",
               fixed=T)
})


test_that("(pca:error): catches all invalid values of scale", {
  
  data(multidrug)
  X <- multidrug$ABC.trans[1:10, 1:10]
  
  expect_error(pca(X, scale = "number"),
               "'scale' should be either a logical value or a numeric vector of length equal to the number of columns of 'X'.",
               fixed=T)
  
  expect_error(pca(X, scale = c(2,3,4)),
               "'scale' should be either a logical value or a numeric vector of length equal to the number of columns of 'X'.",
               fixed=T)
})


test_that("(pca:error): catches when multilevel has differing length to nrow(X)", {
  
  data(vac18)
  X <- vac18$genes[1:9, 1:10] # nrow(X) = 9
  ml <- c(rep(1,2), rep(2,2), rep(3,2), rep(4,2), rep(5,2)) # length(ml) = 10
  
  expect_error(pca(X, multilevel = ml),
               "unequal number of rows in 'X' and 'multilevel'.",
               fixed=T)
})
