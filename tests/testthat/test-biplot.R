
###############################################################################
### ================================ BASIC ================================ ###
###############################################################################

test_that("(biplot:basic): pca", {
  
  data(nutrimouse)
  X <- nutrimouse$lipid[1:10,1:10]
  
  res.pca <- pca(X)
  
  pca.biplot <- biplot(res.pca)
  
  expect_true("ggplot" %in% class(pca.biplot))
})

test_that("(biplot:basic): spca", {
  
  data(nutrimouse)
  X <- nutrimouse$lipid[1:10,1:10]
  
  res.spca <- spca(X, keepX = c(3,3))
  
  spca.biplot <- biplot(res.spca)
  
  expect_true("ggplot" %in% class(spca.biplot))
})