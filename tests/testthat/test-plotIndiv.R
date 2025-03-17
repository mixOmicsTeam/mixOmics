context("plotIndiv")

## ------------------------------------------------------------------------ ##
## Input checking

test_that("plotIndiv throws error for invalid 'style", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
  expect_error(plotIndiv(nutri.res, ncomp = 5, style = "other"), 
               "'style' must be one of 'ggplot2', 'lattice', 'graphics' or '3d' .")
})

test_that("plotIndiv throws error for invalid 'axes.box", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
  expect_error(plotIndiv(nutri.res, ncomp = 5, axes.box = NA), 
               "'axes.box' should be one of 'box', 'bbox' or 'both'.")
})

test_that("plotIndiv throws error for invalid 'ellipse.level", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
  expect_error(plotIndiv(nutri.res, ncomp = 5, ellipse.level = 1.5), 
               "The value taken by 'ellipse.level' must be between 0 and 1")
})

test_that("plotIndiv throws error for invalid 'legend", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
  expect_error(plotIndiv(nutri.res, ncomp = 5, legend = "legend"), 
               "'legend' must be a logical value.")
})

test_that("plotIndiv throws error for invalid 'alpha", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
  expect_error(plotIndiv(nutri.res, ncomp = 5, alpha = 1.5), 
               "The value taken by 'alpha' must be between 0 and 1")
})

test_that("plotIndiv throws error for invalid 'ncomp' for a 2D plot", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 5, lambda1 = 0.064, lambda2 = 0.008)
  expect_error(plotIndiv(nutri.res, comp = 5, style = "ggplot2"), 
               "'comp' must be a numeric vector of length 2.")
})

test_that("plotIndiv throws error for invalid 'ncomp' for a 3D plot", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 5, lambda1 = 0.064, lambda2 = 0.008)
  expect_error(plotIndiv(nutri.res, comp = 5, style = "3d"), 
               "'comp' must be a numeric vector of length 3.")
})

test_that("plotIndiv throws error for invalid 'ellipse'", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 5, lambda1 = 0.064, lambda2 = 0.008)
  expect_error(plotIndiv(nutri.res, ellipse = "other"), 
               "'ellipse' must be either TRUE or FALSE")
})

test_that("plotIndiv throws error for invalid 'centroid'", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 5, lambda1 = 0.064, lambda2 = 0.008)
  expect_error(plotIndiv(nutri.res, centroid = "other"), 
               "'centroid' must be either TRUE or FALSE")
})

test_that("plotIndiv throws error for invalid 'star'", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 5, lambda1 = 0.064, lambda2 = 0.008)
  expect_error(plotIndiv(nutri.res, star = "other"), 
               "'star' must be either TRUE or FALSE")
})

test_that("plotIndiv throws error for invalid 'abline'", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 5, lambda1 = 0.064, lambda2 = 0.008)
  expect_error(plotIndiv(nutri.res, abline = "other"), 
               "'abline' must be either TRUE or FALSE")
})

test_that("X.label validation works correctly", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
  
  # X.label is valid
  expect_silent(plotIndiv(nutri.res, X.label = "X-axis label"))
  
  # X.label must be a vector of length 1
  expect_error(plotIndiv(nutri.res, X.label = c("X-axis", "Extra label")), 
               "'X.label' must be a vector of length 1")
  
  # X.label must not be a non-vector (e.g., a matrix)
  expect_error(plotIndiv(nutri.res, X.label = matrix("X-axis", nrow=1, ncol=1)), 
               "'X.label' must be a vector of length 1")
})

test_that("Y.label validation works correctly", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
  
  # Y.label is valid
  expect_silent(plotIndiv(nutri.res, Y.label = "Y-axis label"))

  # Y.label must be a vector of length 1
  expect_error(plotIndiv(nutri.res, Y.label = c("Y-axis", "Extra label")), 
               "'Y.label' must be a vector of length 1")
})

test_that("Z.label validation works correctly", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
  
  # Z.label in non-3d style should raise a warning
  expect_warning(plotIndiv(nutri.res, Z.label = "Z-axis label", style = "ggplot2"), 
                 "'Z.label' is not used as style!= '3d'")
  
  # Z.label must be a vector of length 1
  expect_error(plotIndiv(nutri.res, Z.label = c("Z-axis", "Extra label"), style = "3d"), 
               "'Z.label' must be a vector of length 1")
  
  # Z.label must not be a non-vector (e.g., a data frame)
  expect_error(plotIndiv(nutri.res, Z.label = data.frame("Z-axis"), style = "3d"), 
               "'Z.label' must be a vector of length 1")
})

test_that("size.title validation", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
  
  expect_error(plotIndiv(nutri.res, size.title = -1), "'size.title' needs to be a non negative number")
  expect_error(plotIndiv(nutri.res, size.title = c(1, 2)), "'size.title' needs to be a non negative number")
  expect_silent(plotIndiv(nutri.res, size.title = 10))
})

test_that("size.title validation", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
  
  expect_error(plotIndiv(nutri.res, size.title = -1), "'size.title' needs to be a non negative number")
  expect_error(plotIndiv(nutri.res, size.title = c(1, 2)), "'size.title' needs to be a non negative number")
  expect_silent(plotIndiv(nutri.res, size.title = 10))
})

test_that("size.subtitle validation", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
  
  expect_error(plotIndiv(nutri.res, size.subtitle = -1), "'size.subtitle' needs to be a non negative number")
  expect_error(plotIndiv(nutri.res, size.subtitle = c(1, 2)), "'size.subtitle' needs to be a non negative number")
  expect_silent(plotIndiv(nutri.res, size.subtitle = 8))
})

test_that("size.xlabel validation", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
  
  expect_error(plotIndiv(nutri.res, size.xlabel = -1), "'size.xlabel' needs to be a non negative number")
  expect_error(plotIndiv(nutri.res, size.xlabel = c(1, 2)), "'size.xlabel' needs to be a non negative number")
  expect_silent(plotIndiv(nutri.res, size.xlabel = 5))
})

test_that("size.ylabel validation", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
  
  expect_error(plotIndiv(nutri.res, size.ylabel = -1), "'size.ylabel' needs to be a non negative number")
  expect_error(plotIndiv(nutri.res, size.ylabel = c(1, 2)), "'size.ylabel' needs to be a non negative number")
  expect_silent(plotIndiv(nutri.res, size.ylabel = 4))
})

test_that("size.axis validation", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
  
  expect_error(plotIndiv(nutri.res, size.axis = -1), "'size.axis' needs to be a non negative number")
  expect_error(plotIndiv(nutri.res, size.axis = c(1, 2)), "'size.axis' needs to be a non negative number")
  expect_silent(plotIndiv(nutri.res, size.axis = 3))
})

test_that("size.legend validation", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
  
  expect_error(plotIndiv(nutri.res, size.legend = -1), "'size.legend' needs to be a non negative number")
  expect_error(plotIndiv(nutri.res, size.legend = c(1, 2)), "'size.legend' needs to be a non negative number")
  expect_silent(plotIndiv(nutri.res, size.legend = 2))
})

test_that("size.legend.title validation", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
  
  expect_error(plotIndiv(nutri.res, size.legend.title = -1), "'size.legend.title' needs to be a non negative number")
  expect_error(plotIndiv(nutri.res, size.legend.title = c(1, 2)), "'size.legend.title' needs to be a non negative number")
  expect_silent(plotIndiv(nutri.res, size.legend.title = 1.5))
})

test_that("legend.position validation", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
  
  expect_error(plotIndiv(nutri.res, legend.position = "invalid"), '"legend.position" needs to be one of "bottom", "left", "right" or "top"')
  expect_silent(plotIndiv(nutri.res, legend.position = "top"))
})

test_that("point.lwd validation", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
  
  expect_error(plotIndiv(nutri.res, point.lwd = -1), "'point.lwd' needs to be a non negative number")
  expect_error(plotIndiv(nutri.res, point.lwd = c(1, 2)), "'point.lwd' needs to be a non negative number")
  expect_silent(plotIndiv(nutri.res, point.lwd = 1))
})

test_that("ind.names validation", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
  
  expect_error(plotIndiv(nutri.res, ind.names = c("Sample1", "Sample2", "Sample3")), "'ind.names' must be a character vector of length")
  expect_silent(plotIndiv(nutri.res, ind.names = nutri.res$sample_names))
})