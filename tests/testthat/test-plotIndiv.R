context("plotIndiv")

# make pca model
data(multidrug)
X <- multidrug$ABC.trans
pca.res <- pca(X, ncomp = 4, scale = TRUE)

# make rcc model
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
rcc.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

# make diablo model
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = nutrimouse$diet)
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3,
                byrow = TRUE, dimnames = list(names(data), names(data)))
diablo.res = block.plsda(X = data, indY = 3) # indY indicates where the outcome Y is in the list X

# make mint model
data(stemcells)
mint.res <- mint.splsda(X = stemcells$gene, Y = stemcells$celltype, ncomp = 2, keepX = c(10, 5),
                        study = stemcells$study)


## ------------------------------------------------------------------------ ##
## Input checking

test_that("plotIndiv throws error for invalid 'style", {
  expect_error(plotIndiv(rcc.res, ncomp = 5, style = "other"), 
               "'style' must be one of 'ggplot2', 'lattice', 'graphics' or '3d' .")
})

test_that("plotIndiv throws error for invalid 'style for MINT", {
  expect_error(plotIndiv(mint.res, ncomp = 5, style = "3d"), 
               "'style' must be one of 'ggplot2', 'lattice' or 'graphics'.")
})


test_that("plotIndiv throws error for invalid 'axes.box", {
  expect_error(plotIndiv(rcc.res, ncomp = 5, axes.box = NA), 
               "'axes.box' should be one of 'box', 'bbox' or 'both'.")
})

test_that("plotIndiv throws error for invalid 'ellipse.level", {
  expect_error(plotIndiv(rcc.res, ncomp = 5, ellipse.level = 1.5), 
               "The value taken by 'ellipse.level' must be between 0 and 1")
})

test_that("plotIndiv throws error for invalid 'legend", {
  expect_error(plotIndiv(rcc.res, ncomp = 5, legend = "legend"), 
               "'legend' must be a logical value.")
})

test_that("plotIndiv throws error for invalid 'alpha", {
  expect_error(plotIndiv(rcc.res, ncomp = 5, alpha = 1.5), 
               "The value taken by 'alpha' must be between 0 and 1")
})

test_that("plotIndiv throws error for invalid 'ncomp' for a 2D plot", {
  expect_error(plotIndiv(rcc.res, comp = 5, style = "ggplot2"), 
               "'comp' must be a numeric vector of length 2.")
})

test_that("plotIndiv throws error for invalid 'ncomp' for a 3D plot", {
  expect_error(plotIndiv(rcc.res, comp = 5, style = "3d"), 
               "'comp' must be a numeric vector of length 3.")
})

test_that("plotIndiv throws error for invalid 'ellipse'", {
  expect_error(plotIndiv(rcc.res, ellipse = "other"), 
               "'ellipse' must be either TRUE or FALSE")
})

test_that("plotIndiv throws error for invalid 'centroid'", {
  expect_error(plotIndiv(rcc.res, centroid = "other"), 
               "'centroid' must be either TRUE or FALSE")
})

test_that("plotIndiv throws error for invalid 'star'", {
  expect_error(plotIndiv(rcc.res, star = "other"), 
               "'star' must be either TRUE or FALSE")
})

test_that("plotIndiv throws error for invalid 'abline'", {
  expect_error(plotIndiv(rcc.res, abline = "other"), 
               "'abline' must be either TRUE or FALSE")
})

test_that("X.label validation works correctly", {
  
  # X.label is valid
  expect_silent(plotIndiv(rcc.res, X.label = "X-axis label"))
  
  # X.label must be a vector of length 1
  expect_error(plotIndiv(rcc.res, X.label = c("X-axis", "Extra label")), 
               "'X.label' must be a vector of length 1")
  
  # X.label must not be a non-vector (e.g., a matrix)
  expect_error(plotIndiv(rcc.res, X.label = matrix("X-axis", nrow=1, ncol=1)), 
               "'X.label' must be a vector of length 1")
})

test_that("Y.label validation works correctly", {
  
  # Y.label is valid
  expect_silent(plotIndiv(rcc.res, Y.label = "Y-axis label"))

  # Y.label must be a vector of length 1
  expect_error(plotIndiv(rcc.res, Y.label = c("Y-axis", "Extra label")), 
               "'Y.label' must be a vector of length 1")
})

test_that("Z.label validation works correctly", {
  
  # Z.label in non-3d style should raise a warning
  expect_warning(plotIndiv(rcc.res, Z.label = "Z-axis label", style = "ggplot2"), 
                 "'Z.label' is not used as style!= '3d'")
  
  # Z.label must be a vector of length 1
  expect_error(plotIndiv(rcc.res, Z.label = c("Z-axis", "Extra label"), style = "3d"), 
               "'Z.label' must be a vector of length 1")
  
  # Z.label must not be a non-vector (e.g., a data frame)
  expect_error(plotIndiv(rcc.res, Z.label = data.frame("Z-axis"), style = "3d"), 
               "'Z.label' must be a vector of length 1")
})

test_that("size.title validation", {
  
  expect_error(plotIndiv(rcc.res, size.title = -1), "'size.title' needs to be a non negative number")
  expect_error(plotIndiv(rcc.res, size.title = c(1, 2)), "'size.title' needs to be a non negative number")
  expect_silent(plotIndiv(rcc.res, size.title = 10))
})

test_that("size.title validation", {
  
  expect_error(plotIndiv(rcc.res, size.title = -1), "'size.title' needs to be a non negative number")
  expect_error(plotIndiv(rcc.res, size.title = c(1, 2)), "'size.title' needs to be a non negative number")
  expect_silent(plotIndiv(rcc.res, size.title = 10))
})

test_that("size.subtitle validation", {
  
  expect_error(plotIndiv(rcc.res, size.subtitle = -1), "'size.subtitle' needs to be a non negative number")
  expect_error(plotIndiv(rcc.res, size.subtitle = c(1, 2)), "'size.subtitle' needs to be a non negative number")
  expect_silent(plotIndiv(rcc.res, size.subtitle = 8))
})

test_that("size.xlabel validation", {
  
  expect_error(plotIndiv(rcc.res, size.xlabel = -1), "'size.xlabel' needs to be a non negative number")
  expect_error(plotIndiv(rcc.res, size.xlabel = c(1, 2)), "'size.xlabel' needs to be a non negative number")
  expect_silent(plotIndiv(rcc.res, size.xlabel = 5))
})

test_that("size.ylabel validation", {
  
  expect_error(plotIndiv(rcc.res, size.ylabel = -1), "'size.ylabel' needs to be a non negative number")
  expect_error(plotIndiv(rcc.res, size.ylabel = c(1, 2)), "'size.ylabel' needs to be a non negative number")
  expect_silent(plotIndiv(rcc.res, size.ylabel = 4))
})

test_that("size.axis validation", {
  
  expect_error(plotIndiv(rcc.res, size.axis = -1), "'size.axis' needs to be a non negative number")
  expect_error(plotIndiv(rcc.res, size.axis = c(1, 2)), "'size.axis' needs to be a non negative number")
  expect_silent(plotIndiv(rcc.res, size.axis = 3))
})

test_that("size.legend validation", {
  
  expect_error(plotIndiv(rcc.res, size.legend = -1), "'size.legend' needs to be a non negative number")
  expect_error(plotIndiv(rcc.res, size.legend = c(1, 2)), "'size.legend' needs to be a non negative number")
  expect_silent(plotIndiv(rcc.res, size.legend = 2))
})

test_that("size.legend.title validation", {
  
  expect_error(plotIndiv(rcc.res, size.legend.title = -1), "'size.legend.title' needs to be a non negative number")
  expect_error(plotIndiv(rcc.res, size.legend.title = c(1, 2)), "'size.legend.title' needs to be a non negative number")
  expect_silent(plotIndiv(rcc.res, size.legend.title = 1.5))
})

test_that("legend.position validation", {
  
  expect_error(plotIndiv(rcc.res, legend.position = "invalid"), '"legend.position" needs to be one of "bottom", "left", "right" or "top"')
  expect_silent(plotIndiv(rcc.res, legend.position = "top"))
})

test_that("point.lwd validation", {
  
  expect_error(plotIndiv(rcc.res, point.lwd = -1), "'point.lwd' needs to be a non negative number")
  expect_error(plotIndiv(rcc.res, point.lwd = c(1, 2)), "'point.lwd' needs to be a non negative number")
  expect_silent(plotIndiv(rcc.res, point.lwd = 1))
})

test_that("ind.names validation", {
  
  expect_error(plotIndiv(rcc.res, ind.names = c("Sample1", "Sample2", "Sample3")), "'ind.names' must be a character vector of length")
  expect_silent(plotIndiv(rcc.res, ind.names = rcc.res$sample_names))
})


## ------------------------------------------------------------------------ ##
## Input checking for plotIndiv.pca

# Test for inappropriate arguments
test_that("plotIndiv.pca warns for inappropriate arguments", {
  expect_warning(plotIndiv(pca.res, rep.space = "XY"), 
                 "'rep.space' is not used for PCA, sPCA, IPCA, sIPCA")
  
  expect_warning(plotIndiv(pca.res, study = "A", layout = "grid"), 
                 "'study' and 'layout' arguments are only used for MINT models")
  
  expect_warning(plotIndiv(pca.res, blocks = 2), 
                 "'blocks' argument is only used for multiblock models")
})

# Test for deprecated arguments
test_that("plotIndiv.pca warns for deprecated arguments", {
  expect_warning(plotIndiv(pca.res, col.per.group = c("red", "blue")), 
                 "'col.per.group' is deprecated, please use 'col' to specify colours for each group")
  
  expect_warning(plotIndiv(pca.res, pch.levels = c(1, 2)), 
                 "'pch.levels' is deprecated, please use 'pch' to specify point types")
  
  expect_warning(plotIndiv(pca.res, cols = c("red", "blue")), 
                 "'cols' is not a valid argument, did you mean 'col' ?")
})

## ------------------------------------------------------------------------ ##
## Input checking for plotIndiv.pls

# Test for inappropriate arguments
test_that("plotIndiv.pls warns for inappropriate arguments", {
  expect_warning(plotIndiv(rcc.res, study = "A", layout = "grid"), 
                 "'study' and 'layout' arguments are only used for MINT models")
  
  expect_warning(plotIndiv(rcc.res, blocks = 2), 
                 "'blocks' argument is only used for multiblock models")
})

# Test for deprecated arguments
test_that("plotIndiv.pls warns for deprecated arguments", {
  expect_warning(plotIndiv(rcc.res, col.per.group = c("red", "blue")), 
                 "'col.per.group' is deprecated, please use 'col' to specify colours for each group")
  
  expect_warning(plotIndiv(rcc.res, pch.levels = c(1, 2)), 
                 "'pch.levels' is deprecated, please use 'pch' to specify point types")
  
  expect_warning(plotIndiv(rcc.res, cols = c("red", "blue")), 
                 "'cols' is not a valid argument, did you mean 'col' ?")
})

## ------------------------------------------------------------------------ ##
## Input checking for plotIndiv.sgcca

# Test for inappropriate arguments
test_that("plotIndiv.sgcca warns for inappropriate arguments", {
  expect_warning(plotIndiv(diablo.res, rep.space = "XY"), 
                 "'rep.space' is not used for multiblock models, use 'blocks' to specify which data blocks to plot")
  
  expect_warning(plotIndiv(diablo.res, study = "A", layout = "grid"), 
                 "'study' and 'layout' arguments are only used for MINT models")
})

# Test for deprecated arguments
test_that("plotIndiv.sgcca warns for deprecated arguments", {
  expect_warning(plotIndiv(diablo.res, col.per.group = c("red", "blue")), 
                 "'col.per.group' is deprecated, please use 'col' to specify colours for each group")
  
  expect_warning(plotIndiv(diablo.res, pch.levels = c(1, 2)), 
                 "'pch.levels' is deprecated, please use 'pch' to specify point types")
  
  expect_warning(plotIndiv(diablo.res, cols = c("red", "blue")), 
                 "'cols' is not a valid argument, did you mean 'col' ?")
})

## ------------------------------------------------------------------------ ##
## Input checking for plotIndiv.mint

# Test for inappropriate arguments
test_that("plotIndiv.mint warns for inappropriate arguments", {
  expect_warning(plotIndiv(mint.res, rep.space = "XY"), 
                 "'rep.space' is not used for MINT models, use 'study' to specify whether to plot studies together or separately")
  
  expect_warning(plotIndiv(mint.res, blocks = 2), 
                 "'blocks' argument is only used for multiblock models")
})

# Test for deprecated arguments
test_that("plotIndiv.mint warns for deprecated arguments", {
  expect_warning(plotIndiv(mint.res, col.per.group = c("red", "blue")), 
                 "'col.per.group' is deprecated, please use 'col' to specify colours for each group")
  
  expect_warning(plotIndiv(mint.res, pch.levels = c(1, 2)), 
                 "'pch.levels' is deprecated, please use 'pch' to specify point types")
  
  expect_warning(plotIndiv(mint.res, cols = c("red", "blue")), 
                 "'cols' is not a valid argument, did you mean 'col' ?")
})
