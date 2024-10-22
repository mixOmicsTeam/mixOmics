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

## ------------------------------------------------------------------------ ##
## Test that outputs are correct when running default style = "ggplot2"

test_that("plotIndiv works for rcc", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
  
  pl.res <- plotIndiv(nutri.res)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check coordinates
  .expect_numerically_close(pl.res$graph$data$x[1], 0.87088852)
  
  pl.res <- plotIndiv(nutri.res, rep.space= 'XY-variate', group = nutrimouse$genotype,
                      legend = TRUE)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check coordinates
  .expect_numerically_close(pl.res$graph$data$x[1], 0.8270997)
  # check groups
  expect_true(!is.null(pl.res$df$group))
  expect_equal(length(unique(pl.res$df$group)), length(unique(nutrimouse$genotype)))
})

test_that("plotIndiv works for (s)pca", {
  data(srbct)
  X <- srbct$gene
  Y <- srbct$class
  pca.srbct = pca(X, ncomp = 10, center = TRUE, scale = TRUE)
  groups <- factor(srbct$class, levels = c("RMS", "NB", "EWS", "BL"))
  pl.res <- plotIndiv(pca.srbct, group = groups, ind.names = FALSE, # plot the samples projected
            legend = TRUE, title = 'PCA on SRBCT, comp 1 - 2',
            col.per.group = c("red", "blue", "green", "black")) # onto the PCA subspace
  
  # check coordinates
  .expect_numerically_close(pl.res$graph$data$x[1], 10.13857)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check colour assignments are correct
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "RMS"]), "red")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "NB"]), "blue")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "EWS"]), "green")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "BL"]), "black")
  # check right number of samples
  expect_equal(dim(pca.srbct$X)[1], dim(pl.res$df)[1])
  # check groups
  expect_true(!is.null(pl.res$df$group))
  expect_equal(length(unique(pl.res$df$group)), length(unique(groups)))
  
})

test_that("plotIndiv works for (s)pls", {
  data(liver.toxicity)
  X <- liver.toxicity$gene
  Y <- liver.toxicity$clinic
  toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50),
                        keepY = c(10, 10, 10))
  
  pl.res <- plotIndiv(toxicity.spls, rep.space="X-variate", ind.name = FALSE,
                      group = factor(liver.toxicity$treatment$Time.Group),
                      legend.title = 'Time',
                      col.per.group = c("red", "blue", "green", "black"),
                      pch = factor(liver.toxicity$treatment$Dose.Group),
                      legend.title.pch = 'Dose',
                      legend = TRUE)
  
  # check coordinates
  .expect_numerically_close(pl.res$graph$data$x[1], 4.146771)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check colour assignments are correct
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 6]), "red")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 18]), "blue")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 24]), "green")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 48]), "black")
  # check right number of samples
  expect_equal(dim(toxicity.spls$X)[1], dim(pl.res$df)[1])
  # check groups
  expect_true(!is.null(pl.res$df$group))
  expect_equal(length(unique(pl.res$df$group)), length(unique(liver.toxicity$treatment$Time.Group)))
  
})

test_that("plotIndiv works for (s)plsda", {
  data(breast.tumors)
  X <- breast.tumors$gene.exp
  Y <- breast.tumors$sample$treatment
  splsda.breast <- splsda(X, Y,keepX=c(10,10), ncomp=2)
  
  pl.res <- plotIndiv(splsda.breast)
  # check coordinates
  .expect_numerically_close(pl.res$graph$data$x[1], -1.075222)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check colours
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "BE"]), "#F68B33")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "AF"]), "#388ECC")
  # check right number of samples
  expect_equal(dim(splsda.breast$X)[1], dim(pl.res$df)[1])
})

test_that("plotIndiv works for (s)plsda and ellipses", {
  data(srbct)
  X <- srbct$gene
  Y <- srbct$class
  srbct.splsda <- splsda(X, Y, ncomp = 10)
  groups <- factor(srbct$class, levels = c("RMS", "NB", "EWS", "BL"))
  pl.res <- plotIndiv(srbct.splsda , comp = 1:2, col.per.group = c("red", "blue", "green", "black"),
            group = groups, ind.names = FALSE,  # colour points by class
            ellipse = TRUE, # include 95% confidence ellipse for each class
            legend = TRUE, title = '(a) PLSDA with confidence ellipses')
  
  # check coordinates
  .expect_numerically_close(pl.res$graph$data$x[1], -6.83832)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check colours
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "RMS"]), "red")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "NB"]), "blue")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "EWS"]), "green")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "BL"]), "black")
  # check right number of samples
  expect_equal(dim(srbct.splsda$X)[1], dim(pl.res$df)[1])
  # check ellipses
  expect_false(is.null(pl.res$df.ellipse))
  # check groups
  expect_true(!is.null(pl.res$df$group))
  expect_equal(length(unique(pl.res$df$group)), length(unique(groups)))
})

test_that("plotIndiv works for (s)plsda and backgrounds", {
  data(srbct)
  X <- srbct$gene
  Y <- srbct$class
  srbct.splsda <- splsda(X, Y, ncomp = 10)
  groups <- factor(srbct$class, levels = c("RMS", "NB", "EWS", "BL"))
  background = background.predict(srbct.splsda, comp.predicted=2, dist = "max.dist")
  pl.res <- plotIndiv(srbct.splsda, comp = 1:2,
            group = srbct$class, ind.names = FALSE, # colour points by class
            background = background, # include prediction background for each class
            legend = TRUE, title = " (b) PLSDA with prediction background")
  
  # check coordinates
  .expect_numerically_close(pl.res$graph$data$x[1], -6.83832)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check right number of samples
  expect_equal(dim(srbct.splsda$X)[1], dim(pl.res$df)[1])
  # check groups
  expect_true(!is.null(pl.res$df$group))
  expect_equal(length(unique(pl.res$df$group)), length(unique(groups)))
})

test_that("plotIndiv works for mint.(s)plsda", {
  data(stemcells)
  res = mint.splsda(X = stemcells$gene, Y = stemcells$celltype, ncomp = 2, keepX = c(10, 5),
                    study = stemcells$study)
  pl.res <-   plotIndiv(res)
  # check coordinates
  .expect_numerically_close(pl.res$graph$data$x[1], -1.543685)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "graph"))
  # check right number of samples
  expect_equal(dim(res$X)[1], dim(pl.res$df)[1])
})

test_that("plotIndiv works for sgcca and rgcca", {
  data(nutrimouse)
  Y = unmap(nutrimouse$diet)
  data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)
  design1 = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
  nutrimouse.sgcca <- wrapper.sgcca(X = data,
                                    design = design1,
                                    penalty = c(0.3, 0.5, 1),
                                    ncomp = 3,
                                    scheme = "horst")
  
  pl.res <-  plotIndiv(nutrimouse.sgcca)
  # check coordinates
  .expect_numerically_close(pl.res$graph$data$x[1], 3.319955)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check right number of samples - here have 40 samples across 3 modalities (gene, lipid, Y)
  expect_equal(dim(nutrimouse.sgcca$X$gene)[1] + dim(nutrimouse.sgcca$X$lipid)[1] + dim(nutrimouse.sgcca$X$Y)[1], dim(pl.res$df)[1])
})

test_that("plotIndiv works for sgccda", {
  data(nutrimouse)
  Y = nutrimouse$diet
  data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
  design1 = matrix(c(0,1,0,1), ncol = 2, nrow = 2, byrow = TRUE)
  
  nutrimouse.sgccda1 <- wrapper.sgccda(X = data,
                                       Y = Y,
                                       design = design1,
                                       ncomp = 2,
                                       keepX = list(gene = c(10,10), lipid = c(15,15)),
                                       scheme = "centroid")
  pl.res <- plotIndiv(nutrimouse.sgccda1)
  # check coordinates
  .expect_numerically_close(pl.res$graph$data$x[1], 2.448086)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check right number of samples - here have 40 samples across 2 modalities (gene, lipid)
  expect_equal(dim(nutrimouse.sgccda1$X$gene)[1] + dim(nutrimouse.sgccda1$X$lipid)[1], dim(pl.res$df)[1])
  
})

## ------------------------------------------------------------------------ ##

test_that("plotIndiv.rcc works without ind.names", code = {
    data(nutrimouse)
    X <- nutrimouse$lipid
    Y <- nutrimouse$gene
    nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
    
    plotIndiv.res <- plotIndiv(nutri.res, group = nutrimouse$genotype, ind.names = FALSE, legend = TRUE)
    
    expect_is(plotIndiv.res$graph, "ggplot")
})

test_that("plotIndiv.sgcca(..., blocks = 'average') works", code = {
    data(nutrimouse)
    Y = unmap(nutrimouse$diet)
    data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)
    design1 = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
    nutrimouse.sgcca <- wrapper.sgcca(X = data,
                                      design = design1,
                                      penalty = c(0.3, 0.5, 1),
                                      ncomp = 2,
                                      scheme = "horst")
    
    # default style: one panel for each block
    plotindiv_res <- plotIndiv(nutrimouse.sgcca, blocks = c("lipid","average"))
    
    expect_true(any(grepl(pattern = "average", x = unique(plotindiv_res$df$Block))))
})

## ------------------------------------------------------------------------ ##
## Edge cases

test_that("plotIndiv.sgccda(..., blocks = 'average') works with ind.names and ell", code = {
    data("breast.TCGA")
    data = list(mrna = breast.TCGA$data.train$mrna, mirna = breast.TCGA$data.train$mirna,
                protein = breast.TCGA$data.train$protein)
    design = matrix(1, ncol = length(data), nrow = length(data),
                    dimnames = list(names(data), names(data)))
    diag(design) =  0
    # set number of variables to select, per component and per data set (this is set arbitrarily)
    list.keepX = list(mrna = rep(4, 2), mirna = rep(5,2), protein = rep(5, 2))
    TCGA.block.splsda = block.splsda(X = data, Y = breast.TCGA$data.train$subtype,
                                     ncomp = 2, keepX = list.keepX, design = design)
    blocks <- c("average", "mrna", "weighted.average")
    diablo_plot <- plotIndiv(TCGA.block.splsda, ind.names = FALSE, blocks = blocks)
    expect_true(all(unique(diablo_plot$df$Block) %in% c('average', 'Block: mrna', 'average (weighted)')))
})

test_that("plotIndiv.sgccda(..., blocks = 'average') works with ellipse=TRUE", code = {
    data("breast.TCGA")
    data = list(mrna = breast.TCGA$data.train$mrna, mirna = breast.TCGA$data.train$mirna,
                protein = breast.TCGA$data.train$protein)
    design = matrix(1, ncol = length(data), nrow = length(data),
                    dimnames = list(names(data), names(data)))
    diag(design) =  0
    # set number of variables to select, per component and per data set (this is set arbitrarily)
    list.keepX = list(mrna = rep(4, 2), mirna = rep(5,2), protein = rep(5, 2))
    TCGA.block.splsda = block.splsda(X = data, Y = breast.TCGA$data.train$subtype,
                                     ncomp = 2, keepX = list.keepX, design = design)
    blocks <- c("average", "mrna", "weighted.average")
    diablo_plot <- plotIndiv(TCGA.block.splsda, ind.names = TRUE, blocks = blocks, ellipse = TRUE)
    expect_true(all(unique(diablo_plot$df.ellipse$Block) %in% c('average', 'Block: mrna', 'average (weighted)')))
})


test_that("mint.splsda can be visualised with predicted background", {
  data(stemcells)
  X <- stemcells$gene
  Y <- stemcells$celltype
  S <- stemcells$study

  model <- mint.splsda(X = X,
                       Y = Y,
                       study = S)

  bgp <- background.predict(model, comp.predicted = 2,
                            resolution = 10)

  output <- plotIndiv(model, background = bgp, style = "ggplot2")

  expect_equal(dim(output$df), c(125, 9))
})

## ------------------------------------------------------------------------ ##
## Plotting with 'lattice' style

test_that("plotIndiv works for rcc (lattice style)", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
  
  pl.res <- plotIndiv(nutri.res, style = "lattice")
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check coordinates
  .expect_numerically_close(pl.res$df[1,1], 0.87088852)
  
  pl.res <- plotIndiv(nutri.res, rep.space= 'XY-variate', group = nutrimouse$genotype,
                      legend = TRUE, style = "lattice")
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check coordinates
  .expect_numerically_close(pl.res$df[1,1], 0.8270997)
  # check groups
  expect_true(!is.null(pl.res$df$group))
  expect_equal(length(unique(pl.res$df$group)), length(unique(nutrimouse$genotype)))
})

test_that("plotIndiv works for (s)pca (lattice style)", {
  data(srbct)
  X <- srbct$gene
  Y <- srbct$class
  pca.srbct = pca(X, ncomp = 10, center = TRUE, scale = TRUE)
  groups <- factor(srbct$class, levels = c("RMS", "NB", "EWS", "BL"))
  pl.res <- plotIndiv(pca.srbct, group = groups, ind.names = FALSE, # plot the samples projected
                      legend = TRUE, title = 'PCA on SRBCT, comp 1 - 2',
                      col.per.group = c("red", "blue", "green", "black"),
                      style = "lattice") # onto the PCA subspace
  
  # check coordinates
  .expect_numerically_close(pl.res$df[1,1], 10.13857)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check colour assignments are correct
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "RMS"]), "red")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "NB"]), "blue")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "EWS"]), "green")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "BL"]), "black")
  # check right number of samples
  expect_equal(dim(pca.srbct$X)[1], dim(pl.res$df)[1])
  # check groups
  expect_true(!is.null(pl.res$df$group))
  expect_equal(length(unique(pl.res$df$group)), length(unique(groups)))
  
})

test_that("plotIndiv works for (s)pls (lattice style)", {
  data(liver.toxicity)
  X <- liver.toxicity$gene
  Y <- liver.toxicity$clinic
  toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50),
                        keepY = c(10, 10, 10))
  
  pl.res <- plotIndiv(toxicity.spls, rep.space="X-variate", ind.name = FALSE,
                      group = factor(liver.toxicity$treatment$Time.Group),
                      legend.title = 'Time',
                      col.per.group = c("red", "blue", "green", "black"),
                      pch = factor(liver.toxicity$treatment$Dose.Group),
                      legend.title.pch = 'Dose',
                      legend = TRUE, style = "lattice")
  
  # check coordinates
  .expect_numerically_close(pl.res$df[1,1], 4.146771)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check colour assignments are correct
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 6]), "red")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 18]), "blue")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 24]), "green")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 48]), "black")
  # check right number of samples
  expect_equal(dim(toxicity.spls$X)[1], dim(pl.res$df)[1])
  # check groups
  expect_true(!is.null(pl.res$df$group))
  expect_equal(length(unique(pl.res$df$group)), length(unique(liver.toxicity$treatment$Time.Group)))
  
})

test_that("plotIndiv works for (s)plsda (lattice style)", {
  data(breast.tumors)
  X <- breast.tumors$gene.exp
  Y <- breast.tumors$sample$treatment
  splsda.breast <- splsda(X, Y,keepX=c(10,10), ncomp=2)
  
  pl.res <- plotIndiv(splsda.breast, style = "lattice")
  # check coordinates
  .expect_numerically_close(pl.res$df[1,1], -1.075222)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check colours
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "BE"]), "#F68B33")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "AF"]), "#388ECC")
  # check right number of samples
  expect_equal(dim(splsda.breast$X)[1], dim(pl.res$df)[1])
})

## ------------------------------------------------------------------------ ##
## Plotting with 'graphics' style

test_that("plotIndiv works for rcc (graphics style)", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
  
  pl.res <- plotIndiv(nutri.res, style = "graphics")
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check coordinates
  .expect_numerically_close(pl.res$df[1,1], 0.87088852)
  
  pl.res <- plotIndiv(nutri.res, rep.space= 'XY-variate', group = nutrimouse$genotype,
                      legend = TRUE, style = "lattice")
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check coordinates
  .expect_numerically_close(pl.res$df[1,1], 0.8270997)
  # check groups
  expect_true(!is.null(pl.res$df$group))
  expect_equal(length(unique(pl.res$df$group)), length(unique(nutrimouse$genotype)))
})

test_that("plotIndiv works for (s)pca (graphics style)", {
  data(srbct)
  X <- srbct$gene
  Y <- srbct$class
  pca.srbct = pca(X, ncomp = 10, center = TRUE, scale = TRUE)
  groups <- factor(srbct$class, levels = c("RMS", "NB", "EWS", "BL"))
  pl.res <- plotIndiv(pca.srbct, group = groups, ind.names = FALSE, # plot the samples projected
                      legend = TRUE, title = 'PCA on SRBCT, comp 1 - 2',
                      col.per.group = c("red", "blue", "green", "black"),
                      style = "graphics") # onto the PCA subspace
  
  # check coordinates
  .expect_numerically_close(pl.res$df[1,1], 10.13857)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check colour assignments are correct
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "RMS"]), "red")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "NB"]), "blue")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "EWS"]), "green")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "BL"]), "black")
  # check right number of samples
  expect_equal(dim(pca.srbct$X)[1], dim(pl.res$df)[1])
  # check groups
  expect_true(!is.null(pl.res$df$group))
  expect_equal(length(unique(pl.res$df$group)), length(unique(groups)))
  
})

test_that("plotIndiv works for (s)pls (graphics style)", {
  data(liver.toxicity)
  X <- liver.toxicity$gene
  Y <- liver.toxicity$clinic
  toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50),
                        keepY = c(10, 10, 10))
  
  pl.res <- plotIndiv(toxicity.spls, rep.space="X-variate", ind.name = FALSE,
                      group = factor(liver.toxicity$treatment$Time.Group),
                      legend.title = 'Time',
                      col.per.group = c("red", "blue", "green", "black"),
                      pch = factor(liver.toxicity$treatment$Dose.Group),
                      legend.title.pch = 'Dose',
                      legend = TRUE, style = "graphics")
  
  # check coordinates
  .expect_numerically_close(pl.res$df[1,1], 4.146771)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check colour assignments are correct
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 6]), "red")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 18]), "blue")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 24]), "green")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 48]), "black")
  # check right number of samples
  expect_equal(dim(toxicity.spls$X)[1], dim(pl.res$df)[1])
  # check groups
  expect_true(!is.null(pl.res$df$group))
  expect_equal(length(unique(pl.res$df$group)), length(unique(liver.toxicity$treatment$Time.Group)))
  
})

test_that("plotIndiv works for (s)plsda (graphics style)", {
  data(breast.tumors)
  X <- breast.tumors$gene.exp
  Y <- breast.tumors$sample$treatment
  splsda.breast <- splsda(X, Y,keepX=c(10,10), ncomp=2)
  
  pl.res <- plotIndiv(splsda.breast, style = "graphics")
  # check coordinates
  .expect_numerically_close(pl.res$df[1,1], -1.075222)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check colours
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "BE"]), "#F68B33")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "AF"]), "#388ECC")
  # check right number of samples
  expect_equal(dim(splsda.breast$X)[1], dim(pl.res$df)[1])
})

## ------------------------------------------------------------------------ ##
## Plotting with '3d' style

test_that("plotIndiv works for rcc (3d style)", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
  
  pl.res <- plotIndiv(nutri.res, style = "3d")
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check coordinates
  .expect_numerically_close(pl.res$df[1,1], 0.87088852)
  
  pl.res <- plotIndiv(nutri.res, rep.space= 'XY-variate', group = nutrimouse$genotype,
                      legend = TRUE, style = "3d")
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check coordinates
  .expect_numerically_close(pl.res$df[1,1], 0.8270997)
  # check groups
  expect_true(!is.null(pl.res$df$group))
  expect_equal(length(unique(pl.res$df$group)), length(unique(nutrimouse$genotype)))
})

test_that("plotIndiv works for (s)pca (3d style)", {
  data(srbct)
  X <- srbct$gene
  Y <- srbct$class
  pca.srbct = pca(X, ncomp = 10, center = TRUE, scale = TRUE)
  groups <- factor(srbct$class, levels = c("RMS", "NB", "EWS", "BL"))
  pl.res <- plotIndiv(pca.srbct, group = groups, ind.names = FALSE, # plot the samples projected
                      legend = TRUE, title = 'PCA on SRBCT, comp 1 - 2',
                      col.per.group = c("red", "blue", "green", "black"),
                      style = "3d") # onto the PCA subspace
  
  # check coordinates
  .expect_numerically_close(pl.res$df[1,1], 10.13857)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check colour assignments are correct
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "RMS"]), "red")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "NB"]), "blue")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "EWS"]), "green")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "BL"]), "black")
  # check right number of samples
  expect_equal(dim(pca.srbct$X)[1], dim(pl.res$df)[1])
  # check groups
  expect_true(!is.null(pl.res$df$group))
  expect_equal(length(unique(pl.res$df$group)), length(unique(groups)))
  
})

test_that("plotIndiv works for (s)pls (3d style)", {
  data(liver.toxicity)
  X <- liver.toxicity$gene
  Y <- liver.toxicity$clinic
  toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50),
                        keepY = c(10, 10, 10))
  
  # expect error when passing numbers as pch argument, informative message about what pch levels can be passed for 3d plot
  expect_error(plotIndiv(toxicity.spls, rep.space="X-variate", ind.name = FALSE,
                      group = factor(liver.toxicity$treatment$Time.Group),
                      legend.title = 'Time',
                      col.per.group = c("red", "blue", "green", "black"),
                      pch = factor(liver.toxicity$treatment$Dose.Group),
                      legend.title.pch = 'Dose',
                      legend = TRUE, style = "3d"),
               "pch' must be a simple character or character vector from {'sphere', 'tetra', 'cube', 'octa', 'icosa', 'dodeca'}.",
               fixed = TRUE)
  
  # plot runs when the correct pch values are used for 3d plot, though not sure if it will be used much as not very informative plot
  pchs <- factor(liver.toxicity$treatment$Dose.Group)
  levels(pchs) <- c("sphere", "tetra", "octa", "icosa")
  pl.res <- plotIndiv(toxicity.spls, rep.space="X-variate", ind.name = FALSE,
                      group = factor(liver.toxicity$treatment$Time.Group),
                      legend.title = 'Time',
                      col.per.group = c("red", "blue", "green", "black"),
                      pch = pchs,
                      legend.title.pch = 'Dose',
                      legend = TRUE, style = "3d")
  
  # check coordinates
  .expect_numerically_close(pl.res$df[1,1], 4.146771)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check colour assignments are correct
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 6]), "red")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 18]), "blue")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 24]), "green")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 48]), "black")
  # check right number of samples
  expect_equal(dim(toxicity.spls$X)[1], dim(pl.res$df)[1])
  # check groups
  expect_true(!is.null(pl.res$df$group))
  expect_equal(length(unique(pl.res$df$group)), length(unique(liver.toxicity$treatment$Time.Group)))
  
})

test_that("plotIndiv works for (s)plsda (3d style)", {
  data(breast.tumors)
  X <- breast.tumors$gene.exp
  Y <- breast.tumors$sample$treatment
  splsda.breast <- splsda(X, Y,keepX=c(10,10), ncomp=3)
  
  pl.res <- plotIndiv(splsda.breast, style = "3d")
  # check coordinates
  .expect_numerically_close(pl.res$df[1,1], -1.075222)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check colours
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "BE"]), "#F68B33")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "AF"]), "#388ECC")
  # check right number of samples
  expect_equal(dim(splsda.breast$X)[1], dim(pl.res$df)[1])
})

unlink(list.files(pattern = "*.pdf"))

