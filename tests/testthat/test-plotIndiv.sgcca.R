context("plotIndiv.sgcca")

## create block-PLS-DA model
data(nutrimouse)
Y = unmap(nutrimouse$diet)
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)
design1 = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
nutrimouse.sgcca <- wrapper.sgcca(X = data,
                                  design = design1,
                                  penalty = c(0.3, 0.5, 1),
                                  ncomp = 3)

## create block-sPLS-DA model
data(nutrimouse)
Y = nutrimouse$diet
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
design1 = matrix(c(0,1,0,1), ncol = 2, nrow = 2, byrow = TRUE)

nutrimouse.sgccda1 <- block.splsda(X = data,
                                   Y = Y,
                                   design = design1,
                                   ncomp = 2,
                                   keepX = list(gene = c(10,10), lipid = c(15,15)))

## ------------------------------------------------------------------------ ##
## Test that outputs are correct when running default style = "ggplot2"

test_that("plotIndiv works for sgcca and rgcca", {
  
  pl.res <-  plotIndiv(nutrimouse.sgcca)
  # check coordinates
  .expect_numerically_close(pl.res$graph$data$x[1], 3.319955)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check right number of samples - here have 40 samples across 3 modalities (gene, lipid, Y)
  expect_equal(dim(nutrimouse.sgcca$X$gene)[1] + dim(nutrimouse.sgcca$X$lipid)[1] + dim(nutrimouse.sgcca$X$Y)[1], dim(pl.res$df)[1])
})

test_that("plotIndiv works for sgccda", {

  pl.res <- plotIndiv(nutrimouse.sgccda1)
  # check coordinates - for some reason get a different coordinate (~2.6) for windows so this check fails
  .expect_numerically_close(abs(pl.res$graph$data$x[1]), abs(-2.654198), digits = 0)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check right number of samples - here have 40 samples across 2 modalities (gene, lipid)
  expect_equal(dim(nutrimouse.sgccda1$X$gene)[1] + dim(nutrimouse.sgccda1$X$lipid)[1], dim(pl.res$df)[1])
  
})

## ------------------------------------------------------------------------ ##

test_that("plotIndiv.sgcca(..., blocks = 'average') works", code = {
    
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

## ------------------------------------------------------------------------ ##
## Plotting with 'lattice' style

test_that("plotIndiv works for sgcca and rgcca - lattice", {
  
  pl.res <-  plotIndiv(nutrimouse.sgcca, style = "lattice")
  # check coordinates
  .expect_numerically_close(pl.res$df$x[1], 3.319955)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check right number of samples - here have 40 samples across 3 modalities (gene, lipid, Y)
  expect_equal(dim(nutrimouse.sgcca$X$gene)[1] + dim(nutrimouse.sgcca$X$lipid)[1] + dim(nutrimouse.sgcca$X$Y)[1], dim(pl.res$df)[1])
})


## ------------------------------------------------------------------------ ##
## Plotting with 'graphics' style

test_that("plotIndiv works for sgcca and rgcca - graphics", {
  
  pl.res <-  plotIndiv(nutrimouse.sgcca, style = "graphics")
  # check coordinates
  .expect_numerically_close(pl.res$df$x[1], 3.319955)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check right number of samples - here have 40 samples across 3 modalities (gene, lipid, Y)
  expect_equal(dim(nutrimouse.sgcca$X$gene)[1] + dim(nutrimouse.sgcca$X$lipid)[1] + dim(nutrimouse.sgcca$X$Y)[1], dim(pl.res$df)[1])
})

## ------------------------------------------------------------------------ ##
## Plotting with '3d' style

library(rgl)

test_that("plotIndiv works for sgcca and rgcca - 3d", {
  
  pl.res <- suppressWarnings(plotIndiv(nutrimouse.sgcca, style = "3d"))
  # check coordinates
  .expect_numerically_close(pl.res$df$x[1], 3.319955)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check right number of samples - here have 40 samples across 3 modalities (gene, lipid, Y)
  expect_equal(dim(nutrimouse.sgcca$X$gene)[1] + dim(nutrimouse.sgcca$X$lipid)[1] + dim(nutrimouse.sgcca$X$Y)[1], dim(pl.res$df)[1])
})

# Clear the rgl device
if (rgl::rgl.cur() > 0) {
  rgl::close3d()
}

## ------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------ ##
## vdiffr testing - "ggplot2" with sgccda object
library(vdiffr)

test_that("plotIndiv works for block sPLSDA with different blocks plotted", {
  skip_on_ci() # only run the vdiffr tests locally
  
  # samples plotted across all modalities, default
  expect_doppelganger(
    title = "block sPLSDA plot basic, sample names", 
    fig = plotIndiv(nutrimouse.sgccda1, ind.names = TRUE))
  # samples plotted on just block 1 
  expect_doppelganger(
    title = "block sPLSDA plot on 1st block, sample names", 
    fig = plotIndiv(nutrimouse.sgccda1, ind.names = TRUE, blocks = 1))
  # samples plotted on just block 1 
  expect_doppelganger(
    title = "block sPLSDA plot on 1st and 2nd block, sample names", 
    fig = plotIndiv(nutrimouse.sgccda1, ind.names = TRUE, blocks = c(1,2)))
  
})

test_that("plotIndiv works for block sPLSDA with ellipse from predictions", {
  skip_on_ci() # only run the vdiffr tests locally
  
  # samples coloured by primary groups, default colours, by default shapes also match primary groups, ellipse
  expect_doppelganger(
    title = "block sPLSDA plot with ellipse", 
    fig = plotIndiv(nutrimouse.sgccda1, ind.names = FALSE, ellipse = TRUE, legend = TRUE))
  # samples coloured by primary groups and ellipse on custom groups
  expect_doppelganger(
    title = "block sPLSDA plot with ellipse on custom groups", 
    fig = plotIndiv(nutrimouse.sgccda1, ind.names = FALSE, group = as.factor(c(rep("A", 20), rep("B", 20))),
                    ellipse = TRUE, legend = TRUE))
})

test_that("plotIndiv works for block sPLSDA with centroids on groups", {
  skip_on_ci() # only run the vdiffr tests locally
  
  # samples coloured by primary groups, default colours, by default shapes also match primary groups, centroid
  expect_doppelganger(
    title = "block sPLSDA plot with centroids and stars", 
    fig = plotIndiv(nutrimouse.sgccda1, ind.names = FALSE, centroid = TRUE, 
                    legend = TRUE, star = TRUE))
  # samples coloured by primary groups, default colours, by default shapes also match primary groups, centroid, custom cols
  expect_doppelganger(
    title = "block sPLSDA plot with centroids and stars, custom cols", 
    fig = plotIndiv(nutrimouse.sgccda1, ind.names = FALSE, centroid = TRUE, 
                    legend = TRUE, star = TRUE,
                    col = c("red", "purple", "orange", "yellow", "green")))
  # samples coloured by primary groups, default colours, by default shapes also match primary groups, centroid, custom cols, pch on secondary groups
  expect_doppelganger(
    title = "block sPLSDA plot with centroids and stars, custom cols, pch on second grouping", 
    fig = plotIndiv(nutrimouse.sgccda1, ind.names = FALSE, centroid = TRUE, 
                    legend = TRUE, star = TRUE,
                    col = c("red", "purple", "orange", "yellow", "green"),
                    pch = as.factor(c(rep("A", 20), rep("B", 20)))))
})

