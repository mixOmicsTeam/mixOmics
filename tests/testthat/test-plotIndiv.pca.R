context("plotIndiv.pca")
library(vdiffr)

## ------------------------------------------------------------------------ ##

## Full data object and grouping
data("srbct")
pca.obj.full <- pca(srbct$gene, ncomp = 3)
groups <- srbct$class

## Smaller data object and grouping
X <- srbct$gene[1:6, ]
rownames(X) <- c(paste0("Sample_", 1:6))
pca.obj <- pca(X, ncomp = 3)

primary_groups <- as.factor(c(rep("Group_1", 2), rep("Group_2", 2), rep("Group_3", 2)))
# [1] Group_1 Group_1 Group_2 Group_2 Group_3 Group_3
# Levels: Group_1 Group_2 Group_3
secondary_groups <- as.factor(c(rep("A", 3), rep("B", 2), rep("C", 1)))
# [1] A A A B B C
# Levels: A B C

## ------------------------------------------------------------------------ ##
## Test that outputs are correct when running default style = "ggplot2"

test_that("plotIndiv works for pca with sample names (default)", {
  # simple plot showing sample names
  pl.res <- expect_doppelganger(
    title = "Basic PCA plot", 
    fig = plotIndiv(pca.obj))
  # samples coloured by primary groups, sample names, default colours
  pl.res <- expect_doppelganger(
    title = "Basic PCA plot coloured by primary groups",
    fig = plotIndiv(pca.obj, group = primary_groups))
  # samples coloured by primary groups, sample names, user-defined colours
  pl.res <- expect_doppelganger(
    title = "Basic PCA plot coloured by primary groups custom cols",
    fig = plotIndiv(pca.obj, group = primary_groups, col = c("red", "blue", "black")))
  
})

test_that("plotIndiv works for pca without sample names", {
  # simple plot
  pl.res <- plotIndiv(pca.obj, ind.names = FALSE)
  # samples coloured by primary groups, default colours, by default shapes also match primary groups
  pl.res <- plotIndiv(pca.obj, ind.names = FALSE,
                      group = primary_groups, legend = TRUE)
  # samples coloured by primary groups, user-defined colours, by default shapes also match primary groups
  pl.res <- plotIndiv(pca.obj, ind.names = FALSE,
                      group = primary_groups, col = c("red", "blue", "black"), legend = TRUE)
  # samples coloured by primary groups, user-defined colours, groups reordered, by default shapes also match primary groups
  pl.res <- plotIndiv(pca.obj, ind.names = FALSE,
                      group = factor(primary_groups, levels = c("Group_2", "Group_3", "Group_1")),
                      col = c("red", "blue", "black"), legend = TRUE)

})

test_that("plotIndiv works for pca with ellipse on groups", {
  # use full data so can make ellipse

  # samples coloured by primary groups, default colours, by default shapes also match primary groups, ellipse
  pl.res <- plotIndiv(pca.obj.full, ind.names = FALSE,
                      group = groups, ellipse = TRUE, legend = TRUE)
  # samples coloured by primary groups, user-defined colours, by default shapes also match primary groups, ellipse
  pl.res <- plotIndiv(pca.obj.full, ind.names = FALSE,
                      group = groups, ellipse = TRUE, col = c("red", "blue", "black", "yellow"), legend = TRUE)
  # samples coloured by primary groups, user-defined colours, by default shapes also match primary groups, ellipse, sample names
  pl.res <- plotIndiv(pca.obj.full, ind.names = TRUE,
                      group = groups, ellipse = TRUE, col = c("red", "blue", "black", "yellow"), legend = TRUE)
  # samples coloured by primary groups, default colours, by default shapes also match primary groups, ellipse confidence
  pl.res <- plotIndiv(pca.obj.full, ind.names = FALSE, ellipse = TRUE, ellipse.level = 0.5,
                      group = groups, legend = TRUE)


})


test_that("plotIndiv works for pca with centroids on groups", {
  # use full data so can make centroids

  # samples coloured by primary groups, default colours, by default shapes also match primary groups, ellipse
  pl.res <- plotIndiv(pca.obj.full, ind.names = FALSE,
                      group = groups, centroid = TRUE, legend = TRUE)
  # samples coloured by primary groups, user-defined colours, by default shapes also match primary groups, ellipse
  pl.res <- plotIndiv(pca.obj.full, ind.names = FALSE,
                      group = groups, centroid = TRUE, col = c("red", "blue", "black", "yellow"), legend = TRUE)
  # samples coloured by primary groups, user-defined colours, by default shapes also match primary groups, ellipse, sample names
  pl.res <- plotIndiv(pca.obj.full, ind.names = TRUE,
                      group = groups, centroid = TRUE, col = c("red", "blue", "black", "yellow"), legend = TRUE)
  # samples coloured by primary groups, default colours, by default shapes also match primary groups, ellipse confidence
  pl.res <- plotIndiv(pca.obj.full, ind.names = FALSE, centroid = TRUE, star = TRUE,
                      group = groups, legend = TRUE)


})

test_that("plotIndiv controlling colours and pch", {
  # force pch to be the same for all samples, so groups only coloured differently
  plotIndiv(pca.obj, group = primary_groups, legend = TRUE,
            pch = 1)
  # force pch to be the same for all samples, so groups only coloured differently - different shape
  plotIndiv(pca.obj, group = primary_groups, legend = TRUE,
            pch = 3)
  # control the pch of each of the primary groups by setting pch with length = levels(primary_groups)
  plotIndiv(pca.obj, group = primary_groups, legend = TRUE,
            pch = c(1, 2, 3))
  # use pch to show a whole new grouping (secondary groups)
  plotIndiv(pca.obj, group = primary_groups, legend = TRUE,
            pch = secondary_groups)
  # change order of secondary grouping
  plotIndiv(pca.obj, group = primary_groups, legend = TRUE,
            pch = factor(secondary_groups, levels = c("B", "C", "A")))

  # !!!force specific pch for specific secondary groupings - this doesn't work!
  plotIndiv(pca.obj, group = primary_groups, legend = TRUE,
            pch = c("10", "10", "10", "5", "5", "6"))


})



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

## ------------------------------------------------------------------------ ##
## Plotting with 'lattice' style

test_that("plotIndiv works for (s)pca (lattice style)", {
  data(srbct)
  X <- srbct$gene
  Y <- srbct$class
  pca.srbct = pca(X, ncomp = 10, center = TRUE, scale = TRUE)
  groups <- factor(srbct$class, levels = c("RMS", "NB", "EWS", "BL"))
  pl.res <- plotIndiv(pca.srbct, group = groups, ind.names = FALSE, # plot the samples projected
                      legend = TRUE, title = 'PCA on SRBCT, comp 1 - 2',
                      col = c("red", "blue", "green", "black"),
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

## ------------------------------------------------------------------------ ##
## Plotting with 'graphics' style

test_that("plotIndiv works for (s)pca (graphics style)", {
  data(srbct)
  X <- srbct$gene
  Y <- srbct$class
  pca.srbct = pca(X, ncomp = 10, center = TRUE, scale = TRUE)
  groups <- factor(srbct$class, levels = c("RMS", "NB", "EWS", "BL"))
  pl.res <- plotIndiv(pca.srbct, group = groups, ind.names = FALSE, # plot the samples projected
                      legend = TRUE, title = 'PCA on SRBCT, comp 1 - 2',
                      col = c("red", "blue", "green", "black"),
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

## ------------------------------------------------------------------------ ##
## Plotting with '3d' style

library(rgl)

test_that("plotIndiv works for (s)pca (3d style)", {
  data(srbct)
  X <- srbct$gene
  Y <- srbct$class
  pca.srbct = pca(X, ncomp = 10, center = TRUE, scale = TRUE)
  groups <- factor(srbct$class, levels = c("RMS", "NB", "EWS", "BL"))

  clear3d()
  pl.res <- suppressWarnings(suppressMessages(plotIndiv(pca.srbct, group = groups, ind.names = FALSE,
                                                        legend = TRUE, title = 'PCA on SRBCT, comp 1 - 2',
                                                        col = c("red", "blue", "green", "black"),
                                                        style = "3d", pch = "sphere")))

  # Check coordinates
  .expect_numerically_close(pl.res$df[1,1], 10.13857)
  # Check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # Check colour assignments are correct
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "RMS"]), "red")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "NB"]), "blue")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "EWS"]), "green")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "BL"]), "black")
  # Check right number of samples
  expect_equal(dim(pca.srbct$X)[1], dim(pl.res$df)[1])
  # Check groups
  expect_true(!is.null(pl.res$df$group))
  expect_equal(length(unique(pl.res$df$group)), length(unique(groups)))
})
