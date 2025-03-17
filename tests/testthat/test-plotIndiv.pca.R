context("plotIndiv.pca")


## ------------------------------------------------------------------------ ##
## Test that outputs are correct when running default style = "ggplot2"

test_that("plotIndiv works for (s)pca", {
  data(srbct)
  X <- srbct$gene
  Y <- srbct$class
  pca.srbct = pca(X, ncomp = 10, center = TRUE, scale = TRUE)
  groups <- factor(srbct$class, levels = c("RMS", "NB", "EWS", "BL"))
  pl.res <- plotIndiv(pca.srbct, group = groups, ind.names = FALSE, # plot the samples projected
            legend = TRUE, title = 'PCA on SRBCT, comp 1 - 2',
            col = c("red", "blue", "green", "black")) # onto the PCA subspace
  
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
