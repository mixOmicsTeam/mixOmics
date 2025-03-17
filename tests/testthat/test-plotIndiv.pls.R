context("plotIndiv.pls")


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

test_that("plotIndiv works for (s)pls", {
  data(liver.toxicity)
  X <- liver.toxicity$gene
  Y <- liver.toxicity$clinic
  toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50),
                        keepY = c(10, 10, 10))
  
  pl.res <- plotIndiv(toxicity.spls, rep.space="X-variate", ind.name = FALSE,
                      group = factor(liver.toxicity$treatment$Time.Group),
                      legend.title = 'Time',
                      col = c("red", "blue", "green", "black"),
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
  pl.res <- plotIndiv(srbct.splsda , comp = 1:2, col = c("red", "blue", "green", "black"),
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

## ------------------------------------------------------------------------ ##

test_that("plotIndiv.rcc works without ind.names", code = {
    data(nutrimouse)
    X <- nutrimouse$lipid
    Y <- nutrimouse$gene
    nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
    
    plotIndiv.res <- plotIndiv(nutri.res, group = nutrimouse$genotype, ind.names = FALSE, legend = TRUE)
    
    expect_is(plotIndiv.res$graph, "ggplot")
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

test_that("plotIndiv works for (s)pls (lattice style)", {
  data(liver.toxicity)
  X <- liver.toxicity$gene
  Y <- liver.toxicity$clinic
  toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50),
                        keepY = c(10, 10, 10))
  
  pl.res <- plotIndiv(toxicity.spls, rep.space="X-variate", ind.name = FALSE,
                      group = factor(liver.toxicity$treatment$Time.Group),
                      legend.title = 'Time',
                      col = c("red", "blue", "green", "black"),
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

test_that("plotIndiv works for (s)pls (graphics style)", {
  data(liver.toxicity)
  X <- liver.toxicity$gene
  Y <- liver.toxicity$clinic
  toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50),
                        keepY = c(10, 10, 10))
  
  pl.res <- plotIndiv(toxicity.spls, rep.space="X-variate", ind.name = FALSE,
                      group = factor(liver.toxicity$treatment$Time.Group),
                      legend.title = 'Time',
                      col = c("red", "blue", "green", "black"),
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

library(rgl)

test_that("plotIndiv works for rcc (3d style)", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
  
  # Clear any existing rgl plots and disable plot rendering for non-interactive environments
  options(rgl.useNULL = TRUE)
  clear3d()
  pl.res <- suppressWarnings(suppressMessages(plotIndiv(nutri.res, style = "3d", pch = "sphere")))
  
  # Check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # Check coordinates
  .expect_numerically_close(pl.res$df[1,1], 0.87088852)
  
  clear3d()
  pl.res <- suppressWarnings(suppressMessages(plotIndiv(nutri.res, rep.space = 'XY-variate', group = nutrimouse$genotype,
                                                        legend = TRUE, style = "3d", pch = "sphere")))
  
  # Check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # Check coordinates
  .expect_numerically_close(pl.res$df[1,1], 0.8270997)
  # Check groups
  expect_true(!is.null(pl.res$df$group))
  expect_equal(length(unique(pl.res$df$group)), length(unique(nutrimouse$genotype)))
})

test_that("plotIndiv works for (s)pls (3d style)", {
  data(liver.toxicity)
  X <- liver.toxicity$gene
  Y <- liver.toxicity$clinic
  toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50),
                        keepY = c(10, 10, 10))
  
  # Expect error when passing numbers as pch argument
  expect_error(suppressWarnings(suppressMessages(plotIndiv(toxicity.spls, rep.space="X-variate", ind.name = FALSE,
                                                           group = factor(liver.toxicity$treatment$Time.Group),
                                                           legend.title = 'Time',
                                                           col = c("red", "blue", "green", "black"),
                                                           pch = factor(liver.toxicity$treatment$Dose.Group),
                                                           legend.title.pch = 'Dose',
                                                           legend = TRUE, style = "3d"))),
               "pch' must be a simple character or character vector from {'sphere', 'tetra', 'cube', 'octa', 'icosa', 'dodeca'}.",
               fixed = TRUE)
  
  # Plot with correct pch values
  pchs <- factor(liver.toxicity$treatment$Dose.Group)
  levels(pchs) <- c("sphere", "tetra", "octa", "icosa")
  
  clear3d()
  pl.res <- suppressWarnings(suppressMessages(plotIndiv(toxicity.spls, rep.space="X-variate", ind.name = FALSE,
                                                        group = factor(liver.toxicity$treatment$Time.Group),
                                                        legend.title = 'Time',
                                                        col = c("red", "blue", "green", "black"),
                                                        pch = pchs,
                                                        legend.title.pch = 'Dose',
                                                        legend = TRUE, style = "3d")))
  
  # Check coordinates
  .expect_numerically_close(pl.res$df[1,1], 4.146771)
  # Check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # Check colour assignments are correct
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 6]), "red")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 18]), "blue")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 24]), "green")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 48]), "black")
  # Check right number of samples
  expect_equal(dim(toxicity.spls$X)[1], dim(pl.res$df)[1])
  # Check groups
  expect_true(!is.null(pl.res$df$group))
  expect_equal(length(unique(pl.res$df$group)), length(unique(liver.toxicity$treatment$Time.Group)))
})

test_that("plotIndiv works for (s)plsda (3d style)", {
  data(breast.tumors)
  X <- breast.tumors$gene.exp
  Y <- breast.tumors$sample$treatment
  splsda.breast <- splsda(X, Y, keepX=c(10,10), ncomp=3)
  
  clear3d()
  pl.res <- suppressWarnings(suppressMessages(plotIndiv(splsda.breast, style = "3d", pch = "cube")))
  
  # Check coordinates
  .expect_numerically_close(pl.res$df[1,1], -1.075222)
  # Check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # Check colours
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "BE"]), "#F68B33")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "AF"]), "#388ECC")
  # Check right number of samples
  expect_equal(dim(splsda.breast$X)[1], dim(pl.res$df)[1])
})

# Clear the rgl device
if (rgl::rgl.cur() > 0) {
  rgl::close3d()
}
