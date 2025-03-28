context("plotIndiv.mint")

#set up model
data(stemcells)
mint.res <- mint.splsda(X = stemcells$gene, Y = stemcells$celltype, ncomp = 2, keepX = c(10, 5),
                  study = stemcells$study)

## ------------------------------------------------------------------------ ##
## Test that outputs are correct when running default style = "ggplot2"

test_that("plotIndiv works for mint.(s)plsda", {
  pl.res <- plotIndiv(mint.res)
  # check coordinates
  .expect_numerically_close(pl.res$graph$data$x[1], -1.543685)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "graph"))
  # check right number of samples
  expect_equal(nrow(mint.res$X), nrow(pl.res$df))
})

## ------------------------------------------------------------------------ ##
## Edge cases

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

test_that("plotIndiv works for mint.(s)plsda - lattice", {
  pl.res <-   plotIndiv(mint.res, style = "lattice")
  # check coordinates
  .expect_numerically_close(pl.res$df$x[1], -1.543685)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "graph"))
  # check right number of samples
  expect_equal(dim(mint.res$X)[1], dim(pl.res$df)[1])
})

## ------------------------------------------------------------------------ ##
## Plotting with 'graphics' style

test_that("plotIndiv works for mint.(s)plsda - graphics", {
  pl.res <-   plotIndiv(mint.res, style = "graphics")
  # check coordinates
  .expect_numerically_close(pl.res$df$x[1], -1.543685)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "graph"))
  # check right number of samples
  expect_equal(dim(mint.res$X)[1], dim(pl.res$df)[1])
})

## ------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------ ##
## vdiffr testing - "ggplot2" with MINT sPLS-DA object
library(vdiffr)

test_that("plotIndiv works for mint sPLSDA plotting studies separately and together", {
  skip_on_ci() # only run the vdiffr tests locally

  invisible(capture.output(
  # all studies in one plot
  expect_doppelganger(
    title = "mint sPLSDA plot global",
    fig = plotIndiv(mint.res, study = "global"))
  ))
  invisible(capture.output(
  # studies in separate plots
  expect_doppelganger(
    title = "mint sPLSDA plot studies facetted",
    fig = plotIndiv(mint.res, study = "all.partial"))
  ))
  invisible(capture.output(
  # studies in separate plots - different layout
  expect_doppelganger(
    title = "mint sPLSDA plot studies facetted different layout",
    fig = plotIndiv(mint.res, study = "all.partial", layout = c(3,4)))
  ))
})

test_that("plotIndiv works for mint sPLSDA plotting controlled groups cols, ellipse, etc", {
  skip_on_ci() # only run the vdiffr tests locally

  invisible(capture.output(
  # default groups custom cols
  expect_doppelganger(
    title = "mint sPLSDA plot custom cols default groups",
    fig = plotIndiv(mint.res, col = c("red", "yellow", "blue")))
  ))
  invisible(capture.output(
  # default cols custom groups
  expect_doppelganger(
    title = "mint sPLSDA plot default cols custom groups",
    fig = plotIndiv(mint.res, group = as.factor(c(rep("A", 50), rep("B", 50), rep("C", 25))),
    legend = TRUE))
  ))
  invisible(capture.output(
  # ellipse
  expect_doppelganger(
    title = "mint sPLSDA plot with ellipse",
    fig = plotIndiv(mint.res, ellipse = TRUE, legend = TRUE))
  ))
  invisible(capture.output(
  # star
  expect_doppelganger(
    title = "mint sPLSDA plot with star",
    fig = plotIndiv(mint.res, centroid = TRUE, star = TRUE, legend = TRUE))
  ))

})

test_that("plotIndiv works for mint sPLSDA with background", {
  skip_on_ci() # only run the vdiffr tests locally

  background = background.predict(mint.res, comp.predicted=2, dist = "max.dist")
  invisible(capture.output(
  expect_doppelganger(
    title = "mint sPLSDA plot with background",
    fig = plotIndiv(mint.res, background = background, legend = TRUE))
  ))
})

## vdiffr testing - "ggplot2" with MINT sPLS object
res = mint.pls(X = stemcells$gene[,-1], Y = stemcells$gene[,1], ncomp = 3,
               study = stemcells$study)

test_that("plotIndiv works for mint PLS plotting studies separately and together", {
  skip_on_ci() # only run the vdiffr tests locally
  
  invisible(capture.output(
    # all studies in one plot
    expect_doppelganger(
      title = "mint sPLS plot global",
      fig = plotIndiv(res, study = "global"))
  ))
  invisible(capture.output(
    # studies in separate plots
    expect_doppelganger(
      title = "mint PLS plot studies facetted",
      fig = suppressWarnings(plotIndiv(res, study = "all.partial")))
  ))
  invisible(capture.output(
    # studies in separate plots - different layout
    expect_doppelganger(
      title = "mint sPLS plot studies facetted different layout",
      fig = suppressWarnings(plotIndiv(res, study = "all.partial", layout = c(3,4))))
  ))
})

test_that("plotIndiv works for mint PLS plotting controlled groups cols, ellipse, etc", {
  skip_on_ci() # only run the vdiffr tests locally
  
  invisible(capture.output(
    # default groups custom cols
    expect_doppelganger(
      title = "mint PLS plot custom cols default groups",
      fig = plotIndiv(res, col = c("red", "yellow", "blue")))
  ))
  invisible(capture.output(
    # default cols custom groups
    expect_doppelganger(
      title = "mint PLS plot default cols custom groups",
      fig = plotIndiv(res, group = as.factor(c(rep("A", 50), rep("B", 50), rep("C", 25))),
                      legend = TRUE))
  ))
  
})