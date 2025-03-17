context("plotIndiv.mint")

## ------------------------------------------------------------------------ ##
## Test that outputs are correct when running default style = "ggplot2"

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


## ------------------------------------------------------------------------ ##
## Plotting with 'graphics' style


## ------------------------------------------------------------------------ ##
## Plotting with '3d' style

library(rgl)

# Clear the rgl device
if (rgl::rgl.cur() > 0) {
  rgl::close3d()
}
