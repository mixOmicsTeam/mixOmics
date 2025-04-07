context("plotLoadings.pca")

# Create a PCA object for testing
data(srbct)
X <- srbct$gene[1:6, 1:10]
rownames(X) <- c(paste0("Sample_", 1:6))
pca.obj <- pca(X, ncomp = 3)

# Unit test 1: Test default behavior
test_that("Test default behavior with graphics style", {
  png(tempfile(), width = 1200, height = 1000, res = 150)
  result <- plotLoadings(pca.obj, comp = 1, style = "graphics")
  dev.off()
  expect_s3_class(result, "data.frame")
  expect_equal(ncol(result), 2)  # Expect 2 columns in output df
})

# Unit test 3: Test incorrect 'col' value
test_that("Test invalid 'col' argument", {
  expect_error(plotLoadings(pca.obj, comp = 1, col = "nonexistent_color"),
               "'col' must be a single valid color.")
})

# Unit test 4: Test invalid 'comp' value
test_that("Test invalid 'comp' argument", {
  expect_error(plotLoadings.pca(pca.obj, comp = -1, style = "graphics"),
               "'comp' must be a positive integer.")
})

# Unit test 5: Test 'ndisplay' argument
test_that("Test invalid 'ndisplay' argument", {
  expect_error(plotLoadings(pca.obj, comp = 1, ndisplay = -5, style = "graphics"),
               "'ndisplay' must be a positive integer.")
})

## ------------------------------------------------------------------------ ##
## vdiffr testing - "graphics"
library(vdiffr)

test_that("plotLoadings works for pca graphics", {
  skip_on_ci() # only run the vdiffr tests locally
  
  # simple plot
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot simple - graphics", 
      fig = plotLoadings(pca.obj))
  ))
  # only top 3 genes, change gene names
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot change gene names and plot top 3  - graphics",
      fig = plotLoadings(pca.obj, ndisplay = 3, size.name = 1.5, name.var = c(paste0("Gene_", 1:10))))
  ))
  # change labels and label sizes
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot change labels and label sizes - graphics",
      fig = plotLoadings(pca.obj, X.label = "X", Y.label = "Y", title = "Title Test",
                         size.labs = 2, size.title = 3, size.axis = 1))
  ))
  # change colours and borders
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot change cols and borders - graphics",
      fig = plotLoadings(pca.obj, col = "red", border = "grey", xlim = c(-1, 1)))
  ))
  
})

## ------------------------------------------------------------------------ ##
## vdiffr testing - "ggplot2"

test_that("plotLoadings works for pca ggplot2", {
  skip_on_ci() # only run the vdiffr tests locally
  
  # simple plot
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot simple - ggplot2", 
      fig = plotLoadings(pca.obj, style = "ggplot2"))
  ))
  # only top 3 genes, change gene names
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot change gene names and plot top 3  - ggplot2",
      fig = plotLoadings(pca.obj, style = "ggplot2", ndisplay = 3, size.name = 1.5, name.var = c(paste0("Gene_", 1:10))))
  ))
  # change labels and label sizes
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot change labels and label sizes - ggplot2",
      fig = plotLoadings(pca.obj, style = "ggplot2", X.label = "X", Y.label = "Y", title = "Title Test",
                         size.labs = 2, size.title = 3, size.axis = 1))
  ))
  # change colours and borders
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot change cols and borders - ggplot2",
      fig = plotLoadings(pca.obj, style = "ggplot2", col = "red", border = "grey", xlim = c(-1, 1)))
  ))
  
})
