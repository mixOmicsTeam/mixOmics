context("plotLoadings.pls")

# Create a sPLS object for testing
data(liver.toxicity)
X = liver.toxicity$gene
Y = liver.toxicity$clinic
spls.obj = spls(X, Y, ncomp = 2, keepX = c(10, 10), keepY = c(5, 5))

# does not run due to 'plotLoadings encountered margin errors. Ensure feature names are not too long and the "Plots" pane is enlarged.'
# # Unit test 1: Test default behavior
# test_that("Test default behavior with graphics style", {
#   result <- plotLoadings(spls.obj, comp = 1, style = "graphics")
#   expect_equal(class(result), "list")
#   expect_equal(length(result), 2)  # Expect 2 elements in list
# })

# Unit test 2: Test incorrect 'col' value
test_that("Test invalid 'col' argument", {
  expect_error(plotLoadings(spls.obj, comp = 1, col = "nonexistent_color"),
               "'col' must be a single valid color.")
})

# Unit test 3: Test invalid 'comp' value
test_that("Test invalid 'comp' argument", {
  expect_error(plotLoadings.pca(spls.obj, comp = -1, style = "graphics"),
               "'comp' must be a positive integer.")
})

# Unit test 4: Test 'ndisplay' argument
test_that("Test invalid 'ndisplay' argument", {
  expect_error(plotLoadings(spls.obj, comp = 1, ndisplay = -5, style = "graphics"),
               "'ndisplay' must be a positive integer.")
})

## ------------------------------------------------------------------------ ##
## vdiffr testing - "graphics"
library(vdiffr)

test_that("plotLoadings works for spls graphics", {
  skip_on_ci() # only run the vdiffr tests locally
  
  # simple plot
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot spls - graphics", 
      fig = plotLoadings(spls.obj))
  ))
  # only top 3 genes, change gene names
  new_names <- list(c(paste0("Gene_", 1:3116)), c(paste0("Clinical_", 1:10)))
  names(new_names) <- c("X", "Y")
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot spls change gene names and plot top 3  - graphics",
      fig = plotLoadings.mixo_pls(spls.obj, style = "graphics", ndisplay = 3, size.name = 1.5, name.var = new_names))
  ))
  # change labels and label sizes
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot spls change labels and label sizes - graphics",
      fig = plotLoadings(spls.obj, X.label = "X", Y.label = "Y", title = "Title Test",
                         size.labs = 2, size.title = 3, size.subtitle = 1, size.axis = 1))
  ))
  # change colours and borders
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot spls change cols and borders - graphics",
      fig = plotLoadings(spls.obj, col = "red", border = "grey", xlim = c(-1, 1)))
  ))
  # change layout
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot spls change layout - graphics",
      fig = plotLoadings(spls.obj, layout = c(2,1)))
  ))
  # plot only one block
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot spls change block - graphics",
      fig = plotLoadings(spls.obj, block = "Y"))
  ))
  
})

## ------------------------------------------------------------------------ ##
## vdiffr testing - "ggplot2"

test_that("plotLoadings works for pls ggplot2", {
  skip_on_ci() # only run the vdiffr tests locally
  
  # simple plot
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot spls simple - ggplot2", 
      fig = plotLoadings(spls.obj, style = "ggplot2"))
  ))
  # only top 3 genes, change gene names
  new_names <- list(c(paste0("Gene_", 1:3116)), c(paste0("Clinical_", 1:10)))
  names(new_names) <- c("X", "Y")
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot spls change gene names and plot top 3  - ggplot2",
      fig = plotLoadings.mixo_pls(spls.obj, style = "ggplot2", ndisplay = 3, size.name = 1.5, name.var = new_names))
  ))
  # change labels and label sizes
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot spls change labels and label sizes - ggplot2",
      fig = plotLoadings(spls.obj, style = "ggplot2", X.label = "X", Y.label = "Y", title = "Title Test",
                         size.labs = 2, size.title = 3, size.axis = 1))
  ))
  # change colours and borders
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot spls change cols and borders - ggplot2",
      fig = plotLoadings(spls.obj, style = "ggplot2", col = "green", border = "grey", xlim = c(-1, 1)))
  ))
  # change layout
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot spls change layout - ggplot2",
      fig = plotLoadings(spls.obj, layout = c(2,1), style = "ggplot2"))
  ))
  # plot only one block
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot spls change block - ggplot2",
      fig = plotLoadings(spls.obj, block = "Y", style = "ggplot2"))
  ))
  
})
