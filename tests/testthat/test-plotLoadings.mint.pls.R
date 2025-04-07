context("plotLoadings.mint.pls")

# Create a MINT PLS object for testing
data(stemcells)
studies <- stemcells$study
levels(studies) <- c("Study_A", "Study_B", "Study_C", "Study_D")
mint.spls.obj <- mint.pls(X = stemcells$gene[,20:25], Y = stemcells$gene[,1:10], ncomp = 3, study = studies)

# Unit test 1: Test default behavior
test_that("Test default behavior with graphics style", {
  result <- plotLoadings(mint.spls.obj, comp = 1, style = "graphics")
  expect_equal(class(result), "list")
  expect_true(all(sapply(result, is.data.frame)))
})

# Unit test 2: Test incorrect 'col' value
test_that("Test invalid 'col' argument", {
  expect_error(plotLoadings(mint.spls.obj, comp = 1, col = "nonexistent_color"),
               "'col' must be a single valid color.")
})

# Unit test 3: Test invalid 'comp' value
test_that("Test invalid 'comp' argument", {
  expect_error(plotLoadings(mint.spls.obj, comp = 0),
               "'comp' must be a positive integer.")
})

# Unit test 4: Test 'ndisplay' argument
test_that("Test invalid 'ndisplay' argument", {
  expect_error(plotLoadings(mint.spls.obj, comp = 1, ndisplay = -5),
               "'ndisplay' must be a positive integer.")
})

# Unit test 5: Test invalid 'study' argument
test_that("Test invalid 'study' argument", {
  expect_error(plotLoadings(mint.spls.obj, comp = 1, study = "Invalid_Study"),
               "study' must from one of 'object\\$study', 'global' or 'all.partial', see help file.")
})

# Unit test 6: Test invalid 'block' argument with all.partial
test_that("Test invalid 'block' argument with all.partial", {
  expect_error(plotLoadings(mint.spls.obj, comp = 1, study = "all.partial", block = c("X", "Y")),
               "When study = 'all.partial' or a specific study is specified, only one block can be plotted at a time")
})

## ------------------------------------------------------------------------ ##
## vdiffr testing - "graphics"
library(vdiffr)

test_that("plotLoadings works for mint.pls graphics", {
  skip_on_ci() # only run the vdiffr tests locally
  
  # simple plot
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot mint.pls - graphics", 
      fig = plotLoadings(mint.spls.obj))
  ))
  
  # plot with all.partial
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot mint.pls all.partial - graphics",
      fig = plotLoadings(mint.spls.obj, study = "all.partial"))
  ))
  
  # plot specific study
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot mint.pls specific study - graphics",
      fig = plotLoadings(mint.spls.obj, study = "Study_A"))
  ))
  
  # change gene names
  new_names <- list(c(paste0("Gene_", 1:1:ncol(mint.spls.obj$X))), c(paste0("Clinical_", 1:1:ncol(mint.spls.obj$Y))))
  names(new_names) <- c("X", "Y")
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot mint.pls change gene names - graphics",
      fig = plotLoadings(mint.spls.obj, name.var = new_names))
  ))
  
  # change labels and label sizes
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot mint.pls change labels and label sizes - graphics",
      fig = plotLoadings(mint.spls.obj, X.label = "X", Y.label = "Y", title = "Title Test",
                         size.labs = 2, size.title = 3, size.subtitle = 1, size.axis = 1))
  ))
  
  # change colours and borders
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot mint.pls change cols and borders - graphics",
      fig = plotLoadings(mint.spls.obj, col = "red", border = "grey", xlim = c(-1, 1)))
  ))
  
  # change layout
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot mint.pls change layout - graphics",
      fig = plotLoadings(mint.spls.obj, study = c("Study_A", "Study_B"), layout = c(2,1)))
  ))
})

## ------------------------------------------------------------------------ ##
## vdiffr testing - "ggplot2"

test_that("plotLoadings works for mint.pls ggplot2", {
  skip_on_ci() # only run the vdiffr tests locally

  # simple plot
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot mint.pls simple - ggplot2",
      fig = plotLoadings(mint.spls.obj, style = "ggplot2"))
  ))

  # plot with all.partial
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot mint.pls all.partial - ggplot2",
      fig = plotLoadings(mint.spls.obj, study = "all.partial", style = "ggplot2"))
  ))

  # plot specific study
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot mint.pls specific study - ggplot2",
      fig = plotLoadings(mint.spls.obj, study = "Study_A", style = "ggplot2"))
  ))

  # change gene names
  new_names <- list(c(paste0("Gene_", 1:1:ncol(mint.spls.obj$X))), c(paste0("Clinical_", 1:1:ncol(mint.spls.obj$Y))))
  names(new_names) <- c("X", "Y")
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot mint.pls change gene names - ggplot2",
      fig = plotLoadings(mint.spls.obj, name.var = new_names, style = "ggplot2"))
  ))

  # change labels and label sizes
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot mint.pls change labels and label sizes - ggplot2",
      fig = plotLoadings(mint.spls.obj, style = "ggplot2", X.label = "X", Y.label = "Y", title = "Title Test",
                         size.labs = 2, size.title = 3, size.axis = 1))
  ))

  # change colours and borders
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot mint.pls change cols and borders - ggplot2",
      fig = plotLoadings(mint.spls.obj, style = "ggplot2", col = "green", border = "grey", xlim = c(-1, 1)))
  ))

  # change layout
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot mint.pls change layout - ggplot2",
      fig = plotLoadings(mint.spls.obj, study = c("Study_A", "Study_B"), layout = c(2,1), style = "ggplot2"))
  ))
})