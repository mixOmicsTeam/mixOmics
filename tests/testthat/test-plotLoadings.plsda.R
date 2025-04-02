context("plotLoadings.plsda")

# Create a sPLS-DA object for testing
data(liver.toxicity)
X = as.matrix(liver.toxicity$gene)
Y = as.factor(paste0('treatment_', liver.toxicity$treatment[, 4]))
splsda.obj = splsda(X, Y, ncomp = 2, keepX = c(20, 20))

# Unit test 1: Test default behavior
test_that("Test default behavior with graphics style", {
  result <- plotLoadings(splsda.obj, comp = 1, method = 'median')
  expect_equal(class(result), "list")
  expect_equal(length(result), 1)  # Expect 1 element in list for sPLS-DA
})

# Unit test 2: Test contrib argument
test_that("Test contrib argument", {
  # Test with contrib = "max"
  result_max <- plotLoadings(splsda.obj, comp = 1, method = 'median', contrib = "max")
  expect_equal(class(result_max), "list")
  
  # Test with contrib = "min"
  result_min <- plotLoadings(splsda.obj, comp = 1, method = 'median', contrib = "min")
  expect_equal(class(result_min), "list")
  
  # Test with contrib = NULL (should switch to classical plotLoadings)
  result_null <- plotLoadings(splsda.obj, comp = 1, contrib = NULL)
  expect_equal(class(result_null), "list")
})

# Unit test 3: Test method argument
test_that("Test method argument", {
  # Test with method = "mean"
  result_mean <- plotLoadings(splsda.obj, comp = 1, method = 'mean', contrib = "max")
  expect_equal(class(result_mean), "list")
  
  # Test with method = "median"
  result_median <- plotLoadings(splsda.obj, comp = 1, method = 'median', contrib = "max")
  expect_equal(class(result_median), "list")
  
  # Test with invalid method (should default to median)
  expect_warning(plotLoadings(splsda.obj, comp = 1, method = 'invalid', contrib = "max"),
                "'method' should be either 'mean' or 'median', set to 'median' by default")
})

# Unit test 4: Test show.ties argument
test_that("Test show.ties argument", {
  # Test with show.ties = TRUE (default)
  result_show <- plotLoadings(splsda.obj, comp = 1, method = 'median', contrib = "max", show.ties = TRUE)
  expect_equal(class(result_show), "list")
  
  # Test with show.ties = FALSE
  result_hide <- plotLoadings(splsda.obj, comp = 1, method = 'median', contrib = "max", show.ties = FALSE)
  expect_equal(class(result_hide), "list")
})

# Unit test 5: Test legend arguments
test_that("Test legend arguments", {
  # Test with legend = TRUE (default)
  result_legend <- plotLoadings(splsda.obj, comp = 1, method = 'median', contrib = "max", legend = TRUE)
  expect_equal(class(result_legend), "list")
  
  # Test with legend = FALSE
  result_no_legend <- plotLoadings(splsda.obj, comp = 1, method = 'median', contrib = "max", legend = FALSE)
  expect_equal(class(result_no_legend), "list")
  
  # Test with custom legend colors
  result_custom_colors <- plotLoadings(splsda.obj, comp = 1, method = 'median', contrib = "max", 
                                     legend.color = c(1:4))
  expect_equal(class(result_custom_colors), "list")
})

# Unit test 6: Test name.var argument
test_that("Test name.var argument", {
  # Test with name.var
  name.var = liver.toxicity$gene.ID[, 'geneBank']
  names(name.var) = rownames(liver.toxicity$gene.ID)
  result_namevar <- plotLoadings(splsda.obj, comp = 1, method = 'median', contrib = "max", 
                                name.var = name.var)
  expect_equal(class(result_namevar), "list")
})

# Unit test 7: Test style argument
test_that("Test style argument", {
  # Test with style = "graphics" (default)
  result_graphics <- plotLoadings(splsda.obj, comp = 1, method = 'median', contrib = "max", style = "graphics")
  expect_equal(class(result_graphics), "list")
  
  # Test with style = "ggplot2"
  result_ggplot2 <- plotLoadings(splsda.obj, comp = 1, method = 'median', contrib = "max", style = "ggplot2")
  expect_equal(class(result_ggplot2), "list")
})

## ------------------------------------------------------------------------ ##
## vdiffr testing - "graphics"
library(vdiffr)

test_that("plotLoadings works for splsda graphics", {
  skip_on_ci() # only run the vdiffr tests locally
  
  # simple plot
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot splsda - graphics", 
      fig = plotLoadings(splsda.obj, contrib = "max"))
  ))
  
  # with custom names
  name.var = liver.toxicity$gene.ID[, 'geneBank']
  names(name.var) = rownames(liver.toxicity$gene.ID)
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot splsda with custom names - graphics",
      fig = plotLoadings(splsda.obj, contrib = "max", name.var = name.var))
  ))
  
  # with custom legend
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot splsda with custom legend - graphics",
      fig = plotLoadings(splsda.obj, contrib = "max", legend.title = "Treatment", 
                        legend.color = c(1:4)))
  ))
  
  # with custom title and labels
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot splsda with custom title and labels - graphics",
      fig = plotLoadings(splsda.obj, contrib = "max", title = "Liver Data", 
                        X.label = "Loading", Y.label = "Genes"))
  ))
})

## ------------------------------------------------------------------------ ##
## vdiffr testing - "ggplot2"

test_that("plotLoadings works for splsda ggplot2", {
  skip_on_ci() # only run the vdiffr tests locally
  
  # simple plot
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot splsda - ggplot2", 
      fig = plotLoadings(splsda.obj, contrib = "max", style = "ggplot2", ndisplay = 5))
  ))
  
  # with custom names
  new_names <- c(paste0("Gene_", 1:3116))
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot splsda with custom names - ggplot2",
      fig = plotLoadings(splsda.obj, contrib = "max", style = "ggplot2", name.var = new_names))
  ))
  
  # with custom legend
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot splsda with custom legend - ggplot2",
      fig = plotLoadings(splsda.obj, contrib = "max", style = "ggplot2", 
                        legend.title = "Treatment", legend.color = c(1:4)))
  ))
  
  # with custom title and labels
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot splsda with custom title and labels - ggplot2",
      fig = plotLoadings(splsda.obj, contrib = "max", style = "ggplot2", 
                        title = "Liver Data", X.label = "Loading", Y.label = "Genes"))
  ))
}) 
  