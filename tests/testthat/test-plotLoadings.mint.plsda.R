context("plotLoadings.mint.plsda")

# Load required libraries
library(mixOmics)
library(ggplot2)
library(gridExtra)

# Load data
data(stemcells)
data = stemcells$gene
type.id = stemcells$celltype
study.id = stemcells$study

# Create test object
mint.splsda.obj = mint.splsda(X = data, Y = type.id, study = study.id, ncomp = 2, keepX = c(10, 5))

# Basic Unit Tests
test_that("plotLoadings.mint.plsda works", {
    # Test default behavior with graphics style
    expect_silent(invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, style = "graphics"))))
    
    # Test with ggplot2 style
    expect_silent(invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, style = "ggplot2"))))
    
    # Test with invalid style
    expect_error(plotLoadings(mint.splsda.obj, comp = 1, style = "invalid"),
                 "'style' must be either 'graphics' or 'ggplot2'")
    
    # Test with invalid component
    expect_error(plotLoadings(mint.splsda.obj, comp = 0),
                 "'comp' must be a positive integer")
    
    # Test with invalid ndisplay
    expect_error(plotLoadings(mint.splsda.obj, comp = 1, ndisplay = 0),
                 "'ndisplay' must be a positive integer")
    
    # Test with invalid study
    expect_error(plotLoadings(mint.splsda.obj, comp = 1, study = "invalid"),
                 "'study' must from one of 'object\\$study', 'global' or 'all.partial', see help file")
    
    # Test with invalid block argument
    expect_warning(invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, block = "X"))),
                  "'block' argument is not used for mint.plsda or mint.splsda objects")
})

# Test study-specific functionality
test_that("plotLoadings.mint.plsda handles different study options", {
    # Test single study
    expect_silent(invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, study = 2))))
    
    # Test multiple studies
    expect_silent(invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, study = c(2, 3)))))
    
    # Test all.partial option
    expect_silent(invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, study = "all.partial"))))
    
    # Test global option
    expect_silent(invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, study = "global"))))
})

# Test customization options
test_that("plotLoadings.mint.plsda handles customization options", {
    # Test with custom title
    expect_silent(invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, title = "Custom Title"))))
    
    # Test with custom subtitle
    expect_silent(invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, study = c(1, 2), 
                             subtitle = c("Subtitle 1", "Subtitle 2")))))
    
    # Test with custom sizes
    expect_silent(invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, 
                             size.title = 2, 
                             size.subtitle = 1.6,
                             size.name = 0.7,
                             size.axis = 0.7,
                             size.labs = 1))))
    
    # Test with custom labels
    expect_silent(invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, 
                             X.label = "Custom X Label",
                             Y.label = "Custom Y Label"))))
    
    # Test with custom colors
    expect_silent(invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, 
                             legend.color = c("red", "blue", "green")))))
    
    # Test with custom legend title
    expect_silent(invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, 
                             legend.title = "Custom Legend"))))
})

# Test display options
test_that("plotLoadings.mint.plsda handles display options", {
    # Test with ndisplay
    expect_silent(invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, ndisplay = 5))))
    
    # Test with xlim
    expect_silent(invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, xlim = c(-1, 1)))))
    
    # Test with layout
    expect_silent(invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, study = c(3, 4), 
                             layout = c(1, 2)))))
})

# Test contribution options
test_that("plotLoadings.mint.plsda handles contribution options", {
    # Test with max contribution
    expect_silent(invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, contrib = "max"))))
    
    # Test with min contribution
    expect_silent(invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, contrib = "min"))))
    
    # Test with mean method
    expect_silent(invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, contrib = "max", method = "mean"))))
    
    # Test with median method
    expect_silent(invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, contrib = "max", method = "median"))))
    
    # Test with ties handling
    expect_silent(invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, contrib = "max", 
                             show.ties = TRUE, col.ties = "white"))))
})

# Test return value
test_that("plotLoadings.mint.plsda returns correct structure", {
    # Test return value for single study
    result <- plotLoadings(mint.splsda.obj, comp = 1, study = 1)
    expect_true(is.data.frame(result[[1]]))
    expect_true(all(c("importance", "color", "names") %in% names(result[[1]])))
    
    # Test return value for multiple studies
    result <- plotLoadings(mint.splsda.obj, comp = 1, study = c(1, 2))
    expect_true(is.list(result))
    expect_equal(length(result), 2)
    expect_true(all(sapply(result, is.data.frame)))
})

# vdiffr tests for graphics style
test_that("plotLoadings.mint.plsda graphics style plots match", {
    skip_on_ci()
    
    # Simple plot
    vdiffr::expect_doppelganger(
        "mint.plsda.loadings.simple",
        invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, style = "graphics")))
    )
    
    # Plot with all.partial
    vdiffr::expect_doppelganger(
        "mint.plsda.loadings.all.partial",
        invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, study = "all.partial", style = "graphics")))
    )
    
    # Plot for specific study
    vdiffr::expect_doppelganger(
        "mint.plsda.loadings.specific.study",
        invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, study = 1, style = "graphics")))
    )
    
    # Change gene names
    vdiffr::expect_doppelganger(
        "mint.plsda.loadings.custom.names",
        invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, name.var = paste0("Gene_", 1:ncol(data)), style = "graphics")))
    )
    
    # Change labels and label sizes
    vdiffr::expect_doppelganger(
        "mint.plsda.loadings.custom.labels",
        invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, 
                    X.label = "Custom X", Y.label = "Custom Y",
                    size.name = 1, size.axis = 1, size.labs = 1.5, style = "graphics")))
    )
    
    # Change colors and borders
    vdiffr::expect_doppelganger(
        "mint.plsda.loadings.custom.colors",
        invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, 
                    legend.color = c("red", "blue", "green"),
                    border = "black", style = "graphics")))
    )
    
    # Change layout
    vdiffr::expect_doppelganger(
        "mint.plsda.loadings.custom.layout",
        invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, study = c(1, 2), 
                    layout = c(1, 2), style = "graphics")))
    )
})

# vdiffr tests for ggplot2 style
test_that("plotLoadings.mint.plsda ggplot2 style plots match", {
    skip_on_ci()
    
    # Simple plot
    vdiffr::expect_doppelganger(
        "mint.plsda.loadings.ggplot2.simple",
        invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, style = "ggplot2")))
    )
    
    # Plot with all.partial
    vdiffr::expect_doppelganger(
        "mint.plsda.loadings.ggplot2.all.partial",
        invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, study = "all.partial", style = "ggplot2")))
    )
    
    # Plot for specific study
    vdiffr::expect_doppelganger(
        "mint.plsda.loadings.ggplot2.specific.study",
        invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, study = 3, style = "ggplot2")))
    )
    
    # Change gene names
    vdiffr::expect_doppelganger(
        "mint.plsda.loadings.ggplot2.custom.names",
        invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, name.var = paste0("Gene_", 1:ncol(data)), style = "ggplot2")))
    )
    
    # Change labels and label sizes
    vdiffr::expect_doppelganger(
        "mint.plsda.loadings.ggplot2.custom.labels",
        invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, 
                    X.label = "Custom X", Y.label = "Custom Y",
                    size.name = 1, size.axis = 1, size.labs = 1.5, style = "ggplot2")))
    )
    
    # Change colors and borders
    vdiffr::expect_doppelganger(
        "mint.plsda.loadings.ggplot2.custom.colors",
        invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, 
                    legend.color = c("red", "blue", "green"),
                    border = "black", style = "ggplot2")))
    )
    
    # Change layout
    vdiffr::expect_doppelganger(
        "mint.plsda.loadings.ggplot2.custom.layout",
        invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, study = c(2, 3), 
                    layout = c(1, 2), style = "ggplot2")))
    )
}) 