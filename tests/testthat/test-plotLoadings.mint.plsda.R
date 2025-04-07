context("plotLoadings.mint.plsda")

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
    png(tempfile(), width = 1200, height = 1000, res = 150)
    old_par <- par(no.readonly = TRUE)  # Save current par settings
    par(mar = c(8, 4, 4, 2))  # Increase bottom margin to fit long names
    expect_silent(invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, style = "graphics"))))
    par(old_par)  # Restore original par settings
    dev.off()
    
    # Test with ggplot2 style
    png(tempfile(), width = 1200, height = 1000, res = 150)
    old_par <- par(no.readonly = TRUE)  # Save current par settings
    par(mar = c(8, 4, 4, 2))  # Increase bottom margin to fit long names
    expect_silent(invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, style = "ggplot2"))))
    par(old_par)  # Restore original par settings
    dev.off()
    
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
    png(tempfile(), width = 1200, height = 1000, res = 150)
    old_par <- par(no.readonly = TRUE)  # Save current par settings
    par(mar = c(8, 4, 4, 2))  # Increase bottom margin to fit long names
    expect_warning(invisible(capture.output(plotLoadings(mint.splsda.obj, comp = 1, block = "X"))),
                  "'block' argument is not used for mint.plsda or mint.splsda objects")
    par(old_par)  # Restore original par settings
    dev.off()
})

# Test return value
test_that("plotLoadings.mint.plsda returns correct structure", {
    # Test return value for single study
    png(tempfile(), width = 1200, height = 1000, res = 150)
    old_par <- par(no.readonly = TRUE)  # Save current par settings
    par(mar = c(8, 4, 4, 2))  # Increase bottom margin to fit long names
    result <- plotLoadings(mint.splsda.obj, comp = 1, study = 1, contrib = "max")
    par(old_par)  # Restore original par settings
    dev.off()
    expect_true(is.data.frame(result[[1]]))
    expect_true(all(c("importance", "color", "GroupContrib") %in% names(result[[1]])))
    
    # Test return value for multiple studies
    png(tempfile(), width = 1200, height = 1000, res = 150)
    old_par <- par(no.readonly = TRUE)  # Save current par settings
    par(mar = c(8, 4, 4, 2))  # Increase bottom margin to fit long names
    result <- plotLoadings(mint.splsda.obj, comp = 1, study = c(1, 2))
    par(old_par)  # Restore original par settings
    dev.off()
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