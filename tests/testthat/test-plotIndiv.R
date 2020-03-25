context("plotIndiv")

test_that("plotIndiv.rcc works without ind.names", code = {
    data(nutrimouse)
    X <- nutrimouse$lipid
    Y <- nutrimouse$gene
    nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
    
    plotIndiv.res <- tryCatch(expr = {plotIndiv(nutri.res, group = nutrimouse$genotype, ind.names = FALSE, legend = TRUE)},
                              error = function(e) e,
                              warning = function(w) w
    )
    
    expect_is(plotIndiv.res$graph, "ggplot")
})

## ------------------------------------------------------------------------ ##

test_that("plotIndiv.sgcca(..., blocks = 'consensus') works", code = {
    data(nutrimouse)
    Y = unmap(nutrimouse$diet)
    data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)
    design1 = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
    nutrimouse.sgcca <- wrapper.sgcca(X = data,
                                      design = design1,
                                      penalty = c(0.3, 0.5, 1),
                                      ncomp = 2,
                                      scheme = "horst")
    
    # default style: one panel for each block
    plotindiv_res <- plotIndiv(nutrimouse.sgcca, blocks = c("lipid","consensus"))
    
    expect_true(any(grepl(pattern = "Consensus", x = unique(plotindiv_res$df$Block))))
})

## ------------------------------------------------------------------------ ##

test_that("plotIndiv.sgccda(..., blocks = 'consensus') works with ind.names and ell", code = {
    data("breast.TCGA")
    data = list(mrna = breast.TCGA$data.train$mrna, mirna = breast.TCGA$data.train$mirna,
                protein = breast.TCGA$data.train$protein)
    design = matrix(1, ncol = length(data), nrow = length(data),
                    dimnames = list(names(data), names(data)))
    diag(design) =  0
    # set number of variables to select, per component and per data set (this is set arbitrarily)
    list.keepX = list(mrna = rep(4, 2), mirna = rep(5,2), protein = rep(5, 2))
    TCGA.block.splsda = block.splsda(X = data, Y = breast.TCGA$data.train$subtype,
                                     ncomp = 2, keepX = list.keepX, design = design)
    blocks <- c("consensus", "mrna", "weighted.consensus")
    diablo_plot <- plotIndiv(TCGA.block.splsda, ind.names = FALSE, blocks = blocks)
    expect_true(all(unique(diablo_plot$df$Block) %in% c('Consensus', 'Block: mrna', 'Consensus (weighted)')))
})

## ------------------------------------------------------------------------ ##

test_that("plotIndiv.sgccda(..., blocks = 'consensus') works with ellipse=TRUE", code = {
    data("breast.TCGA")
    data = list(mrna = breast.TCGA$data.train$mrna, mirna = breast.TCGA$data.train$mirna,
                protein = breast.TCGA$data.train$protein)
    design = matrix(1, ncol = length(data), nrow = length(data),
                    dimnames = list(names(data), names(data)))
    diag(design) =  0
    # set number of variables to select, per component and per data set (this is set arbitrarily)
    list.keepX = list(mrna = rep(4, 2), mirna = rep(5,2), protein = rep(5, 2))
    TCGA.block.splsda = block.splsda(X = data, Y = breast.TCGA$data.train$subtype,
                                     ncomp = 2, keepX = list.keepX, design = design)
    blocks <- c("consensus", "mrna", "weighted.consensus")
    diablo_plot <- plotIndiv(TCGA.block.splsda, ind.names = TRUE, blocks = blocks, ellipse = TRUE)
    expect_true(all(unique(diablo_plot$df.ellipse$Block) %in% c('Consensus', 'Block: mrna', 'Consensus (weighted)')))
})