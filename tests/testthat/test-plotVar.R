
###############################################################################
### ============================ MISSING TESTS ============================ ###
###############################################################################

# edge cases
## "We detected negative correlation between the variates of some blocks, which means that some clusters of variables observed on the correlation circle plot are not necessarily positively correlated."

###############################################################################
### ============================ GROUND TRUTHS ============================ ###
###############################################################################

Test.Data <- readRDS(system.file("testdata", "testdata-plotVar.rda", package = "mixOmics"))
Testable.Components <- Test.Data$tc
Ground.Truths <- Test.Data$gt

###############################################################################
### ================================ BASIC ================================ ###
###############################################################################


test_that("(plotVar:basic): pca", {
    
    testable.components <- Testable.Components$basic.pca
    GT <- Ground.Truths$basic.pca
    
    data(multidrug)
    X <- multidrug$ABC.trans[,1:10]
    
    res.pca <- pca(X)
    
    pca.plotVar = plotVar(res.pca)
    
    invisible(capture.output(TT <- dput(pca.plotVar[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(plotVar:basic): spca", {
    
    testable.components <- Testable.Components$basic.spca
    GT <- Ground.Truths$basic.spca
    
    data(multidrug)
    X <- multidrug$ABC.trans[,1:10]
    
    choice.keepX = c(3,3)
    
    res.spca <- spca(X, keepX = choice.keepX)
    
    spca.plotVar = plotVar(res.spca)
    
    invisible(capture.output(TT <- dput(spca.plotVar[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(plotVar:basic): pls", {
    
    testable.components <- Testable.Components$basic.pls
    GT <- Ground.Truths$basic.pls
    
    data(liver.toxicity)
    X <- liver.toxicity$gene[,1:10]
    Y <- liver.toxicity$clinic[,1:10]
    
    res.pls <- pls(X, Y)
    
    pls.plotVar = plotVar(res.pls)
    
    invisible(capture.output(TT <- dput(pls.plotVar[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(plotVar:basic): spls", {
    
    testable.components <- Testable.Components$basic.spls
    GT <- Ground.Truths$basic.spls
    
    data(liver.toxicity)
    X <- liver.toxicity$gene[,1:10]
    Y <- liver.toxicity$clinic[,1:10]
    
    choice.keepX = c(3,3)
    
    res.spls <- spls(X, Y, keepX = choice.keepX)
    
    spls.plotVar = plotVar(res.spls)
    
    invisible(capture.output(TT <- dput(spls.plotVar[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(plotVar:basic): rcc", {
    
    testable.components <- Testable.Components$basic.rcc
    GT <- Ground.Truths$basic.rcc
    
    data(nutrimouse)
    X <- nutrimouse$lipid[,1:10]
    Y <- nutrimouse$gene[,1:10]
    
    res.rcc <- rcc(X, Y)
    
    rcc.plotVar = plotVar(res.rcc)
    
    invisible(capture.output(TT <- dput(rcc.plotVar[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(plotVar:basic): plsda", {
    
    testable.components <- Testable.Components$basic.plsda
    GT <- Ground.Truths$basic.plsda
    
    data(srbct)
    X <- srbct$gene[,1:10]
    Y <- srbct$class
    
    res.plsda <- plsda(X, Y)
    
    plsda.plotVar = plotVar(res.plsda)
    
    invisible(capture.output(TT <- dput(plsda.plotVar[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(plotVar:basic): splsda", {
    
    testable.components <- Testable.Components$basic.splsda
    GT <- Ground.Truths$basic.splsda
    
    data(srbct)
    X <- srbct$gene[,1:10]
    Y <- srbct$class
    
    choice.keepX = c(3,3)
    
    res.splsda <- splsda(X, Y, keepX = choice.keepX)
    
    splsda.plotVar = plotVar(res.splsda)
    
    invisible(capture.output(TT <- dput(splsda.plotVar[testable.components])))
    
    expect_equal(TT, GT)
})


# test_that("(plotVar:basic): block.pls", {
#     
#     testable.components <- Testable.Components$basic.block.pls
#     GT <- Ground.Truths$basic.block.pls
#     
#     data(breast.TCGA)
#     X = list(miRNA = breast.TCGA$data.train$mirna[,1:10],
#             mRNA = breast.TCGA$data.train$mrna[,1:10],
#             proteomics = breast.TCGA$data.train$protein[,1:10])
#     
#     res.block.pls <- block.pls(X, indY=3)
#     
#     block.pls.plotVar = plotVar(res.block.pls)
#     
#     invisible(capture.output(TT <- dput(block.pls.plotVar[testable.components])))
#     
#     expect_equal(TT, GT)
# })
# 
# 
# test_that("(plotVar:basic): block.spls", {
#     
#     testable.components <- Testable.Components$basic.block.spls
#     GT <- Ground.Truths$basic.block.spls
#     
#     data(breast.TCGA)
#     X = list(miRNA = breast.TCGA$data.train$mirna[,1:10],
#             mRNA = breast.TCGA$data.train$mrna[,1:10],
#             proteomics = breast.TCGA$data.train$protein[,1:10])
#     
#     choice.keepX = list(miRNA = c(3,3),
#                         mRNA = c(3,3),
#                         proteomics = c(3,3))
#     
#     res.block.spls <- block.spls(X, indY=3, keepX = choice.keepX)
#     
#     block.spls.plotVar = plotVar(res.block.spls)
#     
#     invisible(capture.output(TT <- dput(block.spls.plotVar[testable.components])))
#     
#     expect_equal(TT, GT)
# })
# 
# 
# test_that("(plotVar:basic): block.plsda", {
#     
#     testable.components <- Testable.Components$basic.block.plsda
#     GT <- Ground.Truths$basic.block.plsda
#     
#     data(breast.TCGA)
#     X = list(miRNA = breast.TCGA$data.train$mirna[,1:10],
#              mRNA = breast.TCGA$data.train$mrna[,1:10],
#              proteomics = breast.TCGA$data.train$protein[,1:10])
#     Y = breast.TCGA$data.train$subtype
#     
#     res.block.plsda <- block.plsda(X, Y)
#     
#     block.plsda.plotVar = plotVar(res.block.plsda)
#     
#     invisible(capture.output(TT <- dput(block.plsda.plotVar[testable.components])))
#     
#     expect_equal(TT, GT)
# })
# 
# 
# test_that("(plotVar:basic): block.splsda", {
#     
#     testable.components <- Testable.Components$basic.block.splsda
#     GT <- Ground.Truths$basic.block.splsda
#     
#     data(breast.TCGA)
#     X = list(miRNA = breast.TCGA$data.train$mirna[,1:10],
#              mRNA = breast.TCGA$data.train$mrna[,1:10],
#              proteomics = breast.TCGA$data.train$protein[,1:10])
#     Y = breast.TCGA$data.train$subtype
#     
#     choice.keepX = list(miRNA = c(3,3),
#                         mRNA = c(3,3),
#                         proteomics = c(3,3))
#     
#     res.block.splsda <- block.splsda(X, Y, keepX = choice.keepX)
#     
#     block.splsda.plotVar = plotVar(res.block.splsda)
#     
#     invisible(capture.output(TT <- dput(block.splsda.plotVar[testable.components])))
#     
#     expect_equal(TT, GT)
# })


test_that("(plotVar:basic): mint.pls", {
    
    testable.components <- Testable.Components$basic.mint.pls
    GT <- Ground.Truths$basic.mint.pls
    
    data(stemcells)
    samples <- c(1:5,60:64)
    X <- stemcells$gene[samples, 1:10]
    Y <- stemcells$gene[samples+5, 1:10]
    S <- as.character(stemcells$study[samples])
    
    res.mint.pls <- mint.pls(X, Y, study = S)
    
    mint.pls.plotVar = plotVar(res.mint.pls)
    
    invisible(capture.output(TT <- dput(mint.pls.plotVar[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(plotVar:basic): mint.spls", {
    
    testable.components <- Testable.Components$basic.mint.spls
    GT <- Ground.Truths$basic.mint.spls
    
    data(stemcells)
    samples <- c(1:5,60:64)
    X <- stemcells$gene[samples, 1:10]
    Y <- stemcells$gene[samples+5, 1:10]
    S <- as.character(stemcells$study[samples])
    
    choice.keepX = c(3,3)
    choice.keepY = c(3,3)
    
    res.mint.spls <- mint.spls(X, Y, study = S, keepX = choice.keepX, keepY = choice.keepY)
    
    mint.spls.plotVar = plotVar(res.mint.spls)
    
    invisible(capture.output(TT <- dput(mint.spls.plotVar[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(plotVar:basic): mint.plsda", {
    
    testable.components <- Testable.Components$basic.mint.plsda
    GT <- Ground.Truths$basic.mint.plsda
    
    data(stemcells)
    samples <- c(1:5,60:64)
    X <- stemcells$gene[samples, 1:10]
    Y <- rep(c("hESC", "hiPSC"),5)
    S <- as.character(stemcells$study[samples])
    
    res.mint.plsda <- mint.plsda(X, Y, study = S)
    
    mint.plsda.plotVar = plotVar(res.mint.plsda)
    
    invisible(capture.output(TT <- dput(mint.plsda.plotVar[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(plotVar:basic): mint.splsda", {
    
    testable.components <- Testable.Components$basic.mint.splsda
    GT <- Ground.Truths$basic.mint.splsda
    
    data(stemcells)
    samples <- c(1:5,60:64)
    X <- stemcells$gene[samples, 1:10]
    Y <- rep(c("hESC", "hiPSC"),5)
    S <- as.character(stemcells$study[samples])
    
    choice.keepX = c(3,3)
    
    res.mint.splsda <- mint.splsda(X, Y, study = S, keepX = choice.keepX)
    
    mint.splsda.plotVar = plotVar(res.mint.splsda)
    
    invisible(capture.output(TT <- dput(mint.splsda.plotVar[testable.components])))
    
    expect_equal(TT, GT)
})


###############################################################################
### ============================== PARAMETER ============================== ###
###############################################################################


test_that("(plotVar:parameter): comp", {
    
    testable.components <- Testable.Components$comp.pca
    GT <- Ground.Truths$comp.pca
    
    data(multidrug)
    X <- multidrug$ABC.trans[,1:10]
    
    res.pca <- pca(X, ncomp=3)
    
    comp.plotVar = plotVar(res.pca, comp = c(2,3))
    
    invisible(capture.output(TT <- dput(comp.plotVar[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(plotVar:parameter): comp.select", {
    
    testable.components <- Testable.Components$comp.select.pca
    GT <- Ground.Truths$comp.select.pca
    
    data(multidrug)
    X <- multidrug$ABC.trans[,1:10]
    
    res.pca <- pca(X, ncomp=3)
    
    comp.select.plotVar = plotVar(res.pca, comp.select = c(2,3))
    
    invisible(capture.output(TT <- dput(comp.select.plotVar[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(plotVar:parameter): blocks", {
    
    testable.components <- Testable.Components$blocks.block.splsda
    GT <- Ground.Truths$blocks.block.splsda
    
    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.train$mirna[,1:10],
             mRNA = breast.TCGA$data.train$mrna[,1:10],
             proteomics = breast.TCGA$data.train$protein[,1:10])
    Y = breast.TCGA$data.train$subtype
    
    choice.keepX = list(miRNA = c(3,3),
                        mRNA = c(3,3),
                        proteomics = c(3,3))
    
    res.block.splsda <- block.splsda(X, Y, keepX = choice.keepX)
    
    blocks.plotVar = plotVar(res.block.splsda, blocks = c("miRNA"))
    
    invisible(capture.output(TT <- dput(blocks.plotVar[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(plotVar:parameter): cutoff", {
    
    testable.components <- Testable.Components$cutoff.pca
    GT <- Ground.Truths$cutoff.pca
    
    data(multidrug)
    X <- multidrug$ABC.trans[,1:10]
    
    res.pca <- pca(X, ncomp=3)
    
    cutoff.plotVar = plotVar(res.pca, cutoff = 0.5)
    
    invisible(capture.output(TT <- dput(cutoff.plotVar[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(plotVar:parameter): col", {
    
    testable.components <- Testable.Components$col.pls
    GT <- Ground.Truths$col.pls
    
    data(liver.toxicity)
    X <- liver.toxicity$gene[,1:5]
    Y <- liver.toxicity$clinic[,1:5]
    
    res.pls <- pls(X, Y)
    
    col.plotVar = plotVar(res.pls, col = c("purple", "green"))
    
    invisible(capture.output(TT <- dput(col.plotVar[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(plotVar:parameter): pch", {
    
    testable.components <- Testable.Components$pch.pls
    GT <- Ground.Truths$pch.pls
    
    data(liver.toxicity)
    X <- liver.toxicity$gene[,1:5]
    Y <- liver.toxicity$clinic[,1:5]
    
    res.pls <- pls(X, Y)
    
    pch.plotVar = plotVar(res.pls, pch = c(1, 2))
    
    invisible(capture.output(TT <- dput(pch.plotVar[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(plotVar:parameter): cex", {
    
    testable.components <- Testable.Components$cex.pls
    GT <- Ground.Truths$cex.pls
    
    data(liver.toxicity)
    X <- liver.toxicity$gene[,1:5]
    Y <- liver.toxicity$clinic[,1:5]
    
    res.pls <- pls(X, Y)
    
    cex.plotVar = plotVar(res.pls, cex = c(1, 2))
    
    invisible(capture.output(TT <- dput(cex.plotVar[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(plotVar:parameter): font", {
    
    testable.components <- Testable.Components$font.pls
    GT <- Ground.Truths$font.pls
    
    data(liver.toxicity)
    X <- liver.toxicity$gene[,1:5]
    Y <- liver.toxicity$clinic[,1:5]
    
    res.pls <- pls(X, Y)
    
    f <- list(X=rep(4,5), 
              Y=rep(3,5))
    
    font.plotVar = plotVar(res.pls, font = f)
    
    invisible(capture.output(TT <- dput(font.plotVar[testable.components])))
    
    expect_equal(TT, GT)
})


###############################################################################
### ================================ ERROR ================================ ###
###############################################################################


test_that("(plotVar:error): catches invalid `style` values", {
    
    data(multidrug)
    X <- multidrug$ABC.trans[,1:10]
    
    res.pca <- pca(X)
    
    expect_error(plotVar(res.pca, style = "random.style"),
                 "'style' must be one of 'ggplot2', '3d' , lattice' or 'graphics'.",
                 fixed=T)
})


test_that("(plotVar:error): catches invalid `plot`` values", {
    
    data(multidrug)
    X <- multidrug$ABC.trans[,1:10]
    
    res.pca <- pca(X)
    
    expect_error(plotVar(res.pca, plot = "random.value"),
                 "'plot' must be logical.",
                 fixed=T)
})


test_that("(plotVar:error): catches mismatching `block` names", {
    
    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.train$mirna[,1:10],
             mRNA = breast.TCGA$data.train$mrna[,1:10],
             proteomics = breast.TCGA$data.train$protein[,1:10])
    Y = breast.TCGA$data.train$subtype
    
    res.block.plsda <- block.plsda(X, Y)
    
    expect_error(plotVar(res.block.plsda, blocks = "random.block"),
                 "One element of 'blocks' does not match with the names of the blocks",
                 fixed=T)
})


test_that("(plotVar:error): a given block has only one component selected", {
    
    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.train$mirna[,1:10],
             mRNA = breast.TCGA$data.train$mrna[,1:10],
             proteomics = breast.TCGA$data.train$protein[,1:10])
    Y = breast.TCGA$data.train$subtype
    
    choice.keepX = list(miRNA = 3,
                        mRNA = 3,
                        proteomics = 3)
    
    res.block.splsda <- block.splsda(X, Y, ncomp = 1, keepX = choice.keepX)
    
    expect_error(plotVar(res.block.splsda),
                 "The number of components for one selected block ' miRNA - mRNA - proteomics ' is 1. The number of components must be superior or equal to 2.",
                 fixed=T)
})


test_that("(plotVar:error): catches invalid `rad.in` values", {
    
    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.train$mirna[,1:10],
             mRNA = breast.TCGA$data.train$mrna[,1:10],
             proteomics = breast.TCGA$data.train$protein[,1:10])
    Y = breast.TCGA$data.train$subtype
    
    choice.keepX = list(miRNA = 3,
                        mRNA = 3,
                        proteomics = 3)
    
    res.block.splsda <- block.splsda(X, Y, keepX = choice.keepX)
    
    expect_error(plotVar(res.block.splsda, rad.in = "random.value"),
                 "The value taken by 'rad.in' must be between 0 and 1",
                 fixed=T)
    
    expect_error(plotVar(res.block.splsda, rad.in = -1),
                 "The value taken by 'rad.in' must be between 0 and 1",
                 fixed=T)
    
    expect_error(plotVar(res.block.splsda, rad.in = 2),
                 "The value taken by 'rad.in' must be between 0 and 1",
                 fixed=T)
})


test_that("(plotVar:error): catches invalid `cutoff` values", {
    
    data(multidrug)
    X <- multidrug$ABC.trans[,1:10]
    
    res.pca <- pca(X)
    
    expect_error(plotVar(res.pca, cutoff = "random.value"),
                 "The value taken by 'cutoff' must be between 0 and 1",
                 fixed=T)
    
    expect_error(plotVar(res.pca, cutoff = -1),
                 "The value taken by 'cutoff' must be between 0 and 1",
                 fixed=T)
    
    expect_error(plotVar(res.pca, cutoff = 2),
                 "The value taken by 'cutoff' must be between 0 and 1",
                 fixed=T)
})


test_that("(plotVar:error): catches invalid `comp` values", {
    
    data(multidrug)
    X <- multidrug$ABC.trans[,1:10]
    
    res.pca <- pca(X, ncomp = 3)
    
    expect_error(plotVar(res.pca, comp = c(3,4,5)),
                 "'comp' must be a numeric vector of length 2.",
                 fixed=T)
    
    expect_error(plotVar(res.pca, style = "3d", comp = c(3,4)),
                 "'comp' must be a numeric vector of length 3.",
                 fixed=T)
    
    expect_error(plotVar(res.pca, comp = c(4,5)),
                 "Each element of 'comp' must be positive <= 3.",
                 fixed=T)
    
    expect_error(plotVar(res.pca, comp = c(-1,-2)),
                 "Each element of 'comp' must be positive <= 3.",
                 fixed=T)
})


test_that("(plotVar:error): catches invalid `comp.selected` values", {
    
    data(multidrug)
    X <- multidrug$ABC.trans[,1:10]
    
    res.pca <- pca(X, ncomp = 3)
    
    expect_error(plotVar(res.pca, comp.select = c(4,5)),
                 "Each element of 'comp.select' must be positive and <= 3.",
                 fixed=T)
    
    expect_error(plotVar(res.pca, comp.select = c(-1,-2)),
                 "Each element of 'comp.select' must be positive and <= 3.",
                 fixed=T)
})


###############################################################################
### =============================== WARNINGS ============================== ###
###############################################################################


test_that("(plotVar:warning): notifying user of negative correlation between variates of differing blocks", {
    
    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.train$mirna[,1:10],
             mRNA = breast.TCGA$data.train$mrna[,1:10],
             proteomics = breast.TCGA$data.train$protein[,1:10])
    Y = breast.TCGA$data.train$subtype
    
    choice.keepX = list(miRNA = 3,
                        mRNA = 3,
                        proteomics = 3)
    
    res.block.splsda <- block.splsda(X, Y, ncomp = 3, keepX = choice.keepX)
    
    # induce negative correlation between variates of differing blocks
    res.block.splsda$variates$miRNA[,1] <- -res.block.splsda$variates$mRNA[,1]
    
    expect_warning(plotVar(res.block.splsda),
                   "We detected negative correlation between the variates of some blocks, which means that some clusters of variables observed on the correlation circle plot are not necessarily positively correlated.",
                   fixed=T)
})


###############################################################################

dev.off()
unlink(list.files(pattern = "*.pdf"))
