

###############################################################################
### ============================ MISSING TESTS ============================ ###
###############################################################################

# edge cases
## "'col' is ignored as 'group' has been set."

###############################################################################
### ============================ GROUND TRUTHS ============================ ###
###############################################################################

Test.Data <- readRDS(system.file("testdata", "testdata-plotIndiv.rda", package = "mixOmics"))
Testable.Components <- Test.Data$tc
Ground.Truths <- Test.Data$gt

###############################################################################
### ================================ BASIC ================================ ###
###############################################################################


test_that("(plotIndiv:basic): pca", {
    
    testable.components <- Testable.Components$basic.pca
    GT <- Ground.Truths$basic.pca
    
    data(multidrug)
    X <- multidrug$ABC.trans[1:10, 1:10]
    
    res.pca <- pca(X)
    
    pca.plotIndiv = plotIndiv(res.pca)$graph
    
    invisible(capture.output(TT <- dput(pca.plotIndiv[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(plotIndiv:basic): spca", {
    
    testable.components <- Testable.Components$basic.spca
    GT <- Ground.Truths$basic.spca
    
    data(multidrug)
    X <- multidrug$ABC.trans[1:10, 1:10]
    
    choice.keepX <- c(3,3)
    
    res.spca <- spca(X, keepX = choice.keepX)
    
    spca.plotIndiv = plotIndiv(res.spca)$graph
    
    invisible(capture.output(TT <- dput(spca.plotIndiv[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(plotIndiv:basic): pls", {
    
    testable.components <- Testable.Components$basic.pls
    GT <- Ground.Truths$basic.pls
    
    data(liver.toxicity)
    X <- liver.toxicity$gene[1:10, ]
    Y <- liver.toxicity$clinic[1:10, ]
    
    res.pls <- pls(X, Y)
    
    pls.plotIndiv = plotIndiv(res.pls)$graph
    
    invisible(capture.output(TT <- dput(pls.plotIndiv[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(plotIndiv:basic): spls", {
    
    testable.components <- Testable.Components$basic.spls
    GT <- Ground.Truths$basic.spls
    
    data(liver.toxicity)
    X <- liver.toxicity$gene[1:10, ]
    Y <- liver.toxicity$clinic[1:10, ]
    
    choice.keepX <- c(3,3)
    
    res.spls <- spls(X, Y, keepX = choice.keepX)
    
    spls.plotIndiv = plotIndiv(res.spls)$graph
    
    invisible(capture.output(TT <- dput(spls.plotIndiv[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(plotIndiv:basic): rcc", {
    
    testable.components <- Testable.Components$basic.rcc
    GT <- Ground.Truths$basic.rcc
    
    data(nutrimouse)
    X <- nutrimouse$lipid[1:10, ]
    Y <- nutrimouse$gene[1:10, ]
    
    res.rcc <- rcc(X, Y, method = "shrinkage")
    
    rcc.plotIndiv = plotIndiv(res.rcc)$graph
    
    invisible(capture.output(TT <- dput(rcc.plotIndiv[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(plotIndiv:basic): plsda", {
    
    testable.components <- Testable.Components$basic.plsda
    GT <- Ground.Truths$basic.plsda
    
    data(srbct)
    samples <- c(1:3, 24:26, 41:43)
    X <- srbct$gene[samples,]
    Y <- srbct$class[samples]
    
    res.plsda <- plsda(X, Y)
    
    plsda.plotIndiv = plotIndiv(res.plsda)$graph
    
    invisible(capture.output(TT <- dput(plsda.plotIndiv[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(plotIndiv:basic): splsda", {
    
    testable.components <- Testable.Components$basic.splsda
    GT <- Ground.Truths$basic.splsda
    
    data(srbct)
    samples <- c(1:3, 24:26, 41:43)
    X <- srbct$gene[samples,]
    Y <- srbct$class[samples]
    
    choice.keepX = c(3,3)
    
    res.splsda <- splsda(X, Y, keepX = choice.keepX)
    
    splsda.plotIndiv = plotIndiv(res.splsda)$graph
    
    invisible(capture.output(TT <- dput(splsda.plotIndiv[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(plotIndiv:basic): block.plsda", {
    
    testable.components <- Testable.Components$basic.block.plsda
    GT <- Ground.Truths$basic.block.plsda
    
    data(breast.TCGA)
    samples <- c(1:3, 22:25, 36:38)
    X = list(miRNA = breast.TCGA$data.test$mirna[samples,],
             mRNA = breast.TCGA$data.test$mrna[samples,])
    Y = breast.TCGA$data.test$subtype[samples]
    
    res.block.plsda <- block.plsda(X, Y)
    
    block.plsda.plotIndiv = plotIndiv(res.block.plsda)$graph
    
    invisible(capture.output(TT <- dput(block.plsda.plotIndiv[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(plotIndiv:basic): block.splsda", {
    
    testable.components <- Testable.Components$basic.block.splsda
    GT <- Ground.Truths$basic.block.splsda
    
    data(breast.TCGA)
    samples <- c(1:3, 22:25, 36:38)
    X = list(miRNA = breast.TCGA$data.test$mirna[samples,],
             mRNA = breast.TCGA$data.test$mrna[samples,])
    Y = breast.TCGA$data.test$subtype[samples]
    
    choice.keepX = list(miRNA=c(3,3),
                        mRNA=c(3,3))
    
    res.block.splsda <- block.splsda(X, Y, keepX = choice.keepX)
    
    block.splsda.plotIndiv = plotIndiv(res.block.splsda)$graph
    
    invisible(capture.output(TT <- dput(block.splsda.plotIndiv[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(plotIndiv:basic): mint.pls", {
    
    testable.components <- Testable.Components$basic.mint.pls
    GT <- Ground.Truths$basic.mint.pls
    
    data(stemcells)
    samples <- c(1:5,60:64)
    X <- stemcells$gene[samples, 1:10]
    Y <- stemcells$gene[samples+5, 1:10]
    S <- as.character(stemcells$study[samples])
    
    res.mint.pls <- mint.pls(X, Y, study = S)
    
    mint.pls.plotIndiv = plotIndiv(res.mint.pls)$graph
    
    invisible(capture.output(TT <- dput(mint.pls.plotIndiv[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(plotIndiv:basic): mint.spls", {
    
    testable.components <- Testable.Components$basic.mint.spls
    GT <- Ground.Truths$basic.mint.spls
    
    data(stemcells)
    samples <- c(1:5,60:64)
    X <- stemcells$gene[samples, 1:10]
    Y <- stemcells$gene[samples+5, 1:10]
    S <- as.character(stemcells$study[samples])
    
    choice.keepX <- c(3,3)
    
    res.mint.spls <- mint.spls(X, Y, study = S, keepX = choice.keepX)
    
    mint.spls.plotIndiv = plotIndiv(res.mint.spls)$graph
    
    invisible(capture.output(TT <- dput(mint.spls.plotIndiv[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(plotIndiv:basic): mint.plsda", {
    
    testable.components <- Testable.Components$basic.mint.plsda
    GT <- Ground.Truths$basic.mint.plsda
    
    data(stemcells)
    samples <- c(1:5,60:64)
    X <- stemcells$gene[samples, 1:10]
    Y <- rep(c("hESC", "hiPSC"),5)
    S <- as.character(stemcells$study[samples])
    
    res.mint.plsda <- mint.plsda(X, Y, study = S)
    
    mint.plsda.plotIndiv = plotIndiv(res.mint.plsda)$graph
    
    invisible(capture.output(TT <- dput(mint.plsda.plotIndiv[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(plotIndiv:basic): mint.splsda", {
    
    testable.components <- Testable.Components$basic.mint.splsda
    GT <- Ground.Truths$basic.mint.splsda
    
    data(stemcells)
    samples <- c(1:5,60:64)
    X <- stemcells$gene[samples, 1:10]
    Y <- rep(c("hESC", "hiPSC"),5)
    S <- as.character(stemcells$study[samples])
    
    choice.keepX <- c(3,3)
    
    res.mint.splsda <- mint.splsda(X, Y, study = S, keepX = choice.keepX)
    
    mint.splsda.plotIndiv = plotIndiv(res.mint.splsda)$graph
    
    invisible(capture.output(TT <- dput(mint.splsda.plotIndiv[testable.components])))
    
    expect_equal(TT, GT)
})


###############################################################################
### ============================== PARAMETER ============================== ###
###############################################################################


test_that("(plotIndiv:parameter): comp", {
    
    testable.components <- Testable.Components$comp.pca
    GT <- Ground.Truths$comp.pca
    
    data(multidrug)
    X <- multidrug$ABC.trans[1:10, 1:10]
    
    res.pca <- pca(X, ncomp = 3)
    
    pca.plotIndiv = plotIndiv(res.pca, comp = c(2,3))$graph
    
    invisible(capture.output(TT <- dput(pca.plotIndiv[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(plotIndiv:parameter): study", {
    
    testable.components <- Testable.Components$study.mint.splsda
    GT <- Ground.Truths$study.mint.splsda
    
    data(stemcells)
    samples <- c(1:5,60:64)
    X <- stemcells$gene[samples, 1:10]
    Y <- rep(c("hESC", "hiPSC"),5)
    S <- as.character(stemcells$study[samples])
    
    choice.keepX <- c(3,3)
    
    res.mint.splsda <- mint.splsda(X, Y, study = S, keepX = choice.keepX)
    
    mint.splsda.plotIndiv = plotIndiv(res.mint.splsda, study = "all.partial")$graph
    
    invisible(capture.output(TT <- dput(mint.splsda.plotIndiv[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(plotIndiv:parameter): ellipse", {
    
    testable.components <- Testable.Components$ellipse.splsda
    GT <- Ground.Truths$ellipse.splsda
    
    data(srbct)
    samples <- c(1:3, 24:26, 41:43)
    X <- srbct$gene[samples,]
    Y <- srbct$class[samples]
    
    choice.keepX = c(3,3)
    
    res.splsda <- splsda(X, Y, keepX = choice.keepX)
    
    splsda.plotIndiv = plotIndiv(res.splsda, ellipse = T)
    
    invisible(capture.output(TT <- dput(splsda.plotIndiv[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(plotIndiv:parameter): blocks", {
    
    testable.components <- Testable.Components$blocks.block.plsda
    GT <- Ground.Truths$blocks.block.plsda
    
    data(breast.TCGA)
    samples <- c(1:3, 22:25, 36:38)
    X = list(miRNA = breast.TCGA$data.test$mirna[samples,],
             mRNA = breast.TCGA$data.test$mrna[samples,])
    Y = breast.TCGA$data.test$subtype[samples]
    
    res.block.plsda <- block.plsda(X, Y)
    
    block.plsda.plotIndiv = plotIndiv(res.block.plsda, blocks = "mRNA")$graph
    
    invisible(capture.output(TT <- dput(block.plsda.plotIndiv[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(plotIndiv:parameter): ind.names", {
    
    testable.components <- Testable.Components$ind.names.pca
    GT <- Ground.Truths$ind.names.pca
    
    data(multidrug)
    X <- multidrug$ABC.trans[1:10, 1:10]
    
    res.pca <- pca(X)
    
    pca.plotIndiv = plotIndiv(res.pca, ind.names = F)$graph
    
    invisible(capture.output(TT <- dput(pca.plotIndiv[testable.components])))
    
    expect_equal(TT, GT)
})


###############################################################################
### ================================ ERROR ================================ ###
###############################################################################


test_that("(plotIndiv:error): catches invalid style values", {
    
    data(multidrug)
    X <- multidrug$ABC.trans[1:10, 1:10]
    
    res.pca <- pca(X)
    
    expect_error(plotIndiv(res.pca, style = "random.style"),
                 "'style' must be one of 'ggplot2', 'lattice', 'graphics' or '3d' .",
                 fixed=T)
})


test_that("(plotIndiv:error): catches invalid ellipse/ellipse.level values", {
    
    data(srbct)
    samples <- c(1:3, 24:26, 41:43)
    X <- srbct$gene[samples,]
    Y <- srbct$class[samples]
    
    choice.keepX = c(3,3)
    
    res.splsda <- splsda(X, Y, keepX = choice.keepX)
    
    expect_error(plotIndiv(res.splsda, ellipse = "random.value"),
                 "'ellipse' must be either TRUE or FALSE",
                 fixed=T)
    
    expect_error(plotIndiv(res.splsda, ellipse.level = -1),
                 "The value taken by 'ellipse.level' must be between 0 and 1",
                 fixed=T)
    
    expect_error(plotIndiv(res.splsda, ellipse.level = 2),
                 "The value taken by 'ellipse.level' must be between 0 and 1",
                 fixed=T)
    
    expect_error(plotIndiv(res.splsda, ellipse.level = "random.value"),
                 "The value taken by 'ellipse.level' must be between 0 and 1",
                 fixed=T)
})


test_that("(plotIndiv:error): catches invalid alpha values", {
    
    data(multidrug)
    X <- multidrug$ABC.trans[1:10, 1:10]
    
    res.pca <- pca(X)
    
    expect_error(plotIndiv(res.pca, alpha = "random.value"),
                 "The value taken by 'alpha' must be between 0 and 1",
                 fixed=T)
})


test_that("(plotIndiv:error): catches invalid comp values", {
    
    data(multidrug)
    X <- multidrug$ABC.trans[1:10, 1:10]
    
    res.pca <- pca(X, ncomp = 5)
    
    expect_error(plotIndiv(res.pca, comp = c(1)),
                 "'comp' must be a numeric vector of length 2.",
                 fixed=T)
    
    expect_error(plotIndiv(res.pca, comp = c(1,2,3)),
                 "'comp' must be a numeric vector of length 2.",
                 fixed=T)
    
    expect_error(plotIndiv(res.pca, style = "3d", comp = c(1,2)),
                 "'comp' must be a numeric vector of length 3.",
                 fixed=T)
    
    expect_error(plotIndiv(res.pca, style = "3d", comp = c(1,2,3,4)),
                 "'comp' must be a numeric vector of length 3.",
                 fixed=T)
})


test_that("(plotIndiv:error): catches invalid ind.names values", {
    
    data(multidrug)
    X <- multidrug$ABC.trans[1:10, 1:10]
    
    res.pca <- pca(X, ncomp = 5)
    
    expect_error(plotIndiv(res.pca, ind.names = c("random.name.1", "random.name.2")),
                 "'ind.names' must be a character vector of length 10 or a Logical atomic vector.",
                 fixed=T)
})


test_that("(plotIndiv:error): catches invalid xlim/ylim values", {
    
    data(breast.TCGA)
    samples <- c(1:3, 22:25, 36:38)
    X = list(miRNA = breast.TCGA$data.train$mirna[samples,],
             mRNA = breast.TCGA$data.train$mrna[samples,],
             proteomics = breast.TCGA$data.train$protein[samples,])
    Y = breast.TCGA$data.test$subtype[samples]
    
    res.block.plsda <- block.plsda(X, Y)
    
    expect_error(plotIndiv(res.block.plsda, style="lattice", blocks = c("miRNA", "mRNA", "proteomics"),
                           ylim = list(miRNA = c(-15,15),
                                       mRNA = c(-15,15))),
                 "'ylim' must be a list of 3 vectors of length 2.",
                 fixed=T)
    
    expect_error(plotIndiv(res.block.plsda, blocks = "mRNA",
                           ylim = c(-30, 0, 30)),
                 "'ylim' must be a vector of length 2.",
                 fixed=T)
    
    expect_error(plotIndiv(res.block.plsda, style="lattice", blocks = c("miRNA", "mRNA", "proteomics"),
                           xlim = list(miRNA = c(-15,15),
                                       mRNA = c(-15,15))),
                 "'xlim' must be a list of 3 vectors of length 2.",
                 fixed=T)
    
    expect_error(plotIndiv(res.block.plsda, blocks = "mRNA",
                           xlim = c(-30, 0, 30)),
                 "'xlim' must be a vector of length 2.",
                 fixed=T)
})


test_that("(plotIndiv:error): catches invalid col values", {
    
    data(srbct)
    samples <- c(1:3, 24:26, 41:43)
    X <- srbct$gene[samples,]
    Y <- srbct$class[samples]
    
    res.plsda <- plsda(X, Y)
    
    expect_error(plotIndiv(res.plsda, col=c("#388ECC", "#F68B33")),
                 "Length of 'col' should be of length 3 (the number of groups).",
                 fixed=T)
})


test_that("(plotIndiv:error): catches invalid col.per.group values", {
    
    data(liver.toxicity)
    X <- liver.toxicity$gene[1:10,]
    Y <- liver.toxicity$clinic[1:10,]
    
    res.pls <- pls(X, Y)
    
    group <- c(rep("c1", 5), rep("c2", 5))
    
    expect_error(plotIndiv(res.pls, group = group, col.per.group = c("#388ECC", "#F68B33", "#cc7638")),
                 "Length of 'col.per.group' should be of length 2 (the number of groups).",
                 fixed=T)
})


test_that("(plotIndiv:error): catches invalid pch values when using MINT", {
    
    data(stemcells)
    samples <- c(1:5,60:64)
    X <- stemcells$gene[samples, 1:10]
    Y <- rep(c("hESC", "hiPSC"),5)
    S <- as.character(stemcells$study[samples])
    
    choice.keepX <- c(3,3)
    
    res.mint.splsda <- mint.splsda(X, Y, study = S, keepX = choice.keepX)
    
    expect_error(plotIndiv(res.mint.splsda, pch = c(0,1,2)),
                 "'pch' needs to be of length 'object$study' as each of 'pch' represents a specific study",
                 fixed=T)
})


test_that("(plotIndiv:error): catches use of mint.block.(s)pls(da) as invalid", {
    
    data(breast.TCGA)
    mrna = rbind(breast.TCGA$data.train$mrna, breast.TCGA$data.test$mrna)
    mirna = rbind(breast.TCGA$data.train$mirna, breast.TCGA$data.test$mirna)
    X = list(mrna = mrna, mirna = mirna)
    Y = c(breast.TCGA$data.train$subtype, breast.TCGA$data.test$subtype)
    S = c(rep("study1",150), rep("study2",70))
    
    res.mint.block.plsda <- mint.block.plsda(X, Y, study = S)
    
    expect_error(plotIndiv(res.mint.block.plsda),
                 "No plotIndiv for the following functions at this stage: mint.block.pls, mint.block.spls, mint.block.plsda, mint.block.splsda.",
                 fixed=T)
})


###############################################################################
### =============================== WARNINGS ============================== ###
###############################################################################


test_that("(plotIndiv:warning): warned that col is ignored when group is supplied", {
    
    data(liver.toxicity)
    X <- liver.toxicity$gene[1:10,]
    Y <- liver.toxicity$clinic[1:10,]
    
    res.pls <- pls(X, Y)
    
    group <- c(rep("c1", 5), rep("c2", 5))
    
    expect_warning(plotIndiv(res.pls, group = group,
                             col.per.group = c("#388ECC", "#F68B33"),
                             col = c("#388ECC", "#F68B33")),
                 "'col' is ignored as 'group' has been set.",
                 fixed=T)
})


###############################################################################

dev.off()
unlink(list.files(pattern = "*.pdf"))

