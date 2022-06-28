
###############################################################################
### ============================ MISSING TESTS ============================ ###
###############################################################################

# basic
## mint.pls - add down the line once PR #212 is merged
## mint.spls - add down the line once PR #212 is merged

# parameter
## block - add down the line once PR #212 is merged

# error
## "Duplicate in 'study' not allowed" - add down the line once PR #212 is merged
## "'study' must be one of 'object$study' or 'all'." - add down the line once PR #212 is merged

# edge cases
## 

###############################################################################
### ============================ GROUND TRUTHS ============================ ###
###############################################################################

Test.Data <- readRDS(system.file("testdata", "testdata-plotLoadings.rda", package = "mixOmics"))
Ground.Truths <- Test.Data$gt

###############################################################################
### ================================ BASIC ================================ ###
###############################################################################


test_that("(plotLoadings:basic): pca", {
    
    GT <- Ground.Truths$basic.pca
    
    data(multidrug)
    X <- multidrug$ABC.trans[,1:10]
    
    res.pca <- pca(X)
    
    pca.plotLoadings = plotLoadings(res.pca)
    
    invisible(capture.output(TT <- dput(pca.plotLoadings)))
    
    expect_equal(TT, GT)
})


test_that("(plotLoadings:basic): spca", {
    
    GT <- Ground.Truths$basic.spca
    
    data(multidrug)
    X <- multidrug$ABC.trans[,1:10]
    
    choice.keepX <- c(3,3)
    
    res.spca <- spca(X, keepX = choice.keepX)
    
    spca.plotLoadings = plotLoadings(res.spca)
    
    invisible(capture.output(TT <- dput(spca.plotLoadings)))
    
    expect_equal(TT, GT)
})


test_that("(plotLoadings:basic): pls", {
    
    GT <- Ground.Truths$basic.pls
    
    data(liver.toxicity)
    X <- liver.toxicity$gene[,1:10]
    Y <- liver.toxicity$clinic[,1:10]
    
    res.pls <- pls(X, Y)
    
    pls.plotLoadings = plotLoadings(res.pls)
    
    invisible(capture.output(TT <- dput(pls.plotLoadings)))
    
    expect_equal(TT, GT)
})


test_that("(plotLoadings:basic): spls", {
    
    GT <- Ground.Truths$basic.spls
    
    data(liver.toxicity)
    X <- liver.toxicity$gene[,1:10]
    Y <- liver.toxicity$clinic[,1:10]
    
    choice.keepX <- c(3,3)
    choice.keepY <- c(3,3)
    
    res.spls <- spls(X, Y, keepX = choice.keepX, keepY = choice.keepY)
    
    spls.plotLoadings = plotLoadings(res.spls)
    
    invisible(capture.output(TT <- dput(spls.plotLoadings)))
    
    expect_equal(TT, GT)
})


test_that("(plotLoadings:basic): plsda", {
    
    GT <- Ground.Truths$basic.plsda
    
    data(srbct)
    X <- srbct$gene[,1:10]
    Y <- srbct$class
    
    res.plsda <- plsda(X, Y)
    
    plsda.plotLoadings = plotLoadings(res.plsda)
    
    invisible(capture.output(TT <- dput(plsda.plotLoadings)))
    
    expect_equal(TT, GT)
})


test_that("(plotLoadings:basic): splsda", {
    
    GT <- Ground.Truths$basic.splsda
    
    data(srbct)
    X <- srbct$gene[,1:10]
    Y <- srbct$class
    
    choice.keepX <- c(3,3)
    
    res.splsda <- splsda(X, Y, keepX = choice.keepX)
    
    splsda.plotLoadings = plotLoadings(res.splsda)
    
    invisible(capture.output(TT <- dput(splsda.plotLoadings)))
    
    expect_equal(TT, GT)
})


test_that("(plotLoadings:basic): block.pls", {
    
    GT <- Ground.Truths$basic.block.pls
    
    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.train$mirna[,1:10],
             mRNA = breast.TCGA$data.train$mrna[,1:10],
             proteomics = breast.TCGA$data.train$protein[,1:10])
    
    res.block.pls <- block.pls(X, indY = 3)
    
    block.pls.plotLoadings = plotLoadings(res.block.pls)
    
    invisible(capture.output(TT <- dput(block.pls.plotLoadings)))
    
    expect_equal(TT, GT)
})


test_that("(plotLoadings:basic): block.spls", {
    
    GT <- Ground.Truths$basic.block.spls
    
    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.train$mirna[,1:10],
             mRNA = breast.TCGA$data.train$mrna[,1:10],
             proteomics = breast.TCGA$data.train$protein[,1:10])
    
    choice.keepX = list(miRNA = c(3,3),
                       mRNA = c(3,3),
                       proteomics = c(3,3))
    
    res.block.spls <- block.spls(X, indY = 3, keepX = choice.keepX)
    
    block.spls.plotLoadings = plotLoadings(res.block.spls)
    
    invisible(capture.output(TT <- dput(block.spls.plotLoadings)))
    
    expect_equal(TT, GT)
})


test_that("(plotLoadings:basic): block.plsda", {
    
    GT <- Ground.Truths$basic.block.plsda
    
    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.train$mirna[,1:10],
             mRNA = breast.TCGA$data.train$mrna[,1:10],
             proteomics = breast.TCGA$data.train$protein[,1:10])
    Y = breast.TCGA$data.train$subtype
    
    res.block.plsda <- block.plsda(X, Y)
    
    block.plsda.plotLoadings = plotLoadings(res.block.plsda)
    
    invisible(capture.output(TT <- dput(block.plsda.plotLoadings)))
    
    expect_equal(TT, GT)
})


test_that("(plotLoadings:basic): block.splsda", {
    
    GT <- Ground.Truths$basic.block.splsda
    
    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.train$mirna[,1:10],
             mRNA = breast.TCGA$data.train$mrna[,1:10],
             proteomics = breast.TCGA$data.train$protein[,1:10])
    Y = breast.TCGA$data.train$subtype
    
    choice.keepX = list(miRNA = c(3,3),
                       mRNA = c(3,3),
                       proteomics = c(3,3))
    
    res.block.splsda <- block.splsda(X, Y, keepX = choice.keepX)
    
    block.splsda.plotLoadings = plotLoadings(res.block.splsda)
    
    invisible(capture.output(TT <- dput(block.splsda.plotLoadings)))
    
    expect_equal(TT, GT)
})


test_that("(plotLoadings:basic): mint.plsda", {

    GT <- Ground.Truths$basic.mint.plsda

    data(stemcells)
    samples <- c(1:5,60:64)
    X <- stemcells$gene[samples, 1:10]
    Y <- rep(c("hESC", "hiPSC"),5)
    S <- as.character(stemcells$study[samples])

    res.mint.plsda <- mint.plsda(X, Y, study = S)

    mint.plsda.plotLoadings = plotLoadings(res.mint.plsda)

    invisible(capture.output(TT <- dput(mint.plsda.plotLoadings)))

    expect_equal(TT, GT)
})


test_that("(plotLoadings:basic): mint.splsda", {

    GT <- Ground.Truths$basic.mint.splsda

    samples <- c(1:5,60:64)
    X <- stemcells$gene[samples, 1:10]
    Y <- rep(c("hESC", "hiPSC"),5)
    S <- as.character(stemcells$study[samples])

    choice.keepX = c(3,3)

    res.mint.splsda <- mint.splsda(X, Y, study = S, keepX = choice.keepX)

    mint.splsda.plotLoadings = plotLoadings(res.mint.splsda)

    invisible(capture.output(TT <- dput(mint.splsda.plotLoadings)))

    expect_equal(TT, GT)
})


###############################################################################
### ============================== PARAMETER ============================== ###
###############################################################################


test_that("(plotLoadings:parameter): comp", {
    
    GT <- Ground.Truths$comp.pca
    
    data(multidrug)
    X <- multidrug$ABC.trans[,1:10]
    
    res.pca <- pca(X)
    
    comp.plotLoadings = plotLoadings(res.pca, comp = 2)
    
    invisible(capture.output(TT <- dput(comp.plotLoadings)))
    
    expect_equal(TT, GT)
})


test_that("(plotLoadings:parameter): ndisplay", {
    
    # test ndisplay using single omics type
    
    GT <- Ground.Truths$ndisplay.pca
    
    data(multidrug)
    X <- multidrug$ABC.trans[,1:10]
    
    res.pca <- pca(X)
    
    ndisplay.plotLoadings = plotLoadings(res.pca, ndisplay = 3)
    
    invisible(capture.output(TT <- dput(ndisplay.plotLoadings)))
    
    expect_equal(TT, GT)
    
    # ------------------------------------------------------------- #
    # test ndisplay using multi omics type
    
    GT <- Ground.Truths$ndisplay.pls
    
    data(liver.toxicity)
    X <- liver.toxicity$gene[,1:10]
    Y <- liver.toxicity$clinic[,1:10]
    
    res.pls <- pls(X, Y)
    
    ndisplay.plotLoadings = plotLoadings(res.pls, ndisplay = 3)
    
    invisible(capture.output(TT <- dput(ndisplay.plotLoadings)))
    
    expect_equal(TT, GT)
})


test_that("(plotLoadings:parameter): contrib/method", {
    
    GT <- Ground.Truths$contrib.method.plsda
    
    data(srbct)
    X <- srbct$gene[,1:5]
    Y <- srbct$class
    
    res.plsda <- plsda(X, Y)
    
    contrib.method.plotLoadings = plotLoadings(res.plsda, contrib = "max", 
                                               method = "median")
    
    invisible(capture.output(TT <- dput(contrib.method.plotLoadings)))
    
    expect_equal(TT, GT)
})


test_that("(plotLoadings:parameter): study", {

    GT <- Ground.Truths$study.mint.plsda

    data(stemcells)
    samples <- c(1:5,60:64)
    X <- stemcells$gene[samples, 1:10]
    Y <- rep(c("hESC", "hiPSC"),5)
    S <- as.character(stemcells$study[samples])

    res.mint.plsda <- mint.plsda(X, Y, study = S)

    study.plotLoadings = plotLoadings(res.mint.plsda, study = "1")

    invisible(capture.output(TT <- dput(study.plotLoadings)))

    expect_equal(TT, GT)
})


test_that("(plotLoadings:parameter): name.var", {

    GT <- Ground.Truths$name.var.block.plsda
    
    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.train$mirna[,1:10],
             mRNA = breast.TCGA$data.train$mrna[,1:10],
             proteomics = breast.TCGA$data.train$protein[,1:10])
    Y = breast.TCGA$data.train$subtype
    
    res.block.plsda <- block.plsda(X, Y)
    
    n <- list(miRNA = rep("miRNA.names", 10), 
              mRNA = rep("mRNA.names", 10), 
              proteomics = rep("protein.names", 10))
    
    name.var.plotLoadings = plotLoadings(res.block.plsda, name.var = n)

    invisible(capture.output(TT <- dput(name.var.plotLoadings)))

    expect_equal(TT, GT)
})


###############################################################################
### ================================ ERROR ================================ ###
###############################################################################



test_that("(plotLoadings:error): catches object with no loading values", {
    
    data(multidrug)
    X <- multidrug$ABC.trans[,1:10]
    
    res.pca <- pca(X)
    
    res.pca$loadings <- NULL
    
    expect_error(plotLoadings(res.pca),
                 "'plotLoadings' should be used on object for which object$loadings is present.",
                 fixed=T)
})


test_that("(plotLoadings:error): catches invalid `block` values", {
    
    data(srbct)
    X <- srbct$gene[,1:10]
    Y <- srbct$class
    
    res.plsda <- plsda(X, Y)
    
    expect_error(plotLoadings(res.plsda, block = 2),
                 "'block' can only be 'X' or '1' for plsda and splsda object",
                 fixed=T)
    
    
    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.train$mirna[,1:10],
             mRNA = breast.TCGA$data.train$mrna[,1:10],
             proteomics = breast.TCGA$data.train$protein[,1:10])
    Y = breast.TCGA$data.train$subtype
    
    res.block.plsda <- block.plsda(X, Y)
    
    expect_error(plotLoadings(res.block.plsda, block = 4),
                 "'block' needs to be lower than the number of blocks in the fitted model, which is 3",
                 fixed=T)
    
    expect_error(plotLoadings(res.block.plsda, block = "random.block"),
                 "Incorrect value for 'block', 'block' should be among the blocks used in your object: miRNA, mRNA, proteomics",
                 fixed=T)
})


test_that("(plotLoadings:error): catches invalid `contrib` values", {
    
    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.train$mirna[,1:10],
             mRNA = breast.TCGA$data.train$mrna[,1:10],
             proteomics = breast.TCGA$data.train$protein[,1:10])
    Y = breast.TCGA$data.train$subtype
    
    res.block.plsda <- block.plsda(X, Y)
    
    expect_error(plotLoadings(res.block.plsda, contrib = "random.value"),
                 "'contrib' must be either 'min' or 'max'",
                 fixed=T)
})


test_that("(plotLoadings:error): catches invalid `name.var` values", {
    
    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.train$mirna[,1:10],
             mRNA = breast.TCGA$data.train$mrna[,1:10],
             proteomics = breast.TCGA$data.train$protein[,1:10])
    Y = breast.TCGA$data.train$subtype
    
    res.block.plsda <- block.plsda(X, Y)
    
    expect_error(plotLoadings(res.block.plsda, name.var = list(miRNA = rep("miRNA.names", 10),
                                                               mRNA = rep("mRNA.names", 10),
                                                               proteomics = rep("protein.names", 9))),
                 "For block 'proteomics', 'name.var' should be a vector of length 10",
                 fixed=T)
    
    expect_error(plotLoadings(res.block.plsda, name.var = list(miRNA = rep("miRNA.names", 10),
                                                               mRNA = rep("mRNA.names", 10))),
                 "'names' has to be a list of length the number of block to plot: 3",
                 fixed=T)
})


# test_that("(plotLoadings:error): catches if study has duplicate entries", {
#     
#     data(stemcells)
#     samples <- c(1:5,60:64)
#     X <- stemcells$gene[samples, 1:10]
#     Y <- rep(c("hESC", "hiPSC"),5)
#     S <- as.character(stemcells$study[samples])
# 
#     res.mint.plsda <- mint.plsda(X, Y, study = S)
#     
#     expect_error(plotLoadings(res.mint.plsda, study = c("1", "1")),
#                  "Duplicate in 'study' not allowed",
#                  fixed=T)
# })


###############################################################################
### =============================== WARNINGS ============================== ###
###############################################################################





###############################################################################

dev.off()
unlink(list.files(pattern = "*.pdf"))


# 
# test_that("plotLoadings margin errrors is handled properly", code = {
#     data(nutrimouse)
#     Y = nutrimouse$diet
#     gene = nutrimouse$gene
#     lipid = nutrimouse$lipid
#     ## extend feature names
#     suff <- "-a-long-suffix-from-abolutely-nowhere-which-is-gonna-be-longer-than-margins"
#     colnames(gene) <- paste0(colnames(gene), suff)
#     colnames(lipid) <- paste0(colnames(lipid), suff)
#     data = list(gene = gene, lipid = lipid)
#     design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
#     
#     nutrimouse.sgccda = block.splsda(X = data,
#                                      Y = Y,
#                                      design = design,
#                                      keepX = list(gene = c(10,10), lipid = c(15,15)),
#                                      ncomp = 2,
#                                      scheme = "centroid")
#     expect_error(plotLoadings(nutrimouse.sgccda, contrib = "min"), regexp = "plotLoadings encountered margin errors")
# })
