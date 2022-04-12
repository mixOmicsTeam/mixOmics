
###############################################################################
### ============================ MISSING TESTS ============================ ###
###############################################################################

# basic
## pca - once issue #171 is resolved we can do this one

###############################################################################
### ============================ GROUND TRUTHS ============================ ###
###############################################################################

load(system.file("testdata", "testdata-cim.RData", package = "mixOmics"))
Testable.Components <- Test.Data$tc
Ground.Truths <- Test.Data$gt

###############################################################################
### ================================ BASIC ================================ ###
###############################################################################


test_that("(cim:basic): matrices", {
    
    testable.components <- Testable.Components$basic.matrix
    GT <- Ground.Truths$basic.matrix
    
    data(nutrimouse)
    X <- nutrimouse$lipid[, 1:10]
    Y <- nutrimouse$gene[, 1:10]
    
    matrix.cim <- cim(cor(X, Y), cluster = "none")
    
    invisible(capture.output(TT <- dput(matrix.cim[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(cim:basic): spca", {
    
    testable.components <- Testable.Components$basic.spca
    GT <- Ground.Truths$basic.spca

    data(nutrimouse)
    X <- nutrimouse$lipid[1:10,1:10]

    res.spca <- spca(X, keepX = c(3,3))

    spca.cim <- cim(res.spca)
    
    invisible(capture.output(TT <- dput(spca.cim[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(cim:basic): ipca", {
    
    testable.components <- Testable.Components$basic.ipca
    GT <- Ground.Truths$basic.ipca

    data(nutrimouse)
    X <- nutrimouse$lipid[1:10,1:10]

    res.ipca <- ipca(X)

    ipca.cim <- cim(res.ipca)
    
    invisible(capture.output(TT <- dput(ipca.cim[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(cim:basic): sipca", {
    
    testable.components <- Testable.Components$basic.sipca
    GT <- Ground.Truths$basic.sipca

    data(nutrimouse)
    X <- nutrimouse$lipid[1:10,1:10]

    res.sipca <- sipca(X, ncomp=2, keepX = c(3,3))

    sipca.cim <- cim(res.sipca)
    
    invisible(capture.output(TT <- dput(sipca.cim[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(cim:basic): rcc", {
    
    testable.components <- Testable.Components$basic.rcc
    GT <- Ground.Truths$basic.rcc

    data(nutrimouse)
    X <- nutrimouse$lipid[,1:10]
    Y <- nutrimouse$gene[,1:10]

    res.rcc <- rcc(X, Y, method = "shrinkage")

    rcc.cim <- cim(res.rcc)
    
    invisible(capture.output(TT <- dput(rcc.cim[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(cim:basic): pls", {
    
    testable.components <- Testable.Components$basic.pls
    GT <- Ground.Truths$basic.pls

    data(liver.toxicity)
    X <- liver.toxicity$gene[,1:10]
    Y <- liver.toxicity$clinic[,1:10]

    res.pls <- pls(X, Y)

    pls.cim <- cim(res.pls)
    
    invisible(capture.output(TT <- dput(pls.cim[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(cim:basic): spls", {
    
    testable.components <- Testable.Components$basic.spls
    GT <- Ground.Truths$basic.spls

    data(liver.toxicity)
    X <- liver.toxicity$gene[,1:10]
    Y <- liver.toxicity$clinic[,1:10]

    res.spls <- spls(X, Y, keepX = c(3,3))

    spls.cim <- cim(res.spls)
    
    invisible(capture.output(TT <- dput(spls.cim[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(cim:basic): plsda", {
    
    testable.components <- Testable.Components$basic.plsda
    GT <- Ground.Truths$basic.plsda
    
    set.seed(16)
    samples <- sample(1:63, 10)

    data(srbct)
    X <- srbct$gene[samples,1:10]
    Y <- srbct$class[samples]

    res.plsda <- plsda(X, Y)

    plsda.cim <- cim(res.plsda)
    
    invisible(capture.output(TT <- dput(plsda.cim[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(cim:basic): splsda", {
    
    testable.components <- Testable.Components$basic.splsda
    GT <- Ground.Truths$basic.splsda

    set.seed(16)
    samples <- sample(1:63, 10)

    data(srbct)
    X <- srbct$gene[samples,1:10]
    Y <- srbct$class[samples]

    res.splsda <- splsda(X, Y, keepX = c(3,3))

    splsda.cim <- cim(res.splsda)
    
    invisible(capture.output(TT <- dput(splsda.cim[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(cim:basic): ml-spls", {
    
    testable.components <- Testable.Components$basic.mlspls
    GT <- Ground.Truths$basic.mlspls

    data(vac18)
    X <- vac18$genes[,1:10]
    Y <- vac18$genes[,11:20]
    ml <- vac18$sample

    res.mlpls <- spls(X, Y, keepX = c(3,3), keepY = c(3,3), multilevel = ml)

    mlspls.cim <- cim(res.mlpls)
    
    invisible(capture.output(TT <- dput(mlspls.cim[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(cim:basic): ml-splsda", {
    
    testable.components <- Testable.Components$basic.mlsplsda
    GT <- Ground.Truths$basic.mlsplsda
    
    data(vac18)
    X <- vac18$genes[,1:10]
    Y <- vac18$stimulation
    ml <- vac18$sample

    res.mlplsda <- splsda(X, Y, keepX = c(3,3), multilevel = ml)

    mlsplsda.cim <- cim(res.mlplsda)
    
    invisible(capture.output(TT <- dput(mlsplsda.cim[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(cim:basic): mint.spls", {
    
    testable.components <- Testable.Components$basic.mint.spls
    GT <- Ground.Truths$basic.mint.spls

    data(stemcells)
    X <- stemcells$gene[, 1:10]
    Y <- stemcells$gene[, 1:10]
    s <- stemcells$study

    res.mint.spls <- mint.spls(X, Y, study=s, keepX = c(3,3), keepY = c(3,3))

    mint.spls.cim <- cim(res.mint.spls)
    
    invisible(capture.output(TT <- dput(mint.spls.cim[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(cim:basic): mint.splsda", {
    
    testable.components <- Testable.Components$basic.mint.splsda
    GT <- Ground.Truths$basic.mint.splsda
    
    set.seed(17)
    samples <- sample(1:125, 40)

    data(stemcells)
    X <- stemcells$gene[samples, 1:10]
    Y <- stemcells$celltype[samples]
    s <- c(rep(1,10), rep(2,10), rep(3,10), rep(4,10))

    res.mint.splsda <- mint.splsda(X, Y, study=s, keepX = c(3,3))

    mint.splsda.cim <- cim(res.mint.splsda)
    
    invisible(capture.output(TT <- dput(mint.splsda.cim[testable.components])))
    
    expect_equal(TT, GT)
})


###############################################################################
### ================================ DATA ================================= ###
###############################################################################


test_that("(cim:data): multidrug", {
    
    testable.components <- Testable.Components$multidrug.spca
    GT <- Ground.Truths$multidrug.spca

    data(multidrug)
    X <- multidrug$ABC.trans[1:10, 1:10]

    res.multidrug <- spca(X, keepX = c(3,3))

    multidrug.cim <- cim(res.multidrug)
    
    invisible(capture.output(TT <- dput(multidrug.cim[testable.components])))
    
    expect_equal(TT, GT)
})


###############################################################################
### ============================= PARAMETER =============================== ###
###############################################################################


test_that("(cim:paramter): cutoff", {
    
    testable.components <- Testable.Components$cutoff.spls
    GT <- Ground.Truths$cutoff.spls

    data(liver.toxicity)
    X <- liver.toxicity$gene[,1:10]
    Y <- liver.toxicity$clinic[,1:10]

    res.spls <- spls(X, Y, keepX = c(3,3), keepY = c(3,3))

    cutoff.cim <- cim(res.spls, cutoff = 0.2)
    
    invisible(capture.output(TT <- dput(cutoff.cim[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(cim:paramter): cut.tree", {
    
    testable.components <- Testable.Components$cut.tree.spls
    GT <- Ground.Truths$cut.tree.spls

    data(liver.toxicity)
    X <- liver.toxicity$gene[,1:10]
    Y <- liver.toxicity$clinic[,1:10]

    res.spls <- spls(X, Y, keepX = c(3,3), keepY = c(3,3))

    cut.tree.cim <- cim(res.spls, cut.tree = c(0.5, 0.5))
    
    invisible(capture.output(TT <- dput(cut.tree.cim[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(cim:paramter): cluster", {
    
    testable.components <- Testable.Components$cluster.spls
    GT.none <- Ground.Truths$cluster.none.spls
    GT.row <- Ground.Truths$cluster.row.spls

    data(liver.toxicity)
    X <- liver.toxicity$gene[,1:10]
    Y <- liver.toxicity$clinic[,1:10]

    res.spls <- spls(X, Y, keepX = c(3,3), keepY = c(3,3))

    cluster.none.cim <- cim(res.spls, cluster = "none")
    cluster.row.cim <- cim(res.spls, cluster = "row")
    
    invisible(capture.output(TT <- dput(cluster.none.cim[testable.components])))
    expect_equal(TT, GT.none)
    
    invisible(capture.output(TT <- dput(cluster.row.cim[testable.components])))
    expect_equal(TT, GT.row)
})


test_that("(cim:paramter): dist.method", {
    
    testable.components <- Testable.Components$dist.method.spls
    GT.correlation <- Ground.Truths$dist.correlation.spls
    GT.manhattan <- Ground.Truths$dist.manhattan.spls

    data(liver.toxicity)
    X <- liver.toxicity$gene[,1:10]
    Y <- liver.toxicity$clinic[,1:10]

    res.spls <- spls(X, Y, keepX = c(3,3), keepY = c(3,3))

    dist.correlation.cim <- cim(res.spls, dist.method = c("correlation", "correlation"),
                                cluster = "row")
    dist.manhattan.cim <- cim(res.spls, dist.method = c("manhattan", "manhattan"),
                              cluster = "row")
    
    invisible(capture.output(TT <- dput(dist.correlation.cim[testable.components])))
    expect_equal(TT, GT.correlation)
    
    invisible(capture.output(TT <- dput(dist.manhattan.cim[testable.components])))
    expect_equal(TT, GT.manhattan)
})


test_that("(cim:paramter): clust.method", {
    
    testable.components <- Testable.Components$clust.method.spls
    GT.centroid <- Ground.Truths$clust.centroid.spls
    GT.complete <- Ground.Truths$clust.complete.spls

    data(liver.toxicity)
    X <- liver.toxicity$gene[,1:10]
    Y <- liver.toxicity$clinic[,1:10]

    res.spls <- spls(X, Y, keepX = c(3,3), keepY = c(3,3))

    clust.centroid.cim <- cim(res.spls, clust.method = c("centroid", "centroid"),
                                cluster = "row")
    clust.complete.cim <- cim(res.spls, clust.method = c("complete", "complete"),
                                cluster = "row")
    
    invisible(capture.output(TT <- dput(clust.centroid.cim[testable.components])))
    expect_equal(TT, GT.centroid)
    
    invisible(capture.output(TT <- dput(clust.complete.cim[testable.components])))
    expect_equal(TT, GT.complete)
})


test_that("(cim:paramter): comp",{
    
    testable.components <- Testable.Components$comp.spls
    GT <- Ground.Truths$comp.spls

    data(liver.toxicity)
    X <- liver.toxicity$gene[,1:10]
    Y <- liver.toxicity$clinic[,1:10]

    res.spls <- spls(X, Y, keepX = c(3,3), keepY = c(3,3), ncomp = 4)

    comp.cim <- cim(res.spls, comp = c(2,3))
    
    invisible(capture.output(TT <- dput(comp.cim[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(cim:paramter): mapping",{
    
    testable.components <- Testable.Components$mapping.spls
    GT <- Ground.Truths$mapping.spls

    data(liver.toxicity)
    X <- liver.toxicity$gene[1:10,1:10]
    Y <- liver.toxicity$clinic[1:10,1:10]

    res.spls <- spls(X, Y, keepX = c(3,3), keepY = c(3,3))

    mapping.cim <- cim(res.spls, mapping = "X")
    
    invisible(capture.output(TT <- dput(mapping.cim[testable.components])))
    
    expect_equal(TT, GT)
})


###############################################################################
### ================================ ERROR ================================ ###
###############################################################################


test_that("(cim:error): cluster is valid value", {
    
    data(nutrimouse)
    X <- nutrimouse$lipid
    
    res.spca <- spca(X)
    
    expect_error(cim(res.spca, cluster = "invalid.cluster"),
                 "'cluster' should be one of 'both', 'row', 'column' or 'none'.",
                 fixed=TRUE)
})


test_that("(cim:error): clust.method is valid value", {
    
    data(nutrimouse)
    X <- nutrimouse$lipid
    
    res.spca <- spca(X)
    
    expect_error(cim(res.spca, clust.method = c("ward", "invalid.cluster")),
                 "invalid clustering method.",
                 fixed=TRUE)
})


test_that("(cim:error): dist.method is valid value", {
    
    data(nutrimouse)
    X <- nutrimouse$lipid
    
    res.spca <- spca(X)
    
    expect_error(cim(res.spca, dist.method = c("euclidean", "invalid.cluster")),
                 "invalid distance method.",
                 fixed=TRUE)
})


test_that("(cim:error): cut.tree is valid value", {
    
    data(nutrimouse)
    X <- nutrimouse$lipid
    
    res.spca <- spca(X)
    
    expect_error(cim(res.spca, cut.tree=0.8),
                 "'cut.tree' must be a numeric vector of length 2.",
                 fixed=TRUE)
    
    expect_error(cim(res.spca, cut.tree=c(2,2)),
                 "Components of 'cut.tree' must be between 0 and 1.",
                 fixed=TRUE)
})


test_that("(cim:error): doesn't allow block.splsda objects", {
    
    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.train$mirna,
             mRNA = breast.TCGA$data.train$mrna,
             proteomics = breast.TCGA$data.train$protein)
    Y = breast.TCGA$data.train$subtype
    
    res.block.splsda <- block.splsda(X, Y)
    
    expect_error(cim(res.block.splsda),
                 "Please call the 'cimDiablo' function on your 'block.splsda' object",
                 fixed=TRUE)
})