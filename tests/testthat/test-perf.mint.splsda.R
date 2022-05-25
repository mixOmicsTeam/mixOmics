
###############################################################################
### ============================ MISSING TESTS ============================ ###
###############################################################################


###############################################################################
### ============================ GROUND TRUTHS ============================ ###
###############################################################################

Test.Data <- readRDS(system.file("testdata", "testdata-perf.mint.splsda.rda", package = "mixOmics"))
Testable.Components <- Test.Data$tc
Ground.Truths <- Test.Data$gt

###############################################################################
### ================================ BASIC ================================ ###
###############################################################################


test_that("(perf.mint.splsda:basic): mint.plsda", {
    
    testable.components <- Testable.Components$basic.mint.plsda
    GT <- Ground.Truths$basic.mint.plsda
    
    data(stemcells)
    samples <- c(1:5,60:64)
    X <- stemcells$gene[samples, 1:10]
    Y <- rep(c("hESC", "hiPSC"),5)
    S <- as.character(stemcells$study[samples])

    res.mint.plsda <- mint.plsda(X, Y, study=S)
    
    mint.plsda.perf <- perf(res.mint.plsda)
    
    invisible(capture.output(TT <- dput(mint.plsda.perf[testable.components])))

    expect_equal(TT, GT)
})


test_that("(perf.mint.splsda:basic): mint.splsda", {
    
    testable.components <- Testable.Components$basic.mint.splsda
    GT <- Ground.Truths$basic.mint.splsda
    
    data(stemcells)
    samples <- c(1:5,60:64)
    X <- stemcells$gene[samples, 1:10]
    Y <- rep(c("hESC", "hiPSC"),5)
    S <- as.character(stemcells$study[samples])
    
    choice.keepX <- c(3,3)

    res.mint.splsda <- mint.splsda(X, Y, study=S, keepX = choice.keepX)
    
    mint.splsda.perf <- perf(res.mint.splsda)
    
    invisible(capture.output(TT <- dput(mint.splsda.perf[testable.components])))

    expect_equal(TT, GT)
})


###############################################################################
### ============================== PARAMETER ============================== ###
###############################################################################


test_that("(perf.diablo:parameter): dist", {

    testable.components <- Testable.Components$dist.mint.splsda

    data(stemcells)
    samples <- c(1:5,60:64)
    X <- stemcells$gene[samples, 1:10]
    Y <- rep(c("hESC", "hiPSC"),5)
    S <- as.character(stemcells$study[samples])
    
    choice.keepX <- c(3,3)

    res.mint.splsda <- mint.splsda(X, Y, study=S, keepX = choice.keepX)

    # max.dist
    GT <- Ground.Truths$dist.max.block.splsda
    set.seed(101)
    mint.splsda.perf = perf(res.mint.splsda, folds=3, dist = "max.dist")
    invisible(capture.output(TT <- dput(mint.splsda.perf[testable.components])))
    expect_equal(TT, GT)

    # centroids.dist
    GT <- Ground.Truths$dist.centroids.block.splsda
    set.seed(101)
    mint.splsda.perf = perf(res.mint.splsda, folds=3, dist = "centroids.dist")
    invisible(capture.output(TT <- dput(mint.splsda.perf[testable.components])))
    expect_equal(TT, GT)

    # mahalanobis.dist
    GT <- Ground.Truths$dist.mahalanobis.block.splsda
    set.seed(101)
    mint.splsda.perf = perf(res.mint.splsda, folds=3, dist = "mahalanobis.dist")
    invisible(capture.output(TT <- dput(mint.splsda.perf[testable.components])))
    expect_equal(TT, GT)
})


test_that("(perf.diablo:parameter): auc", {

    testable.components <- Testable.Components$auc.mint.splsda
    GT <- Ground.Truths$auc.mint.splsda

    data(stemcells)
    samples <- c(1:5,60:64)
    X <- stemcells$gene[samples, 1:10]
    Y <- rep(c("hESC", "hiPSC"),5)
    S <- as.character(stemcells$study[samples])
    
    choice.keepX <- c(3,3)

    res.mint.splsda <- mint.splsda(X, Y, study=S, keepX = choice.keepX)

    set.seed(100)
    mint.splsda.perf = perf(res.mint.splsda, folds=3, auc = T)

    invisible(capture.output(TT <- dput(mint.splsda.perf[testable.components])))

    expect_equal(TT, GT)
})