
###############################################################################
### ============================ MISSING TESTS ============================ ###
###############################################################################

# basic
## mint.plsda
## mint.splsda

###############################################################################
### ============================ GROUND TRUTHS ============================ ###
###############################################################################

load(system.file("testdata", "testdata-background.predict.RData", package = "mixOmics"))
Testable.Components <- Test.Data$tc
Ground.Truths <- Test.Data$gt

###############################################################################
### ================================ BASIC ================================ ###
###############################################################################


test_that("(background.predict:basic): plsda", {
    
    testable.components <- Testable.Components$basic.plsda
    GT <- Ground.Truths$basic.plsda
    
    data(breast.tumors)
    X <- breast.tumors$gene.exp
    Y <- breast.tumors$sample$treatment
    
    res.plsda <- plsda(X, Y, ncomp = 2)
    
    plsda.bgp = background.predict(res.plsda, comp.predicted = 2, resolution = 10)
    
    invisible(capture.output(TT <- dput(plsda.bgp[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(background.predict:basic): splsda", {
    
    testable.components <- Testable.Components$basic.splsda
    GT <- Ground.Truths$basic.splsda

    data(breast.tumors)
    X <- breast.tumors$gene.exp
    Y <- breast.tumors$sample$treatment

    choice.keepX <- c(10, 10)

    res.splsda <- splsda(X, Y, ncomp = 2, keepX = choice.keepX)

    splsda.bgp = background.predict(res.splsda, comp.predicted = 2, resolution = 10)
    
    invisible(capture.output(TT <- dput(splsda.bgp[testable.components])))
    
    expect_equal(TT, GT)
})


###############################################################################
### ================================ DATA ================================= ###
###############################################################################


test_that("(background.predict:data): liver.toxicity", {
    
    testable.components <- Testable.Components$liver.toxicity.plsda
    GT <- Ground.Truths$liver.toxicity.plsda

    data(liver.toxicity)
    X = liver.toxicity$gene
    Y = as.factor(liver.toxicity$treatment[, 4])

    res.plsda <- plsda(X, Y, ncomp = 2)
    liver.toxicity.plsda.bgp = background.predict(res.plsda, comp.predicted = 2, resolution = 10)
    
    invisible(capture.output(TT <- dput(liver.toxicity.plsda.bgp[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(background.predict:data): srbct", {
    
    testable.components <- Testable.Components$srbct.plsda
    GT <- Ground.Truths$srbct.plsda

    data(srbct)
    X <- srbct$gene
    Y <- srbct$class

    res.plsda <- plsda(X, Y, ncomp = 2)
    srbct.plsda.bgp = background.predict(res.plsda, comp.predicted = 2, resolution = 10)
    
    invisible(capture.output(TT <- dput(srbct.plsda.bgp[testable.components])))
    
    expect_equal(TT, GT)
})


###############################################################################
### ============================= PARAMETER =============================== ###
###############################################################################


test_that("(background.predict:parameter): comp.predicted", {
    
    testable.components <- Testable.Components$comp.predicted.plsda
    GT <- Ground.Truths$comp.predicted.plsda

    data(breast.tumors)
    X <- breast.tumors$gene.exp
    Y <- breast.tumors$sample$treatment

    res.plsda <- plsda(X, Y, ncomp = 2)

    comp.predicted.bgp = background.predict(res.plsda, comp.predicted = 1, resolution = 10)
    
    invisible(capture.output(TT <- dput(comp.predicted.bgp[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(background.predict:parameter): dist", {
    
    testable.components <- Testable.Components$dist.plsda
    GT.max <- Ground.Truths$dist.max.plsda
    GT.centroids <- Ground.Truths$dist.centroids.plsda
    GT.mahalanobis <- Ground.Truths$dist.mahalanobis.plsda

    data(breast.tumors)
    X <- breast.tumors$gene.exp
    Y <- breast.tumors$sample$treatment

    res.plsda <- plsda(X, Y, ncomp = 2)

    max.dist.bgp = background.predict(res.plsda, comp.predicted = 2,
                                            dist = "max.dist", resolution = 10)

    centroids.dist.bgp = background.predict(res.plsda, comp.predicted = 2,
                                         dist = "centroids.dist", resolution = 10)

    mahalanobis.dist.bgp = background.predict(res.plsda, comp.predicted = 2,
                                         dist = "mahalanobis.dist", resolution = 10)
    
    invisible(capture.output(TT <- dput(max.dist.bgp[testable.components])))
    expect_equal(TT, GT.max)
    
    invisible(capture.output(TT <- dput(centroids.dist.bgp[testable.components])))
    expect_equal(TT, GT.centroids)
    
    invisible(capture.output(TT <- dput(mahalanobis.dist.bgp[testable.components])))
    expect_equal(TT, GT.mahalanobis)
})


test_that("(background.predict:parameter): resolution", {
    
    testable.components <- Testable.Components$resolution.plsda
    GT.20 <- Ground.Truths$resolution.20.plsda
    GT.30 <- Ground.Truths$resolution.30.plsda

    data(breast.tumors)
    X <- breast.tumors$gene.exp
    Y <- breast.tumors$sample$treatment

    res.plsda <- plsda(X, Y, ncomp = 2)

    res.20.bgp = background.predict(res.plsda, comp.predicted = 2,
                                       resolution = 20)

    res.30.bgp = background.predict(res.plsda, comp.predicted = 2,
                                       resolution = 30)
    
    invisible(capture.output(TT <- dput(res.20.bgp[testable.components])))
    expect_equal(TT, GT.20)
    
    invisible(capture.output(TT <- dput(res.30.bgp[testable.components])))
    expect_equal(TT, GT.30)
})


###############################################################################
### ================================ ERROR ================================ ###
###############################################################################



test_that("(background.predict:error): cannot use block.(s)plsda objects", {

    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.train$mirna,
             mRNA = breast.TCGA$data.train$mrna,
             proteomics = breast.TCGA$data.train$protein)
    Y = breast.TCGA$data.train$subtype

    res.block.plsda <- block.plsda(X, Y, design = "full")

    expect_error(background.predict(res.block.plsda),
                 "'background.predict' can only be calculated for 'plsda'
        and 'splsda' objects",
                 fixed = TRUE)

    choice.keepX <- list(miRNA=c(10,10),
                         mRNA=c(10,10),
                         proteomics=c(10,10))

    res.block.splsda <- block.splsda(X, Y, design = "full", keepX = choice.keepX)

    expect_error(background.predict(res.block.splsda),
                 "'background.predict' can only be calculated for 'plsda'
        and 'splsda' objects",
                 fixed = TRUE)
})


test_that("(background.predict:error): ensure dist has valid value", {

    data(breast.tumors)
    X <- breast.tumors$gene.exp
    Y <- breast.tumors$sample$treatment

    res.plsda <- plsda(X, Y, ncomp = 2)

    expect_error(background.predict(res.plsda, dist = "incorrect.dist"),
                 "Choose one of the three following distances: 'max.dist',
        'centroids.dist' or 'mahalanobis.dist'",
                 fixed = TRUE)

})


test_that("(background.predict:error): ensure dist has valid value", {

    data(breast.tumors)
    X <- breast.tumors$gene.exp
    Y <- breast.tumors$sample$treatment

    res.plsda <- plsda(X, Y, ncomp = 2)

    expect_error(background.predict(res.plsda, comp.predicted = 3),
                 "Can only show predicted background for 1 or 2 components",
                 fixed = TRUE)
})


