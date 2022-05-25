
###############################################################################
### ============================ MISSING TESTS ============================ ###
###############################################################################

# edge cases
## At least one class is not represented in one fold, which may unbalance the error rate.\n  Consider a number of folds lower than the minimum in table(Y):

###############################################################################
### ============================ GROUND TRUTHS ============================ ###
###############################################################################

Test.Data <- readRDS(system.file("testdata", "testdata-perf.diablo.rda", package = "mixOmics"))
Testable.Components <- Test.Data$tc
Ground.Truths <- Test.Data$gt

###############################################################################
### ================================ BASIC ================================ ###
###############################################################################


test_that("(perf.diablo:basic): block.plsda", {
    
    testable.components <- Testable.Components$basic.block.plsda
    GT <- Ground.Truths$basic.block.plsda
    
    data(breast.TCGA)
    samples <- c(1:4, 46:49, 79:82)
    X = list(miRNA = breast.TCGA$data.train$mirna[samples, 1:10],
            mRNA = breast.TCGA$data.train$mrna[samples, 1:10])
    Y = breast.TCGA$data.train$subtype[samples]
    
    res.block.plsda <- block.plsda(X, Y)
    
    set.seed(100)
    block.plsda.perf = perf(res.block.plsda, folds=3)
    
    invisible(capture.output(TT <- dput(block.plsda.perf[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(perf.diablo:basic): block.splsda", {
    
    testable.components <- Testable.Components$basic.block.splsda
    GT <- Ground.Truths$basic.block.splsda
    
    data(breast.TCGA)
    samples <- c(1:4, 46:49, 79:82)
    X = list(miRNA = breast.TCGA$data.train$mirna[samples, 1:10],
            mRNA = breast.TCGA$data.train$mrna[samples, 1:10])
    Y = breast.TCGA$data.train$subtype[samples]
    
    keepX <- list(miRNA= c(5,5),
                  mRNA = c(5,5))
    
    res.block.splsda <- block.splsda(X, Y, keepX = keepX)
    
    set.seed(100)
    block.splsda.perf = perf(res.block.splsda, folds=3)
    
    invisible(capture.output(TT <- dput(block.splsda.perf[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(perf.diablo:basic): sgccda", {
    
    testable.components <- Testable.Components$basic.block.splsda
    GT <- Ground.Truths$basic.block.splsda
    
    data(breast.TCGA)
    samples <- c(1:4, 46:49, 79:82)
    X = list(miRNA = breast.TCGA$data.train$mirna[samples, 1:10],
            mRNA = breast.TCGA$data.train$mrna[samples, 1:10])
    Y = breast.TCGA$data.train$subtype[samples]
    
    keepX <- list(miRNA= c(5,5),
                  mRNA = c(5,5))
    
    res.sgccda <- wrapper.sgccda(X, Y, keepX = keepX)
    
    set.seed(100)
    sgccda.perf = perf(res.sgccda, folds=3)
    
    invisible(capture.output(TT <- dput(sgccda.perf[testable.components])))
    
    expect_equal(TT, GT)
})


###############################################################################
### ================================ DATA ================================= ###
###############################################################################


test_that("(perf.diablo:data): breast.tcga test", {
    
    testable.components <- Testable.Components$breast.test.block.splsda
    GT <- Ground.Truths$breast.test.block.splsda
    
    data(breast.TCGA)
    samples <- c(1:4, 22:26, 35:39)
    X = list(miRNA = breast.TCGA$data.test$mirna[samples, 1:10],
            mRNA = breast.TCGA$data.test$mrna[samples, 1:10])
    Y = breast.TCGA$data.test$subtype[samples]
    
    keepX <- list(miRNA= c(5,5),
                  mRNA = c(5,5))
    
    res.block.splsda <- block.splsda(X, Y, keepX = keepX)
    
    set.seed(100)
    block.splsda.perf = perf(res.block.splsda, folds=3)
    
    invisible(capture.output(TT <- dput(block.splsda.perf[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(perf.diablo:data): block.splsda", {
    
    testable.components <- Testable.Components$breast.test.block.splsda
    GT <- Ground.Truths$breast.test.block.splsda
    
    data(breast.TCGA)
    samples <- c(1:4, 22:26, 35:39)
    X = list(miRNA = breast.TCGA$data.test$mirna[samples, 1:10],
            mRNA = breast.TCGA$data.test$mrna[samples, 1:10])
    Y = breast.TCGA$data.test$subtype[samples]
    
    keepX <- list(miRNA= c(5,5),
                  mRNA = c(5,5))
    
    res.block.splsda <- block.splsda(X, Y, keepX = keepX)
    
    set.seed(100)
    block.splsda.perf = perf(res.block.splsda, folds=3)
    
    invisible(capture.output(TT <- dput(block.splsda.perf[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(perf.diablo:data): nutrimouse", {
    
    testable.components <- Testable.Components$nutrimouse.block.splsda
    GT <- Ground.Truths$nutrimouse.block.splsda
    
    data(nutrimouse)
    samples <- c(6, 10, 17, 2, 3, 13, 4, 9, 11, 1, 7, 8, 5, 12, 14)
    X = list(lipid = nutrimouse$lipid[samples, 1:10],
             gene = nutrimouse$gene[samples, 1:10])
    Y = nutrimouse$diet[samples]
    
    keepX <- list(lipid= c(5,5),
                  gene = c(5,5))
    
    res.block.splsda <- block.splsda(X, Y, keepX = keepX)
    
    set.seed(101)
    block.splsda.perf = perf(res.block.splsda, folds=3)
    
    invisible(capture.output(TT <- dput(block.splsda.perf[testable.components])))
    
    expect_equal(TT, GT)
})


###############################################################################
### ============================== PARAMETER ============================== ###
###############################################################################


test_that("(perf.diablo:parameter): dist", {
    
    testable.components <- Testable.Components$dist.block.splsda
    
    data(breast.TCGA)
    samples <- c(1:4, 46:49, 79:82)
    X = list(miRNA = breast.TCGA$data.train$mirna[samples, 1:10],
             mRNA = breast.TCGA$data.train$mrna[samples, 1:10])
    Y = breast.TCGA$data.train$subtype[samples]
    
    keepX <- list(miRNA= c(5,5),
                  mRNA = c(5,5))
    
    res.block.splsda <- block.splsda(X, Y, keepX = keepX)
    
    # max.dist
    GT <- Ground.Truths$dist.max.block.splsda
    set.seed(101)
    block.splsda.perf = perf(res.block.splsda, folds=3, dist = "max.dist")
    invisible(capture.output(TT <- dput(block.splsda.perf[testable.components])))
    expect_equal(TT, GT)
    
    # centroids.dist
    GT <- Ground.Truths$dist.centroids.block.splsda
    set.seed(101)
    block.splsda.perf = perf(res.block.splsda, folds=3, dist = "centroids.dist")
    invisible(capture.output(TT <- dput(block.splsda.perf[testable.components])))
    expect_equal(TT, GT)
    
    # mahalanobis.dist
    GT <- Ground.Truths$dist.mahalanobis.block.splsda
    set.seed(101)
    block.splsda.perf = perf(res.block.splsda, folds=3, dist = "mahalanobis.dist")
    invisible(capture.output(TT <- dput(block.splsda.perf[testable.components])))
    expect_equal(TT, GT)
})


test_that("(perf.diablo:parameter): validation", {
    
    testable.components <- Testable.Components$validation.block.splsda
    GT <- Ground.Truths$validation.block.splsda
    
    data(breast.TCGA)
    samples <- c(1:4, 46:49, 79:82)
    X = list(miRNA = breast.TCGA$data.train$mirna[samples, 1:10],
             mRNA = breast.TCGA$data.train$mrna[samples, 1:10])
    Y = breast.TCGA$data.train$subtype[samples]
    
    keepX <- list(miRNA= c(5,5),
                  mRNA = c(5,5))
    
    res.block.splsda <- block.splsda(X, Y, keepX = keepX)
    
    set.seed(101)
    block.splsda.perf = perf(res.block.splsda, validation = "loo")
    
    invisible(capture.output(TT <- dput(block.splsda.perf[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(perf.diablo:parameter): folds", {
    
    testable.components <- Testable.Components$folds.block.splsda
    GT <- Ground.Truths$folds.block.splsda
    
    data(breast.TCGA)
    samples <- c(1:5, 46:50, 79:83)
    X = list(miRNA = breast.TCGA$data.train$mirna[samples, 1:10],
            mRNA = breast.TCGA$data.train$mrna[samples, 1:10])
    Y = breast.TCGA$data.train$subtype[samples]
    
    keepX <- list(miRNA= c(5,5),
                  mRNA = c(5,5))
    
    res.block.splsda <- block.splsda(X, Y, keepX = keepX)
    
    set.seed(100)
    block.splsda.perf = perf(res.block.splsda, folds=5)
    
    invisible(capture.output(TT <- dput(block.splsda.perf[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(perf.diablo:parameter): nrepeat", {
    
    testable.components <- Testable.Components$nrepeat.block.splsda
    GT <- Ground.Truths$nrepeat.block.splsda
    
    data(breast.TCGA)
    samples <- c(1:4, 46:49, 79:82)
    X = list(miRNA = breast.TCGA$data.train$mirna[samples, 1:10],
            mRNA = breast.TCGA$data.train$mrna[samples, 1:10])
    Y = breast.TCGA$data.train$subtype[samples]
    
    keepX <- list(miRNA= c(5,5),
                  mRNA = c(5,5))
    
    res.block.splsda <- block.splsda(X, Y, keepX = keepX)
    
    set.seed(100)
    block.splsda.perf = perf(res.block.splsda, folds=2, nrepeat = 2)
    
    invisible(capture.output(TT <- dput(block.splsda.perf[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(perf.diablo:parameter): auc", {
    
    testable.components <- Testable.Components$auc.block.splsda
    GT <- Ground.Truths$auc.block.splsda
    
    data(breast.TCGA)
    samples <- c(1:4, 46:49, 79:82)
    X = list(miRNA = breast.TCGA$data.train$mirna[samples, 1:10],
            mRNA = breast.TCGA$data.train$mrna[samples, 1:10])
    Y = breast.TCGA$data.train$subtype[samples]
    
    keepX <- list(miRNA= c(5,5),
                  mRNA = c(5,5))
    
    res.block.splsda <- block.splsda(X, Y, keepX = keepX)
    
    set.seed(100)
    block.splsda.perf = perf(res.block.splsda, folds=3, auc = T)
    
    invisible(capture.output(TT <- dput(block.splsda.perf[testable.components])))
    
    expect_equal(TT, GT)
})


###############################################################################
### ================================ ERROR ================================ ###
###############################################################################


test_that("(perf.diablo:error): catches invalid values of validation", {
  
    data(breast.TCGA)
    samples <- c(1:4, 46:49, 79:82)
    X = list(miRNA = breast.TCGA$data.train$mirna[samples, 1:10],
            mRNA = breast.TCGA$data.train$mrna[samples, 1:10])
    Y = breast.TCGA$data.train$subtype[samples]
    
    keepX <- list(miRNA= c(5,5),
                  mRNA = c(5,5))
    
    res.block.splsda <- block.splsda(X, Y, keepX = keepX)
    
    expect_error(perf(res.block.splsda, validation = "random.method"),
                 "Choose 'validation' among the two following possibilities: 'Mfold' or 'loo'",
                 fixed=T)
})


test_that("(perf.diablo:error): catches invalid values of folds", {
  
    data(breast.TCGA)
    samples <- c(1:4, 46:49, 79:82)
    X = list(miRNA = breast.TCGA$data.train$mirna[samples, 1:10],
            mRNA = breast.TCGA$data.train$mrna[samples, 1:10])
    Y = breast.TCGA$data.train$subtype[samples]
    
    keepX <- list(miRNA= c(5,5),
                  mRNA = c(5,5))
    
    res.block.splsda <- block.splsda(X, Y, keepX = keepX)
    
    expect_error(perf(res.block.splsda, folds = "random.value"),
                 "Invalid number of folds.",
                 fixed=T)
    
    expect_error(perf(res.block.splsda, folds = 1),
                 "Invalid number of folds.",
                 fixed=T)
    
    expect_error(perf(res.block.splsda, folds = 13),
                 "Invalid number of folds.",
                 fixed=T)
    
    expect_warning(perf(res.block.splsda, folds = 5),
                   "At least one class is not represented in one fold, which may unbalance the error rate.\n  Consider a number of folds lower than the minimum in table(Y): 4",
                   fixed=T)
})