
###############################################################################
### ============================ MISSING TESTS ============================ ###
###############################################################################


###############################################################################
### ============================ GROUND TRUTHS ============================ ###
###############################################################################

Test.Data <- readRDS(system.file("testdata", "testdata-predict.rda", package = "mixOmics"))
Testable.Components <- Test.Data$tc
Ground.Truths <- Test.Data$gt

###############################################################################
### ================================ BASIC ================================ ###
###############################################################################


test_that("(predict:basic): pls", {
    
    testable.components <- Testable.Components$basic.pls
    GT <- Ground.Truths$basic.pls
    
    data(liver.toxicity)
    X <- liver.toxicity$gene
    Y <- liver.toxicity$clinic
    
    d <- .minimal_train_test_subset(X, Y, n.tr=6, n.te=3)
    
    res.pls <- pls(d$X.tr, d$Y.tr)
    
    pls.predict = predict(res.pls, d$X.te)
    
    invisible(capture.output(TT <- dput(pls.predict[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(predict:basic): spls", {
    
    testable.components <- Testable.Components$basic.spls
    GT <- Ground.Truths$basic.spls
    
    data(liver.toxicity)
    X <- liver.toxicity$gene
    Y <- liver.toxicity$clinic
    
    d <- .minimal_train_test_subset(X, Y, n.tr=6, n.te=3)
    
    choice.keepX = c(3,3)
    choice.keepY = c(3,3)
    
    res.spls <- spls(d$X.tr, d$Y.tr, keepX = choice.keepX, keepY = choice.keepY)
    
    spls.predict = predict(res.spls, d$X.te)
    
    invisible(capture.output(TT <- dput(spls.predict[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(predict:basic): plsda", {
    
    testable.components <- Testable.Components$basic.plsda
    GT <- Ground.Truths$basic.plsda
    
    data(srbct)
    
    X <- srbct$gene[, 1:10]
    Y <- srbct$class
    
    d <- .minimal_train_test_subset(X, Y)
    
    res.plsda <- plsda(d$X.tr, d$Y.tr)
    
    plsda.predict = predict(res.plsda, d$X.te)
    
    invisible(capture.output(TT <- dput(plsda.predict[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(predict:basic): splsda", {
    
    testable.components <- Testable.Components$basic.splsda
    GT <- Ground.Truths$basic.splsda
    
    data(srbct)
    
    d <- .minimal_train_test_subset(srbct$gene, srbct$class)
    
    choice.keepX = c(3,3)
    
    res.splsda <- splsda(d$X.tr, d$Y.tr, keepX = choice.keepX)
    
    splsda.predict = predict(res.splsda, d$X.te)
    
    invisible(capture.output(TT <- dput(splsda.predict[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(predict:basic): block.pls", {
    
    testable.components <- Testable.Components$basic.block.pls
    GT <- Ground.Truths$basic.block.pls
    
    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.train$mirna[,1:5],
             mRNA = breast.TCGA$data.train$mrna[,1:5])
    Y <- breast.TCGA$data.train$protein[,1:5]
    
    d <- .minimal_train_test_subset(X, Y, n.tr=6, n.te=3)
    
    res.block.pls <- block.pls(d$X.tr, d$Y.tr)
    
    block.pls.predict = predict(res.block.pls, d$X.te)
    
    invisible(capture.output(TT <- dput(block.pls.predict[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(predict:basic): block.spls", {
    
    testable.components <- Testable.Components$basic.block.spls
    GT <- Ground.Truths$basic.block.spls
    
    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.train$mirna[,1:5],
             mRNA = breast.TCGA$data.train$mrna[,1:5])
    Y <- breast.TCGA$data.train$protein[,1:5]
    
    d <- .minimal_train_test_subset(X, Y, n.tr=6, n.te=3)
    
    choice.keepX <- list(miRNA = c(3,3),
                         mRNA = c(3,3))
    
    res.block.spls <- block.spls(d$X.tr, d$Y.tr, keepX = choice.keepX)
    
    block.spls.predict = predict(res.block.spls, d$X.te)
    
    invisible(capture.output(TT <- dput(block.spls.predict[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(predict:basic): block.plsda", {

    testable.components <- Testable.Components$basic.block.plsda
    GT <- Ground.Truths$basic.block.plsda

    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.train$mirna[,1:10],
             mRNA = breast.TCGA$data.train$mrna[,1:10])
    Y = breast.TCGA$data.train$subtype

    d <- .minimal_train_test_subset(X, Y)

    res.block.plsda <- block.plsda(d$X.tr, d$Y.tr)

    block.plsda.predict = predict(res.block.plsda, d$X.te)

    invisible(capture.output(TT <- dput(block.plsda.predict[testable.components])))

    expect_equal(TT, GT)
})


test_that("(predict:basic): block.splsda", {

    testable.components <- Testable.Components$basic.block.splsda
    GT <- Ground.Truths$basic.block.splsda

    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.train$mirna[,1:10],
             mRNA = breast.TCGA$data.train$mrna[,1:10])
    Y = breast.TCGA$data.train$subtype

    d <- .minimal_train_test_subset(X, Y)

    choice.keepX <- list(miRNA = c(3,3),
                         mRNA = c(3,3))

    res.block.splsda <- block.splsda(d$X.tr, d$Y.tr, keepX = choice.keepX)

    block.splsda.predict = predict(res.block.splsda, d$X.te)

    invisible(capture.output(TT <- dput(block.splsda.predict[testable.components])))

    expect_equal(TT, GT)
})


test_that("(predict:basic): mint.pls", {

    testable.components <- Testable.Components$basic.mint.pls
    GT <- Ground.Truths$basic.mint.pls

    data(stemcells)
    X <- stemcells$gene[, 1:10]
    Y <- stemcells$gene[, 11:20]
    S <- as.character(stemcells$study)

    d <- .minimal_train_test_subset(X, Y, S=S)

    res.mint.pls <- suppressWarnings(mint.pls(d$X.tr, d$Y.tr, study = d$S.tr))

    mint.pls.predict = predict(res.mint.pls, d$X.te, d$S.te)

    invisible(capture.output(TT <- dput(mint.pls.predict[testable.components])))

    expect_equal(TT, GT)
})


# test_that("(predict:basic): mint.spls", {
# 
#     testable.components <- Testable.Components$basic.mint.spls
#     GT <- Ground.Truths$basic.mint.spls
# 
#     data(stemcells)
#     X <- stemcells$gene[, 1:10]
#     Y <- stemcells$gene[, 11:20]
#     S <- as.character(stemcells$study)
# 
#     d <- .minimal_train_test_subset(X, Y, S=S)
# 
#     choice.keepX <- c(3,3)
#     choice.keepY <- c(3,3)
#     
#     set.seed(42)
#     res.mint.spls <- suppressWarnings(mint.spls(d$X.tr, d$Y.tr, study = d$S.tr,
#                                                 keepX = choice.keepX,
#                                                 keepY = choice.keepY))
# 
#     mint.spls.predict = predict(res.mint.spls, d$X.te, d$S.te)
# 
#     invisible(capture.output(TT <- dput(mint.spls.predict[testable.components])))
# 
#     expect_equal(TT, GT)
# })


test_that("(predict:basic): mint.plsda", {

    testable.components <- Testable.Components$basic.mint.plsda
    GT <- Ground.Truths$basic.mint.plsda

    data(stemcells)
    X <- stemcells$gene[,1:10]
    Y <- stemcells$celltype
    S <- stemcells$study

    d <- .minimal_train_test_subset(X, Y, S=S)

    res.mint.plsda <- suppressWarnings(mint.plsda(d$X.tr, d$Y.tr, study=d$S.tr))

    mint.plsda.predict = predict(res.mint.plsda, d$X.te, study=d$S.te)

    invisible(capture.output(TT <- dput(mint.plsda.predict[testable.components])))

    expect_equal(TT, GT)
})


test_that("(predict:basic): mint.splsda", {

    testable.components <- Testable.Components$basic.mint.splsda
    GT <- Ground.Truths$basic.mint.splsda

    data(stemcells)
    X <- stemcells$gene[,1:10]
    Y <- stemcells$celltype
    S <- stemcells$study

    d <- .minimal_train_test_subset(X, Y, S=S)

    choice.keepX <- c(3,3)

    res.mint.splsda <- suppressWarnings(mint.splsda(d$X.tr, d$Y.tr, study=d$S.tr, keepX = choice.keepX))

    mint.splsda.predict = predict(res.mint.splsda, d$X.te, d$S.te)

    invisible(capture.output(TT <- dput(mint.splsda.predict[testable.components])))

    expect_equal(TT, GT)
})


###############################################################################
### ================================ DATA ================================= ###
###############################################################################


test_that("(predict:data): breast.test", {

    testable.components <- Testable.Components$breast.test.block.splsda
    GT <- Ground.Truths$breast.test.block.splsda

    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.test$mirna[,1:10],
             mRNA = breast.TCGA$data.test$mrna[,1:10])
    Y = breast.TCGA$data.test$subtype

    d <- .minimal_train_test_subset(X, Y)

    choice.keepX <- list(miRNA = c(3,3),
                         mRNA = c(3,3))

    res.block.splsda <- block.splsda(d$X.tr, d$Y.tr, keepX = choice.keepX)

    breast.test.predict = predict(res.block.splsda, d$X.te)

    invisible(capture.output(TT <- dput(breast.test.predict[testable.components])))

    expect_equal(TT, GT)
})


###############################################################################
### ============================== PARAMETER ============================== ###
###############################################################################


test_that("(predict:parameter): dist", {

    testable.components <- Testable.Components$dist.plsda
    GT <- Ground.Truths$dist.plsda

    data(srbct)

    d <- .minimal_train_test_subset(srbct$gene, srbct$class)

    res.plsda <- plsda(d$X.tr, d$Y.tr)

    # ------------------------------------------------------------------------#
    # centroids.dist

    centroids.predict = predict(res.plsda, d$X.te, dist = "centroids.dist")

    invisible(capture.output(TT <- dput(centroids.predict[testable.components])))
    expect_equal(TT, GT)

    # ------------------------------------------------------------------------#
    # mahalanobis.dist

    mahalanobis.predict = predict(res.plsda, d$X.te, dist = "mahalanobis.dist")

    invisible(capture.output(TT <- dput(mahalanobis.predict[testable.components])))
    expect_equal(TT, GT)
})


# test_that("(predict:parameter): multilevel", {
# 
#     testable.components <- Testable.Components$multilevel.plsda
#     GT <- Ground.Truths$multilevel.plsda
# 
#     data(vac18)
#     X <- vac18$genes[, 1:10]
#     Y <- vac18$stimulation
#     ML <- vac18$sample
# 
#     d <- .minimal_train_test_subset(X, Y, ML=ML)
# 
#     res.plsda <- plsda(d$X.tr, d$Y.tr, multilevel = d$ML.tr)
# 
#     multilevel.predict = predict(res.plsda, d$X.te, multilevel = d$ML.te)
# 
#     invisible(capture.output(TT <- dput(multilevel.predict[testable.components])))
# 
#     expect_equal(TT, GT)
# 
# })


###############################################################################
### ================================ ERROR ================================ ###
###############################################################################


test_that("(predict:error): catches invalid values for 'dist'", {
    
    data(srbct)
    
    X <- srbct$gene[, 1:10]
    Y <- srbct$class
    
    d <- .minimal_train_test_subset(X, Y)
    
    res.plsda <- plsda(d$X.tr, d$Y.tr)
    
    expect_error(predict(res.plsda, d$X.te, dist = "random.dist"),
                 "ERROR : choose one of the four following modes: 'all', 'max.dist', 'centroids.dist' or 'mahalanobis.dist'",
                 fixed=T)
})


test_that("(predict:error): catches invalid usage on 'rcgga' objects", {
    
    data(nutrimouse)
    data = list(gene = nutrimouse$gene, 
                lipid = nutrimouse$lipid, 
                Y = unmap(nutrimouse$diet))
    
    design = matrix(c(0,1,1,
                      1,0,1,
                      1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

    
    res.rgcca <-  wrapper.rgcca(X = data, design = design, tau = c(1, 1, 0),
                                ncomp = 2,
                                scheme = "centroid")

    
    expect_error(predict(res.rgcca, data),
                 "no applicable method for 'predict' applied to an object of class \"c('sparse.rgcca', 'rgcca')\"",
                 fixed=T)
})


test_that("(predict:error): catches invalid values for 'newdata'", {
    
    data(srbct)
    
    X <- srbct$gene[, 1:10]
    Y <- srbct$class
    
    d <- .minimal_train_test_subset(X, Y)
    
    res.plsda <- plsda(d$X.tr, d$Y.tr)
    
    X.te <- data.frame(matrix("val", ncol=ncol(d$X.te), nrow=nrow(d$X.te)))
    colnames(X.te) <- colnames(d$X.te)
    
    expect_error(predict(res.plsda, X.te),
                 "'X[[1]]'  must be a numeric matrix.",
                 fixed=T)
    
    d$X.te <- d$X.te[, 1:9]
    expect_error(suppressWarnings(predict(res.plsda, d$X.te)),
                 "'newdata' must include all the variables of 'object$X'",
                 fixed=T)
})


test_that("(predict:error): catches invalid 'newdata' for block objects", {
    
    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.train$mirna[,1:10],
             mRNA = breast.TCGA$data.train$mrna[,1:10])
    Y = breast.TCGA$data.train$subtype

    d <- .minimal_train_test_subset(X, Y)

    res.block.plsda <- block.plsda(d$X.tr, d$Y.tr)
    
    expect_error(predict(res.block.plsda, d$X.te[[1]]),
                 "'newdata' should be a list",
                 fixed=T)
})


###############################################################################
### ============================== EDGE CASES ============================= ###
###############################################################################