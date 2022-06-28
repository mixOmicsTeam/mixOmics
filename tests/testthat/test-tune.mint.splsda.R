
###############################################################################
### ============================ MISSING TESTS ============================ ###
###############################################################################


###############################################################################
### ============================ GROUND TRUTHS ============================ ###
###############################################################################   

Test.Data <- readRDS(system.file("testdata", "testdata-tune.mint.splsda.rda", package = "mixOmics"))
Testable.Components <- Test.Data$tc
Ground.Truths <- Test.Data$gt

###############################################################################
### ================================ BASIC ================================ ###
###############################################################################


test_that("(tune.mint.splsda:basic): basic", {
    
    testable.components <- Testable.Components$basic.tune.mint.splsda
    GT <- Ground.Truths$basic.tune.mint.splsda
    
    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    S <- stemcells$study

    d <- .minimal_train_test_subset(X, Y, S)
    
    test.keepX <- c(3,6,9)
    
    mint.splsda.tune <- suppressWarnings(tune.mint.splsda(d$X.tr, d$Y.tr,
                                                          study = d$S.tr,
                                                          test.keepX = test.keepX))
    
    invisible(capture.output(TT <- dput(mint.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


###############################################################################
### ================================ DATA ================================= ###
###############################################################################


test_that("(tune.mint.splsda:data): breast tcga train", {
    
    testable.components <- Testable.Components$breast.train.tune.mint.splsda
    GT <- Ground.Truths$breast.train.tune.mint.splsda
    
    data(breast.TCGA)
    X <- breast.TCGA$data.train$mrna
    Y <- breast.TCGA$data.train$subtype
    S <- rep(c(1,2,3), length(Y)/3)

    d <- .minimal_train_test_subset(X, Y, S)
    
    test.keepX <- c(3,6,9)
    
    breast.train.mint.splsda.tune <- suppressWarnings(tune.mint.splsda(d$X.tr, d$Y.tr,
                                                                       study = d$S.tr,
                                                                       test.keepX = test.keepX))
    
    invisible(capture.output(TT <- dput(breast.train.mint.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(tune.mint.splsda:data): breast tcga test", {
    
    testable.components <- Testable.Components$breast.test.tune.mint.splsda
    GT <- Ground.Truths$breast.test.tune.mint.splsda
    
    data(breast.TCGA)
    X <- breast.TCGA$data.test$mrna
    Y <- breast.TCGA$data.test$subtype
    S <- c(rep(c(1,2,3), length(Y)/3), 1)

    d <- .minimal_train_test_subset(X, Y, S)
    
    test.keepX <- c(3,6,9)
    
    breast.test.mint.splsda.tune <- suppressWarnings(tune.mint.splsda(d$X.tr, d$Y.tr,
                                                                       study = d$S.tr,
                                                                       test.keepX = test.keepX))
    
    invisible(capture.output(TT <- dput(breast.test.mint.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


###############################################################################
### ============================== PARAMETER ============================== ###
###############################################################################


test_that("(tune.mint.splsda:parameter): ncomp", {
    
    testable.components <- Testable.Components$ncomp.tune.mint.splsda
    GT <- Ground.Truths$ncomp.tune.mint.splsda
    
    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    S <- stemcells$study

    d <- .minimal_train_test_subset(X, Y, S)
    
    test.keepX <- c(3,6,9)
    
    ncomp.mint.splsda.tune <- suppressWarnings(tune.mint.splsda(d$X.tr, d$Y.tr,
                                                          study = d$S.tr,
                                                          test.keepX = test.keepX,
                                                          ncomp=4))
    
    invisible(capture.output(TT <- dput(ncomp.mint.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(tune.mint.splsda:parameter): already.tested.X", {
    
    testable.components <- Testable.Components$already.tested.X.tune.mint.splsda
    GT <- Ground.Truths$already.tested.X.tune.mint.splsda
    
    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    S <- stemcells$study

    d <- .minimal_train_test_subset(X, Y, S)
    
    test.keepX <- c(3,6,9)
    already.tested.X <- c(15)
    
    already.tested.X.mint.splsda.tune <- suppressWarnings(tune.mint.splsda(d$X.tr, d$Y.tr,
                                                          study = d$S.tr,
                                                          test.keepX = test.keepX,
                                                          ncomp=2, 
                                                          already.tested.X = already.tested.X))
    
    invisible(capture.output(TT <- dput(already.tested.X.mint.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(tune.mint.splsda:parameter): dist", {
    
    testable.components <- Testable.Components$dist.tune.mint.splsda
    
    
    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    S <- stemcells$study

    d <- .minimal_train_test_subset(X, Y, S)
    
    test.keepX <- c(3,6,9)
    
    GT <- Ground.Truths$max.dist.tune.mint.splsda
    
    max.dist.mint.splsda.tune <- suppressWarnings(tune.mint.splsda(d$X.tr, d$Y.tr,
                                                          study = d$S.tr,
                                                          test.keepX = test.keepX,
                                                          dist="max.dist"))
    
    invisible(capture.output(TT <- dput(max.dist.mint.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
    
    # --------- # 
    
    GT <- Ground.Truths$cent.mahal.dist.tune.mint.splsda
    
    centroids.dist.mint.splsda.tune <- suppressWarnings(tune.mint.splsda(d$X.tr, d$Y.tr,
                                                          study = d$S.tr,
                                                          test.keepX = test.keepX,
                                                          dist="centroids.dist"))
    
    invisible(capture.output(TT <- dput(centroids.dist.mint.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
    
    # --------- # 
    
    mahalanobis.dist.mint.splsda.tune <- suppressWarnings(tune.mint.splsda(d$X.tr, d$Y.tr,
                                                          study = d$S.tr,
                                                          test.keepX = test.keepX,
                                                          dist="mahalanobis.dist"))
    
    invisible(capture.output(TT <- dput(mahalanobis.dist.mint.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(tune.mint.splsda:parameter): measure", {
    
    testable.components <- Testable.Components$measure.tune.mint.splsda
    GT <- Ground.Truths$measure.tune.mint.splsda
    
    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    S <- stemcells$study

    d <- .minimal_train_test_subset(X, Y, S)
    
    test.keepX <- c(3,6,9)
    
    measure.mint.splsda.tune <- suppressWarnings(tune.mint.splsda(d$X.tr, d$Y.tr,
                                                                  study = d$S.tr,
                                                                  test.keepX = test.keepX,
                                                                  measure="overall"))
    
    invisible(capture.output(TT <- dput(measure.mint.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(tune.mint.splsda:parameter): auc", {
    
    testable.components <- Testable.Components$auc.tune.mint.splsda
    GT <- Ground.Truths$auc.tune.mint.splsda
    
    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    S <- stemcells$study

    d <- .minimal_train_test_subset(X, Y, S)
    
    test.keepX <- c(3,6,9)
    
    auc.mint.splsda.tune <- suppressWarnings(tune.mint.splsda(d$X.tr, d$Y.tr,
                                                          study = d$S.tr,
                                                          test.keepX = test.keepX,
                                                          auc=T))
    
    invisible(capture.output(TT <- dput(auc.mint.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(tune.mint.splsda:parameter): scale", {
    
    testable.components <- Testable.Components$scale.tune.mint.splsda
    GT <- Ground.Truths$scale.tune.mint.splsda
    
    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    S <- stemcells$study

    d <- .minimal_train_test_subset(X, Y, S)
    
    test.keepX <- c(3,6,9)
    
    scale.mint.splsda.tune <- suppressWarnings(tune.mint.splsda(d$X.tr, d$Y.tr,
                                                          study = d$S.tr,
                                                          test.keepX = test.keepX,
                                                          scale=F))
    
    invisible(capture.output(TT <- dput(scale.mint.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(tune.mint.splsda:parameter): tol", {
    
    testable.components <- Testable.Components$tol.tune.mint.splsda
    GT <- Ground.Truths$tol.tune.mint.splsda
    
    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    S <- stemcells$study

    d <- .minimal_train_test_subset(X, Y, S)
    
    test.keepX <- c(3,6,9)
    
    tol.mint.splsda.tune <- suppressWarnings(tune.mint.splsda(d$X.tr, d$Y.tr,
                                                          study = d$S.tr,
                                                          test.keepX = test.keepX,
                                                          tol=0.1))
    
    invisible(capture.output(TT <- dput(tol.mint.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(tune.mint.splsda:parameter): max.iter", {
    
    testable.components <- Testable.Components$max.iter.tune.mint.splsda
    GT <- Ground.Truths$max.iter.tune.mint.splsda
    
    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    S <- stemcells$study

    d <- .minimal_train_test_subset(X, Y, S)
    
    test.keepX <- c(3,6,9)
    
    max.iter.mint.splsda.tune <- suppressWarnings(tune.mint.splsda(d$X.tr, d$Y.tr,
                                                          study = d$S.tr,
                                                          test.keepX = test.keepX,
                                                          max.iter=2))
    
    invisible(capture.output(TT <- dput(max.iter.mint.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(tune.mint.splsda:parameter): near.zero.var", {

    testable.components <- Testable.Components$near.zero.var.tune.mint.splsda
    GT <- Ground.Truths$near.zero.var.tune.mint.splsda
    
    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    S <- stemcells$study
    
    X[1:115,1:390] <- 0

    d <- .minimal_train_test_subset(X, Y, S)
    
    test.keepX <- c(3,6,9)
    
    near.zero.var.mint.splsda.tune <- suppressWarnings(tune.mint.splsda(d$X.tr, d$Y.tr,
                                                          study = d$S.tr,
                                                          test.keepX = test.keepX,
                                                          near.zero.var=T))
    
    invisible(capture.output(TT <- dput(near.zero.var.mint.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(tune.mint.splsda:parameter): light.output", {

    testable.components <- Testable.Components$light.output.tune.mint.splsda
    GT <- Ground.Truths$light.output.tune.mint.splsda
    
    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    S <- stemcells$study

    d <- .minimal_train_test_subset(X, Y, S)
    
    test.keepX <- c(3,6,9)
    
    light.output.mint.splsda.tune <- suppressWarnings(tune.mint.splsda(d$X.tr, d$Y.tr,
                                                          study = d$S.tr,
                                                          test.keepX = test.keepX,
                                                          light.output=F))
    
    invisible(capture.output(TT <- dput(light.output.mint.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(tune.mint.splsda:parameter): signif.threshold", {
    
    testable.components <- Testable.Components$signif.threshold.tune.mint.splsda
    GT <- Ground.Truths$signif.threshold.tune.mint.splsda
    
    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    S <- stemcells$study

    d <- .minimal_train_test_subset(X, Y, S)
    
    test.keepX <- c(3,6,9)
    
    signif.threshold.mint.splsda.tune <- suppressWarnings(tune.mint.splsda(d$X.tr, d$Y.tr,
                                                          study = d$S.tr,
                                                          test.keepX = test.keepX,
                                                          signif.threshold=1e-12))
    
    invisible(capture.output(TT <- dput(signif.threshold.mint.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


###############################################################################
### ================================ ERROR ================================ ###
###############################################################################


test_that("(tune.mint.splsda:error): catches invalid 'X' objects", {
    
    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    S <- stemcells$study

    d <- .minimal_train_test_subset(X, Y, S)
    
    test.keepX <- c(3,6,9)
    
    expect_error(suppressWarnings(tune.mint.splsda(Y=d$Y.tr,
                                                          study = d$S.tr,
                                                          test.keepX = test.keepX)),
                "X")
    
    d$X.tr <- as.data.frame(apply(d$X.tr, 2, as.character))
    expect_error(suppressWarnings(tune.mint.splsda(X = d$X.tr, Y=d$Y.tr,
                                                          study = d$S.tr,
                                                          test.keepX = test.keepX)),
                "X")
})


test_that("(tune.mint.splsda:error): catches invalid 'Y' objects", {
    
    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    S <- stemcells$study

    d <- .minimal_train_test_subset(X, Y, S)
    
    test.keepX <- c(3,6,9)
    
    expect_error(suppressWarnings(tune.mint.splsda(X=d$X.tr,
                                                          study = d$S.tr,
                                                          test.keepX = test.keepX)),
                "Y")
    
    d$Y.tr <- NULL
    expect_error(suppressWarnings(tune.mint.splsda(X=d$X.tr, Y=d$Y.tr,
                                                          study = d$S.tr,
                                                          test.keepX = test.keepX)),
                "Y")
    
    d$Y.tr <- d$X.tr[,1:2]
    expect_error(suppressWarnings(tune.mint.splsda(X=d$X.tr, Y=d$Y.tr,
                                                          study = d$S.tr,
                                                          test.keepX = test.keepX)),
                "Y")
    
    d$Y.tr <- as.factor(rep("hESC", 12))
    expect_error(suppressWarnings(tune.mint.splsda(X=d$X.tr, Y=d$Y.tr,
                                                          study = d$S.tr,
                                                          test.keepX = test.keepX)),
                "Y")
})


test_that("(tune.mint.splsda:error): catches invalid 'progressBar' value", {
    
    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    S <- stemcells$study

    d <- .minimal_train_test_subset(X, Y, S)
    
    test.keepX <- c(3,6,9)
    
    expect_error(suppressWarnings(tune.mint.splsda(X=d$X.tr, Y = d$Y.tr,
                                                          study = d$S.tr,
                                                          test.keepX = test.keepX,
                                                          progressBar = "random.value")),
                "progressBar")
})


test_that("(tune.mint.splsda:error): catches invalid 'ncomp' value", {
    
    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    S <- stemcells$study

    d <- .minimal_train_test_subset(X, Y, S)
    
    test.keepX <- c(3,6,9)
    
    expect_error(suppressWarnings(tune.mint.splsda(X=d$X.tr, Y = d$Y.tr,
                                                          study = d$S.tr,
                                                          test.keepX = test.keepX,
                                                          ncomp = "random.value")),
                "ncomp")
    
    expect_error(suppressWarnings(tune.mint.splsda(X=d$X.tr, Y = d$Y.tr,
                                                          study = d$S.tr,
                                                          test.keepX = test.keepX,
                                                          ncomp = NULL)),
                "ncomp")
    
    expect_error(suppressWarnings(tune.mint.splsda(X=d$X.tr, Y = d$Y.tr,
                                                          study = d$S.tr,
                                                          test.keepX = test.keepX,
                                                          ncomp = -1)),
                "ncomp")
})


test_that("(tune.mint.splsda:error): catches invalid 'already.tested.X' value", {
    
    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    S <- stemcells$study

    d <- .minimal_train_test_subset(X, Y, S)
    
    test.keepX <- c(3,6,9)
    
    already.tested.X <- list(15)
    expect_error(suppressWarnings(tune.mint.splsda(d$X.tr, d$Y.tr,
                                                  study = d$S.tr,
                                                  test.keepX = test.keepX,
                                                  ncomp=2, 
                                                  already.tested.X = already.tested.X)),
                 "already.tested.X")
    
    already.tested.X <- c(15,15)
    expect_error(suppressWarnings(tune.mint.splsda(d$X.tr, d$Y.tr,
                                                  study = d$S.tr,
                                                  test.keepX = test.keepX,
                                                  ncomp=2, 
                                                  already.tested.X = already.tested.X)),
                 "already.tested.X")
})


test_that("(tune.mint.splsda:error): catches invalid 'study' objects", {
    
    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    S <- stemcells$study

    d <- .minimal_train_test_subset(X, Y, S)
    
    test.keepX <- c(3,6,9)
    
    expect_error(suppressWarnings(tune.mint.splsda(X = d$X.tr, Y=d$Y.tr,
                                                          test.keepX = test.keepX)),
                "study")
    
    d$S.tr <- d$S.tr[1:5]
    expect_error(suppressWarnings(tune.mint.splsda(X = d$X.tr, Y=d$Y.tr,
                                                   study = d$S.tr,
                                                   test.keepX = test.keepX)),
                 "study")
    
    d$S.tr <- as.factor(c(rep(1, 6), rep(2, 5), 3))
    expect_error(suppressWarnings(tune.mint.splsda(X = d$X.tr, Y=d$Y.tr,
                                                   study = d$S.tr,
                                                   test.keepX = test.keepX)),
                "study")
    
    d$S.tr <- as.factor(c(rep(1, 4), rep(2, 4), rep(3, 4)))
    expect_error(suppressWarnings(tune.mint.splsda(X = d$X.tr, Y=d$Y.tr,
                                                   study = d$S.tr,
                                                   test.keepX = test.keepX)),
                "study")
})


test_that("(tune.mint.splsda:error): catches invalid 'light.output' objects", {
    
    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    S <- stemcells$study

    d <- .minimal_train_test_subset(X, Y, S)
    
    test.keepX <- c(3,6,9)
    
    expect_error(suppressWarnings(tune.mint.splsda(X = d$X.tr, Y=d$Y.tr,
                                                   study=d$S.tr,
                                                          test.keepX = test.keepX,
                                                   light.output = "random.value")),
                "light.output")
})


test_that("(tune.mint.splsda:error): catches invalid 'test.keepX' objects", {
    
    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    S <- stemcells$study

    d <- .minimal_train_test_subset(X, Y, S)
    
    test.keepX <- NULL
    expect_error(suppressWarnings(tune.mint.splsda(X = d$X.tr, Y=d$Y.tr,
                                                   study=d$S.tr,
                                                   test.keepX = test.keepX)),
                "test.keepX")
    
    test.keepX <- c(3)
    expect_error(suppressWarnings(tune.mint.splsda(X = d$X.tr, Y=d$Y.tr,
                                                   study=d$S.tr,
                                                   test.keepX = test.keepX)),
                "test.keepX")
    
    test.keepX <- "c(3,6,9)"
    expect_error(suppressWarnings(tune.mint.splsda(X = d$X.tr, Y=d$Y.tr,
                                                   study=d$S.tr,
                                                   test.keepX = test.keepX)),
                "test.keepX")
})


###############################################################################
### =============================== WARNINGS ============================== ###
###############################################################################


test_that("(tune.mint.splsda:warning): warnings for 'study", {
    
    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    S <- stemcells$study

    d <- .minimal_train_test_subset(X, Y, S, n.tr=1)
    
    test.keepX <- c(3,6,9)
    
    expect_warning(tune.mint.splsda(X = d$X.tr, Y=d$Y.tr,
                                     study=d$S.tr,
                                     test.keepX = test.keepX),
                "5 samples")
    
    d$S.tr[12] <- 3
    expect_warning(tune.mint.splsda(X = d$X.tr, Y=d$Y.tr,
                                     study=d$S.tr,
                                     test.keepX = test.keepX),
                "all the levels")
})





