
###############################################################################
### ============================ MISSING TESTS ============================ ###
###############################################################################



###############################################################################
### ============================ GROUND TRUTHS ============================ ###
###############################################################################

Test.Data <- readRDS(system.file("testdata", "testdata-tune.splsda.rda", package = "mixOmics"))
Testable.Components <- Test.Data$tc
Ground.Truths <- Test.Data$gt

###############################################################################
### ================================ BASIC ================================ ###
###############################################################################


test_that("(tune.splsda:basic): basic", {
    
    testable.components <- Testable.Components$basic.tune.splsda
    GT <- Ground.Truths$basic.tune.splsda
    
    data(srbct)
    X <- srbct$gene
    Y <- srbct$class
    
    d <- .minimal_train_test_subset(X, Y)
    
    test.keepX <- c(3,6,9)
    
    splsda.tune <- tune.splsda(X=d$X.tr, Y=d$Y.tr,
                               test.keepX = test.keepX,
                               folds = 2)
    
    invisible(capture.output(TT <- dput(splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


###############################################################################
### ================================ DATA ================================= ###
###############################################################################


test_that("(tune.splsda:data): breast tcga train", {
    
    testable.components <- Testable.Components$breast.train.tune.splsda
    GT <- Ground.Truths$breast.train.tune.splsda
    
    data(breast.TCGA)
    X <- breast.TCGA$data.train$mrna
    Y <- breast.TCGA$data.train$subtype
    
    d <- .minimal_train_test_subset(X, Y)
    
    test.keepX <- c(3,6,9)
    
    breast.train.splsda.tune <- tune.splsda(X=d$X.tr, Y=d$Y.tr,
                                            test.keepX = test.keepX,
                                            folds = 2)
    
    invisible(capture.output(TT <- dput(breast.train.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(tune.splsda:data): breast tcga test", {
    
    testable.components <- Testable.Components$breast.test.tune.splsda
    GT <- Ground.Truths$breast.test.tune.splsda
    
    data(breast.TCGA)
    X <- breast.TCGA$data.test$mrna
    Y <- breast.TCGA$data.test$subtype
    
    d <- .minimal_train_test_subset(X, Y)
    
    test.keepX <- c(3,6,9)
    
    breast.test.splsda.tune <- tune.splsda(X=d$X.tr, Y=d$Y.tr,
                                            test.keepX = test.keepX,
                                            folds = 2)
    
    invisible(capture.output(TT <- dput(breast.test.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(tune.splsda:data): stem cells", {
    
    testable.components <- Testable.Components$stem.cells.tune.splsda
    GT <- Ground.Truths$stem.cells.tune.splsda
    
    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    
    d <- .minimal_train_test_subset(X, Y)
    
    test.keepX <- c(3,6,9)
    
    stem.cells.splsda.tune <- tune.splsda(X=d$X.tr, Y=d$Y.tr,
                                           test.keepX = test.keepX,
                                           folds = 2)
    
    invisible(capture.output(TT <- dput(stem.cells.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


###############################################################################
### ============================== PARAMETER ============================== ###
###############################################################################


test_that("(tune.splsda:parameter): ncomp", {
    
    testable.components <- Testable.Components$ncomp.tune.splsda
    GT <- Ground.Truths$ncomp.tune.splsda
    
    data(srbct)
    X <- srbct$gene
    Y <- srbct$class
    
    d <- .minimal_train_test_subset(X, Y)
    
    test.keepX <- c(3,6,9)
    
    ncomp.splsda.tune <- tune.splsda(X=d$X.tr, Y=d$Y.tr,
                               test.keepX = test.keepX,
                               folds = 2,
                               ncomp=3)
    
    invisible(capture.output(TT <- dput(ncomp.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(tune.splsda:parameter): test.keepX", {
    
    testable.components <- Testable.Components$test.keepX.tune.splsda
    GT <- Ground.Truths$test.keepX.tune.splsda
    
    data(srbct)
    X <- srbct$gene
    Y <- srbct$class
    
    d <- .minimal_train_test_subset(X, Y)
    
    test.keepX <- seq(1, 20, 1)
    
    test.keepX.splsda.tune <- tune.splsda(X=d$X.tr, Y=d$Y.tr,
                               test.keepX = test.keepX,
                               folds = 2)
    
    invisible(capture.output(TT <- dput(test.keepX.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(tune.splsda:parameter): already.tested.X", {
    
    testable.components <- Testable.Components$already.tested.X.tune.splsda
    GT <- Ground.Truths$already.tested.X.tune.splsda
    
    data(srbct)
    X <- srbct$gene
    Y <- srbct$class
    
    d <- .minimal_train_test_subset(X, Y)
    
    test.keepX <- c(3,6,9)
    already.tested.X <- c(15)
    
    already.tested.X.splsda.tune <- tune.splsda(X=d$X.tr, Y=d$Y.tr,
                               test.keepX = test.keepX,
                               folds = 2,
                               already.tested.X = already.tested.X,
                               ncomp=2)
    
    invisible(capture.output(TT <- dput(already.tested.X.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(tune.splsda:parameter): validation", {
    
    testable.components <- Testable.Components$validation.tune.splsda
    GT <- Ground.Truths$validation.tune.splsda
    
    data(srbct)
    X <- srbct$gene
    Y <- srbct$class
    
    d <- .minimal_train_test_subset(X, Y)
    
    test.keepX <- c(3,6,9)
    
    validation.splsda.tune <- tune.splsda(X=d$X.tr, Y=d$Y.tr,
                               test.keepX = test.keepX,
                               validation = "loo")
    
    invisible(capture.output(TT <- dput(validation.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(tune.splsda:parameter): dist", {
    
    testable.components <- Testable.Components$dist.tune.splsda
    
    data(srbct)
    X <- srbct$gene
    Y <- srbct$class
    
    d <- .minimal_train_test_subset(X, Y, n.tr=3)
    
    test.keepX <- c(3,6,9)
    
    set.seed(16)
    centroids.dist.splsda.tune <- tune.splsda(X=d$X.tr, Y=d$Y.tr,
                                    test.keepX = test.keepX,
                                    folds = 2,
                                    dist = "centroids.dist")
    
    invisible(capture.output(TT <- dput(centroids.dist.splsda.tune[testable.components])))
    
    GT <- Ground.Truths$centroids.dist.tune.splsda
    expect_equal(TT, GT)
    
    set.seed(16)
    mahalanobis.dist.splsda.tune <- tune.splsda(X=d$X.tr, Y=d$Y.tr,
                                    test.keepX = test.keepX,
                                    folds = 2,
                                    dist = "mahalanobis.dist")
    
    invisible(capture.output(TT <- dput(mahalanobis.dist.splsda.tune[testable.components])))
    
    GT <- Ground.Truths$mahalanobis.dist.tune.splsda
    expect_equal(TT, GT)
})


test_that("(tune.splsda:parameter): measure", {
    
    testable.components <- Testable.Components$measure.tune.splsda
    
    data(srbct)
    X <- srbct$gene
    Y <- srbct$class
    
    d <- .minimal_train_test_subset(X, Y)
    
    test.keepX <- c(3,6,9)
    
    set.seed(16)
    overall.measure.splsda.tune <- tune.splsda(X=d$X.tr, Y=d$Y.tr,
                               test.keepX = test.keepX,
                               folds = 2,
                               measure="overall")
    
    invisible(capture.output(TT <- dput(overall.measure.splsda.tune[testable.components])))
    
    GT <- Ground.Truths$overall.measure.tune.splsda
    expect_equal(TT, GT)
    
    set.seed(16)
    auc.measure.splsda.tune <- tune.splsda(X=d$X.tr, Y=d$Y.tr,
                               test.keepX = test.keepX,
                               folds = 2,
                               measure="AUC")
    
    invisible(capture.output(TT <- dput(auc.measure.splsda.tune[testable.components])))
    
    GT <- Ground.Truths$auc.measure.tune.splsda
    expect_equal(TT, GT)
})


test_that("(tune.splsda:parameter): scale", {
    
    testable.components <- Testable.Components$scale.tune.splsda
    GT <- Ground.Truths$scale.tune.splsda
    
    data(srbct)
    X <- srbct$gene
    Y <- srbct$class
    
    d <- .minimal_train_test_subset(X, Y)
    
    test.keepX <- c(3,6,9)
    
    scale.splsda.tune <- tune.splsda(X=d$X.tr, Y=d$Y.tr,
                               test.keepX = test.keepX,
                               folds = 2,
                               scale=F)
    
    invisible(capture.output(TT <- dput(scale.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(tune.splsda:parameter): auc", {
    
    testable.components <- Testable.Components$auc.tune.splsda
    GT <- Ground.Truths$auc.tune.splsda
    
    data(srbct)
    X <- srbct$gene
    Y <- srbct$class
    
    d <- .minimal_train_test_subset(X, Y)
    
    test.keepX <- c(3,6,9)
    
    set.seed(16)
    auc.splsda.tune <- tune.splsda(X=d$X.tr, Y=d$Y.tr,
                               test.keepX = test.keepX,
                               folds = 2,
                               auc=T, nrepeat = 3)
    
    invisible(capture.output(TT <- dput(auc.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(tune.splsda:parameter): progressBar", {
    
    testable.components <- Testable.Components$progressBar.tune.splsda
    GT <- Ground.Truths$progressBar.tune.splsda
    
    data(srbct)
    X <- srbct$gene
    Y <- srbct$class
    
    d <- .minimal_train_test_subset(X, Y)
    
    test.keepX <- c(3,6,9)
    
    progressBar.splsda.tune <- quiet(tune.splsda(X=d$X.tr, Y=d$Y.tr,
                                                    test.keepX = test.keepX,
                                                    folds = 2,
                                                    progressBar=T))
    
    invisible(capture.output(TT <- dput(progressBar.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(tune.splsda:parameter): tol", {
    
    testable.components <- Testable.Components$tol.tune.splsda
    GT <- Ground.Truths$tol.tune.splsda
    
    data(srbct)
    X <- srbct$gene
    Y <- srbct$class
    
    d <- .minimal_train_test_subset(X, Y)
    
    test.keepX <- c(3,6,9)
    
    tol.splsda.tune <- tune.splsda(X=d$X.tr, Y=d$Y.tr,
                               test.keepX = test.keepX,
                               folds = 2,
                               tol=1e-12)
    
    invisible(capture.output(TT <- dput(tol.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(tune.splsda:parameter): max.iter", {
    
    testable.components <- Testable.Components$max.iter.tune.splsda
    GT <- Ground.Truths$max.iter.tune.splsda
    
    data(srbct)
    X <- srbct$gene
    Y <- srbct$class
    
    d <- .minimal_train_test_subset(X, Y)
    
    test.keepX <- c(3,6,9)
    
    max.iter.splsda.tune <- tune.splsda(X=d$X.tr, Y=d$Y.tr,
                               test.keepX = test.keepX,
                               folds = 2,
                               max.iter=2)
    
    invisible(capture.output(TT <- dput(max.iter.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(tune.splsda:parameter): near.zero.var", {
    
    testable.components <- Testable.Components$near.zero.var.tune.splsda
    GT <- Ground.Truths$near.zero.var.tune.splsda
    
    data(srbct)
    X <- srbct$gene
    Y <- srbct$class
    
    X[1:60, 1:2000] <- 0
    
    d <- .minimal_train_test_subset(X, Y)
    
    test.keepX <- c(3,6,9)
    
    near.zero.var.splsda.tune <- suppressWarnings(tune.splsda(X=d$X.tr, Y=d$Y.tr,
                                                               test.keepX = test.keepX,
                                                               folds = 2,
                                                               near.zero.var=T))
    
    invisible(capture.output(TT <- dput(near.zero.var.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(tune.splsda:parameter): nrepeat", {
    
    testable.components <- Testable.Components$nrepeat.tune.splsda
    GT <- Ground.Truths$nrepeat.tune.splsda
    
    data(srbct)
    X <- srbct$gene
    Y <- srbct$class
    
    d <- .minimal_train_test_subset(X, Y)
    
    test.keepX <- c(3,6,9)
    
    set.seed(16)
    nrepeat.splsda.tune <- tune.splsda(X=d$X.tr, Y=d$Y.tr,
                               test.keepX = test.keepX,
                               folds = 2,
                               nrepeat=10)
    
    invisible(capture.output(TT <- dput(nrepeat.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(tune.splsda:parameter): logratio", {
    
    testable.components <- Testable.Components$logratio.tune.splsda
    GT <- Ground.Truths$logratio.tune.splsda
    
    data(srbct)
    X <- srbct$gene
    Y <- srbct$class
    
    d <- .minimal_train_test_subset(X, Y)
    
    test.keepX <- c(3,6,9)
    
    set.seed(16)
    logratio.splsda.tune <- tune.splsda(X=d$X.tr, Y=d$Y.tr,
                               test.keepX = test.keepX,
                               folds = 2,
                               logratio="CLR")
    
    invisible(capture.output(TT <- dput(logratio.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(tune.splsda:parameter): multilevel", {
    
    testable.components <- Testable.Components$multilevel.tune.splsda
    GT <- Ground.Truths$multilevel.tune.splsda
    
    data(vac18)
    X <- vac18$genes
    Y <- vac18$stimulation
    ML <- vac18$sample
    
    d <- .minimal_train_test_subset(X, Y, ML=ML, n.tr=4)
    
    test.keepX <- c(3,6,9)
    
    set.seed(16)
    multilevel.splsda.tune <- tune.splsda(X=d$X.tr, Y=d$Y.tr,
                               test.keepX = test.keepX,
                               folds = 2,
                               multilevel=d$ML.tr)
    
    invisible(capture.output(TT <- dput(multilevel.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(tune.splsda:parameter): light.output", {
    
    testable.components <- Testable.Components$light.output.tune.splsda
    GT <- Ground.Truths$light.output.tune.splsda
    
    data(srbct)
    X <- srbct$gene
    Y <- srbct$class
    
    d <- .minimal_train_test_subset(X, Y)
    
    test.keepX <- c(3,6,9)
    
    light.output.splsda.tune <- tune.splsda(X=d$X.tr, Y=d$Y.tr,
                               test.keepX = test.keepX,
                               folds = 2,
                               light.output=F)
    
    invisible(capture.output(TT <- dput(light.output.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(tune.splsda:parameter): signif.threshold", {
    
    testable.components <- Testable.Components$signif.threshold.tune.splsda
    GT <- Ground.Truths$signif.threshold.tune.splsda
    
    data(srbct)
    X <- srbct$gene
    Y <- srbct$class
    
    d <- .minimal_train_test_subset(X, Y)
    
    test.keepX <- c(3,6,9)
    
    signif.threshold.splsda.tune <- tune.splsda(X=d$X.tr, Y=d$Y.tr,
                               test.keepX = test.keepX,
                               folds = 2,
                               signif.threshold=1e-12)
    
    invisible(capture.output(TT <- dput(signif.threshold.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


###############################################################################
### ================================ ERROR ================================ ###
###############################################################################


test_that("(tune.splsda:error): catches invalid 'X' objects", {
    
    data(srbct)
    X <- srbct$gene
    Y <- srbct$class
    
    test.keepX <- c(3,6,9)
    
    expect_error(tune.splsda(Y=Y,
                           test.keepX = test.keepX,
                           folds = 2),
                "X")
    
    X <- as.data.frame(apply(X, 2, as.character))
    expect_error(tune.splsda(X=X, Y=Y,
                               test.keepX = test.keepX,
                               folds = 2),
                "X")
})


test_that("(tune.splsda:error): catches invalid 'Y' objects", {
    
    data(srbct)
    X <- srbct$gene
    Y <- srbct$class
    
    test.keepX <- c(3,6,9)
    
    expect_error(tune.splsda(X=X,
                           test.keepX = test.keepX,
                           folds = 2),
                "Y")
    
    expect_error(tune.splsda(X=X, Y=NULL,
                           test.keepX = test.keepX,
                           folds = 2),
                "Y")
    
    expect_error(tune.splsda(X=X, Y=X,
                           test.keepX = test.keepX,
                           folds = 2),
                "Y")
    
    
    expect_error(tune.splsda(X=X, Y=as.factor(rep("EWS", nrow(X))),
                           test.keepX = test.keepX,
                           folds = 2),
                "Y")
})


test_that("(tune.splsda:error): catches invalid 'multilevel' objects", {
    
    data(vac18)
    X <- vac18$genes
    Y <- vac18$stimulation
    ML <- vac18$sample
    
    test.keepX <- c(3,6,9)
     
    expect_error(tune.splsda(X=X, Y=Y,
                              test.keepX = test.keepX,
                              folds = 2,
                              multilevel=ML[1:10]),
                "multilevel")
    
    expect_error(tune.splsda(X=X, Y=Y,
                              test.keepX = test.keepX,
                              folds = 2,
                              multilevel=cbind(ML, ML)),
                "multilevel")
    
    expect_error(tune.splsda(X=X, Y=cbind(Y,Y,Y),
                           test.keepX = test.keepX,
                           folds = 2,
                           multilevel=ML),
                "Y")
})


test_that("(tune.splsda:error): catches invalid 'progressBar' value", {
    
    data(srbct)
    X <- srbct$gene
    Y <- srbct$class
    
    test.keepX <- c(3,6,9)
    
    expect_error(tune.splsda(X=X, Y=Y,
                           test.keepX = test.keepX,
                           folds = 2,
                           progressBar = "random.value"),
                "progressBar")
})


test_that("(tune.splsda:error): catches invalid 'ncomp' values", {
    
    data(srbct)
    X <- srbct$gene
    Y <- srbct$class
    
    test.keepX <- c(3,6,9)
    
    expect_error(tune.splsda(X=X, Y=Y,
                           test.keepX = test.keepX,
                           folds = 2,
                           ncomp = NULL),
                "ncomp")
    
    expect_error(tune.splsda(X=X, Y=Y,
                           test.keepX = test.keepX,
                           folds = 2,
                           ncomp = "random.value"),
                "ncomp")
    
    expect_error(tune.splsda(X=X, Y=Y,
                           test.keepX = test.keepX,
                           folds = 2,
                           ncomp = -1),
                "ncomp")
})


test_that("(tune.splsda:error): catches invalid 'validation' values", {
    
    data(srbct)
    X <- srbct$gene
    Y <- srbct$class
    
    test.keepX <- c(3,6,9)
    
    expect_error(tune.splsda(X=X, Y=Y,
                           test.keepX = test.keepX,
                           folds = 2,
                           validation = "random.value"),
                "validation")
})


test_that("(tune.mint.splsda:error): catches invalid 'already.tested.X' value", {
    
    data(srbct)
    X <- srbct$gene
    Y <- srbct$class
    
    test.keepX <- c(3,6,9)
    
    already.tested.X <- list()
    expect_error(tune.splsda(X=X, Y=Y,
                           test.keepX = test.keepX,
                           folds = 2,
                           already.tested.X = already.tested.X),
                 "already.tested.X")
    
    already.tested.X <- list(15)
    expect_error(tune.splsda(X=X, Y=Y,
                           test.keepX = test.keepX,
                           folds = 2,
                           already.tested.X = already.tested.X),
                 "already.tested.X")
    
    already.tested.X <- c(15,15)
    expect_error(tune.splsda(X=X, Y=Y,
                           test.keepX = test.keepX,
                           folds = 2,
                           ncomp=2,
                           already.tested.X = already.tested.X),
                 "already.tested.X")
})


test_that("(tune.mint.splsda:error): catches invalid 'measure' value", {
    
    data(srbct)
    X <- srbct$gene
    Y <- srbct$class
    
    test.keepX <- c(3,6,9)
    
    expect_error(tune.splsda(X=X, Y=Y,
                           test.keepX = test.keepX,
                           folds = 2,
                           measure = "random.value"),
                 "measure")
})


test_that("(tune.mint.splsda:error): catches invalid 'folds' value", {
    
    data(srbct)
    X <- srbct$gene
    Y <- srbct$class
    
    test.keepX <- c(3,6,9)
    
    expect_error(tune.splsda(X=X, Y=Y,
                           test.keepX = test.keepX,
                           folds = NULL),
                 "folds")
    
    expect_error(tune.splsda(X=X, Y=Y,
                           test.keepX = test.keepX,
                           folds = 1),
                 "folds")
    
    expect_error(tune.splsda(X=X, Y=Y,
                           test.keepX = test.keepX,
                           folds = 64),
                 "folds")
    
    expect_error(tune.splsda(X=X, Y=Y,
                           test.keepX = test.keepX,
                           folds = "random.value"),
                 "folds")
})


test_that("(tune.mint.splsda:error): catches invalid 'validation' value", {
    
    data(srbct)
    X <- srbct$gene
    Y <- srbct$class
    
    test.keepX <- c(3,6,9)
    
    expect_error(tune.splsda(X=X, Y=Y,
                           test.keepX = test.keepX,
                           folds = 2,
                           validation = "random.value"),
                 "validation")
})


test_that("(tune.mint.splsda:error): catches invalid 'test.keepX' value", {
    
    data(srbct)
    X <- srbct$gene
    Y <- srbct$class
    
    expect_error(tune.splsda(X=X, Y=Y,
                           test.keepX = NULL,
                           folds = 2),
                 "test.keepX")
    
    test.keepX <- c(3,2400, 2500, 2600)
    expect_error(tune.splsda(X=X, Y=Y,
                           test.keepX = test.keepX,
                           folds = 2),
                 "test.keepX")
})


test_that("(tune.mint.splsda:error): catches `X` dataframe with non-unique dimnames", {
    
    data(srbct)
    X <- srbct$gene
    Y <- srbct$class
    
    test.keepX <- c(3,6,9)
    
    rownames(X)[1] <- rownames(X)[2]
    expect_error(tune.splsda(X=X, Y=Y,
                           test.keepX = test.keepX,
                           folds = 2),
                 "unique")
    
    
    X <- srbct$gene
    colnames(X)[1] <- colnames(X)[2]
    expect_error(tune.splsda(X=X, Y=Y,
                           test.keepX = test.keepX,
                           folds = 2),
                 "Unique")
})


test_that("(tune.mint.splsda:error): catches when all features of 'X' have near zero var", {
    
    data(srbct)
    X <- srbct$gene
    Y <- srbct$class
    
    X[1:60,] <- 0
    
    test.keepX <- c(3,6,9)
    
    expect_error(suppressWarnings(tune.splsda(X=X, Y=Y,
                               test.keepX = test.keepX,
                               folds = 2,
                               near.zero.var=T)),
                 "Near Zero Var")
})

 


###############################################################################
### =============================== WARNINGS ============================== ###
###############################################################################


test_that("(tune.splsda:warning): catches when 'validation = loo' and 'nrepeat>1", {
    
    data(srbct)
    X <- srbct$gene
    Y <- srbct$class
    
    test.keepX <- c(3,6,9)
    
    expect_warning(tune.splsda(X=X, Y=Y,
                               test.keepX = test.keepX,
                               folds = 2,
                               validation="loo",
                               nrepeat=2),
                "Leave-One-Out")
})


test_that("(tune.splsda:warning): catches when a class isn't represented in a given fold", {
    
    data(srbct)
    X <- srbct$gene
    Y <- srbct$class
    
    test.keepX <- c(3,6,9)
    
    expect_warning(tune.splsda(X=X, Y=Y,
                               test.keepX = test.keepX,
                               folds = 15),
                "fold")
})

