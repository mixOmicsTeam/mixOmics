
###############################################################################
### ============================ MISSING TESTS ============================ ###
###############################################################################

# parameter
## validation - commented out currently

###############################################################################
### ============================ GROUND TRUTHS ============================ ###
###############################################################################

Test.Data <- readRDS(system.file("testdata", "testdata-tune.block.splsda.rda", package = "mixOmics"))
Testable.Components <- Test.Data$tc
Ground.Truths <- Test.Data$gt

###############################################################################
### ================================ BASIC ================================ ###
###############################################################################


test_that("(tune.block.splsda:basic): basic", {
    
    testable.components <- Testable.Components$basic.tune.block.splsda
    GT <- Ground.Truths$basic.tune.block.splsda
    
    data(breast.TCGA)
    X <- list(miRNA = breast.TCGA$data.train$mirna[,1:10], 
                mRNA = breast.TCGA$data.train$mrna[,1:10], 
                proteomics = breast.TCGA$data.train$protein[,1:10])
    Y <- breast.TCGA$data.train$subtype
    
    d <- .minimal_train_test_subset(X, Y)
    
    test.keepX <- list(miRNA = c(6,9), 
                       mRNA = c(6,9), 
                       proteomics = c(6,9))
    
    block.splsda.tune <- suppressWarnings(suppressMessages(tune.block.splsda(d$X.tr, d$Y.tr,
                                                           test.keepX = test.keepX, 
                                                           folds=3, design = "full")))
    
    invisible(capture.output(TT <- dput(block.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


###############################################################################
### ================================ DATA ================================= ###
###############################################################################


test_that("(tune.block.splsda:data): breast.test", {
    
    testable.components <- Testable.Components$breast.test.tune.block.splsda
    GT <- Ground.Truths$breast.test.tune.block.splsda
    
    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.test$mirna[,1:10],
             mRNA = breast.TCGA$data.test$mrna[,1:10])
    Y = breast.TCGA$data.test$subtype

    d <- .minimal_train_test_subset(X, Y)

    test.keepX <- list(miRNA = c(6,9),
                       mRNA = c(6,9))
    
    breast.test.block.splsda.tune <- suppressWarnings(suppressMessages(tune.block.splsda(d$X.tr, d$Y.tr,
                                                                       test.keepX = test.keepX, 
                                                                       folds=3, design = "full")))
    
    invisible(capture.output(TT <- dput(breast.test.block.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


###############################################################################
### ============================== PARAMETER ============================== ###
###############################################################################


# test_that("(tune.block.splsda:parameter): validation", {
# 
#     testable.components <- Testable.Components$validation.tune.block.splsda
#     GT <- Ground.Truths$validation.tune.block.splsda
# 
#     data(breast.TCGA)
#     X <- list(miRNA = breast.TCGA$data.train$mirna[,1:10],
#                 mRNA = breast.TCGA$data.train$mrna[,1:10],
#                 proteomics = breast.TCGA$data.train$protein[,1:10])
#     Y <- breast.TCGA$data.train$subtype
# 
#     d <- .minimal_train_test_subset(X, Y)
# 
#     test.keepX <- list(miRNA = c(6,9),
#                        mRNA = c(6,9),
#                        proteomics = c(6,9))
# 
#     breast.test.block.splsda.tune <- suppressWarnings(suppressMessages(tune.block.splsda(d$X.tr, d$Y.tr,
#                                                                        test.keepX = test.keepX,
#                                                                        validation = "loo", design = "full")))
# 
#     invisible(capture.output(TT <- dput(breast.test.block.splsda.tune[testable.components])))
# 
#     expect_equal(TT, GT)
# })


test_that("(tune.block.splsda:parameter): dist", {

    testable.components <- Testable.Components$dist.tune.block.splsda
    GT <- Ground.Truths$dist.tune.block.splsda

    data(breast.TCGA)
    X <- list(miRNA = breast.TCGA$data.train$mirna[,1:10],
                mRNA = breast.TCGA$data.train$mrna[,1:10],
                proteomics = breast.TCGA$data.train$protein[,1:10])
    Y <- breast.TCGA$data.train$subtype

    d <- .minimal_train_test_subset(X, Y)

    test.keepX <- list(miRNA = c(6,9),
                       mRNA = c(6,9),
                       proteomics = c(6,9))

    dist.block.splsda.tune <- suppressWarnings(suppressMessages(tune.block.splsda(d$X.tr, d$Y.tr,
                                                                       test.keepX = test.keepX,
                                                                       folds = 3, design = "full",
                                                                       dist="centroids.dist")))

    invisible(capture.output(TT <- dput(dist.block.splsda.tune[testable.components])))

    expect_equal(TT, GT)
})


test_that("(tune.block.splsda:parameter): measure", {

    testable.components <- Testable.Components$basic.tune.block.splsda
    GT <- Ground.Truths$basic.tune.block.splsda

    data(breast.TCGA)
    X <- list(miRNA = breast.TCGA$data.train$mirna[,1:10],
                mRNA = breast.TCGA$data.train$mrna[,1:10],
                proteomics = breast.TCGA$data.train$protein[,1:10])
    Y <- breast.TCGA$data.train$subtype

    d <- .minimal_train_test_subset(X, Y)

    test.keepX <- list(miRNA = c(6,9),
                       mRNA = c(6,9),
                       proteomics = c(6,9))

    measure.block.splsda.tune <- suppressWarnings(suppressMessages(tune.block.splsda(d$X.tr, d$Y.tr,
                                                                       test.keepX = test.keepX,
                                                                       folds = 3, design = "full",
                                                                       measure="overall")))

    invisible(capture.output(TT <- dput(measure.block.splsda.tune[testable.components])))

    expect_equal(TT, GT)
})


test_that("(tune.block.splsda:parameter): BPPARAM", {

    testable.components <- Testable.Components$BPPARAM.tune.block.splsda
    GT <- Ground.Truths$BPPARAM.tune.block.splsda

    library(BiocParallel)

    data(breast.TCGA)
    X <- list(miRNA = breast.TCGA$data.train$mirna[,1:10],
                mRNA = breast.TCGA$data.train$mrna[,1:10],
                proteomics = breast.TCGA$data.train$protein[,1:10])
    Y <- breast.TCGA$data.train$subtype

    d <- .minimal_train_test_subset(X, Y)

    test.keepX <- list(miRNA = c(6,9),
                       mRNA = c(6,9),
                       proteomics = c(6,9))

    BPPARAM.block.splsda.tune <- suppressWarnings(suppressMessages(tune.block.splsda(d$X.tr, d$Y.tr,
                                                                       test.keepX = test.keepX,
                                                                       folds = 3, design = "full",
                                                                       BPPARAM=SnowParam(workers=4))))

    invisible(capture.output(TT <- dput(BPPARAM.block.splsda.tune[testable.components])))

    expect_equal(TT, GT)
})


test_that("(tune.block.splsda:parameter): design", {

    testable.components <- Testable.Components$design.tune.block.splsda
    GT <- Ground.Truths$design.tune.block.splsda

    data(breast.TCGA)
    X <- list(miRNA = breast.TCGA$data.train$mirna[,1:10],
                mRNA = breast.TCGA$data.train$mrna[,1:10],
                proteomics = breast.TCGA$data.train$protein[,1:10])
    Y <- breast.TCGA$data.train$subtype

    d <- .minimal_train_test_subset(X, Y)

    test.keepX <- list(miRNA = c(6,9),
                       mRNA = c(6,9),
                       proteomics = c(6,9))

    design.block.splsda.tune <- suppressWarnings(suppressMessages(tune.block.splsda(d$X.tr, d$Y.tr,
                                                                       test.keepX = test.keepX,
                                                                       folds = 3, design = 0.1)))
    
    invisible(capture.output(TT <- dput(design.block.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(tune.block.splsda:parameter): folds", {

    testable.components <- Testable.Components$folds.tune.block.splsda
    GT <- Ground.Truths$folds.tune.block.splsda

    data(breast.TCGA)
    X <- list(miRNA = breast.TCGA$data.train$mirna[,1:10],
                mRNA = breast.TCGA$data.train$mrna[,1:10],
                proteomics = breast.TCGA$data.train$protein[,1:10])
    Y <- breast.TCGA$data.train$subtype

    d <- .minimal_train_test_subset(X, Y)

    test.keepX <- list(miRNA = c(6,9),
                       mRNA = c(6,9),
                       proteomics = c(6,9))

    folds.block.splsda.tune <- suppressWarnings(suppressMessages(tune.block.splsda(d$X.tr, d$Y.tr,
                                                                       test.keepX = test.keepX,
                                                                       folds = 2, design="full")))
    
    invisible(capture.output(TT <- dput(folds.block.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(tune.block.splsda:parameter): near.zero.var", {

    testable.components <- Testable.Components$near.zero.var.tune.block.splsda
    GT <- Ground.Truths$near.zero.var.tune.block.splsda

    data(breast.TCGA)
    X <- list(miRNA = breast.TCGA$data.train$mirna[,1:10],
                mRNA = breast.TCGA$data.train$mrna[,1:10],
                proteomics = breast.TCGA$data.train$protein[,1:10])
    Y <- breast.TCGA$data.train$subtype
    
    X$miRNA[1:145,1:2] <- 0

    d <- .minimal_train_test_subset(X, Y)

    test.keepX <- list(miRNA = c(6,9),
                       mRNA = c(6,9),
                       proteomics = c(6,9))

    near.zero.var.block.splsda.tune <- suppressWarnings(suppressMessages(tune.block.splsda(d$X.tr, d$Y.tr,
                                                                       test.keepX = test.keepX,
                                                                       folds = 3, design="full",
                                                                       near.zero.var = T)))
    
    invisible(capture.output(TT <- dput(near.zero.var.block.splsda.tune[testable.components])))
    
    expect_equal(TT, GT)
})


###############################################################################
### ================================ ERROR ================================ ###
###############################################################################


test_that("(tune.block.splsda:error): catches invalid use of 'cpus'", {
    
    data(breast.TCGA)
    X <- list(miRNA = breast.TCGA$data.train$mirna[,1:10], 
                mRNA = breast.TCGA$data.train$mrna[,1:10], 
                proteomics = breast.TCGA$data.train$protein[,1:10])
    Y <- breast.TCGA$data.train$subtype
    
    d <- .minimal_train_test_subset(X, Y)
    
    test.keepX <- list(miRNA = c(6,9), 
                       mRNA = c(6,9), 
                       proteomics = c(6,9))
    
    expect_error(suppressWarnings(suppressMessages(tune.block.splsda(d$X.tr, d$Y.tr,
                                                           test.keepX = test.keepX, 
                                                           folds=3, design = "full",
                                                           cpus = 4))),
                 "cpus")
})


test_that("(tune.block.splsda:error): catches invalid 'Y' and 'indY' objects", {
    
    data(breast.TCGA)
    X <- list(miRNA = breast.TCGA$data.train$mirna[,1:10], 
                mRNA = breast.TCGA$data.train$mrna[,1:10], 
                proteomics = breast.TCGA$data.train$protein[,1:10])
    Y <- breast.TCGA$data.train$subtype
    
    d <- .minimal_train_test_subset(X, Y)
    
    test.keepX <- list(miRNA = c(6,9), 
                       mRNA = c(6,9), 
                       proteomics = c(6,9))
    
    d$Y.tr <- matrix(rnorm(6))
    
    expect_error(suppressWarnings(suppressMessages(tune.block.splsda(d$X.tr, d$Y.tr,
                                                           test.keepX = test.keepX, 
                                                           folds=3, design = "full"))),
                 "Y")
    
    d$Y.tr <- factor(rep("LumA", 6))
    
    expect_error(suppressWarnings(suppressMessages(tune.block.splsda(d$X.tr, d$Y.tr,
                                                           test.keepX = test.keepX, 
                                                           folds=3, design = "full"))),
                 "Y")
    
    samples <- c(1:2, 50:51, 79:80)
    X <- list(miRNA = breast.TCGA$data.train$mirna[samples,1:10], 
              mRNA = breast.TCGA$data.train$mrna[samples,1:10], 
              proteomics = breast.TCGA$data.train$protein[samples,1:10],
              Y = factor(rep("LumA", 6)))
    
    expect_error(suppressWarnings(suppressMessages(tune.block.splsda(X, indY=3,
                                                           test.keepX = test.keepX, 
                                                           folds=3, design = "full"))),
                 "Y")
    
    expect_error(suppressWarnings(suppressMessages(tune.block.splsda(X, indY=4,
                                                           test.keepX = test.keepX, 
                                                           folds=3, design = "full"))),
                 "Y")
    
    expect_error(suppressWarnings(suppressMessages(tune.block.splsda(X,
                                                           test.keepX = test.keepX, 
                                                           folds=3, design = "full"))),
                 "Y")

    
})


test_that("(tune.block.splsda:error): catches invalid value of 'progressBar'", {
    
    data(breast.TCGA)
    X <- list(miRNA = breast.TCGA$data.train$mirna[,1:10], 
                mRNA = breast.TCGA$data.train$mrna[,1:10], 
                proteomics = breast.TCGA$data.train$protein[,1:10])
    Y <- breast.TCGA$data.train$subtype
    
    d <- .minimal_train_test_subset(X, Y)
    
    test.keepX <- list(miRNA = c(6,9), 
                       mRNA = c(6,9), 
                       proteomics = c(6,9))
    
    expect_error(suppressWarnings(suppressMessages(tune.block.splsda(d$X.tr, d$Y.tr,
                                                           test.keepX = test.keepX, 
                                                           folds=3, design = "full",
                                                           progressBar = "non-logical"))),
                 "progressBar")
})


test_that("(tune.block.splsda:error): catches invalid value of 'ncomp'", {
    
    data(breast.TCGA)
    X <- list(miRNA = breast.TCGA$data.train$mirna[,1:10], 
                mRNA = breast.TCGA$data.train$mrna[,1:10], 
                proteomics = breast.TCGA$data.train$protein[,1:10])
    Y <- breast.TCGA$data.train$subtype
    
    d <- .minimal_train_test_subset(X, Y)
    
    test.keepX <- list(miRNA = c(6,9), 
                       mRNA = c(6,9), 
                       proteomics = c(6,9))
    
    expect_error(suppressWarnings(suppressMessages(tune.block.splsda(d$X.tr, d$Y.tr,
                                                                     ncomp = -1,
                                                           test.keepX = test.keepX, 
                                                           folds=3, design = "full"))),
                 "ncomp")
    
    expect_error(suppressWarnings(suppressMessages(tune.block.splsda(d$X.tr, d$Y.tr,
                                                                     ncomp = NULL,
                                                           test.keepX = test.keepX, 
                                                           folds=3, design = "full"))),
                 "ncomp")
    
    expect_error(suppressWarnings(suppressMessages(tune.block.splsda(d$X.tr, d$Y.tr,
                                                                     ncomp = "non-numeric",
                                                           test.keepX = test.keepX, 
                                                           folds=3, design = "full"))),
                 "ncomp")
})


test_that("(tune.block.splsda:error): catches invalid value of 'validation'", {
    
    data(breast.TCGA)
    X <- list(miRNA = breast.TCGA$data.train$mirna[,1:10], 
                mRNA = breast.TCGA$data.train$mrna[,1:10], 
                proteomics = breast.TCGA$data.train$protein[,1:10])
    Y <- breast.TCGA$data.train$subtype
    
    d <- .minimal_train_test_subset(X, Y)
    
    test.keepX <- list(miRNA = c(6,9), 
                       mRNA = c(6,9), 
                       proteomics = c(6,9))
    
    expect_error(suppressWarnings(suppressMessages(tune.block.splsda(d$X.tr, d$Y.tr,
                                                           test.keepX = test.keepX, 
                                                           folds=3, design = "full",
                                                           validation = "random-value"))),
                 "validation")
    
    expect_error(suppressWarnings(suppressMessages(tune.block.splsda(d$X.tr, d$Y.tr,
                                                           test.keepX = test.keepX, 
                                                           folds=3, design = "full",
                                                           validation = NA))),
                 "validation")
})


test_that("(tune.block.splsda:error): catches invalid value of 'measure'", {
    
    data(breast.TCGA)
    X <- list(miRNA = breast.TCGA$data.train$mirna[,1:10], 
                mRNA = breast.TCGA$data.train$mrna[,1:10], 
                proteomics = breast.TCGA$data.train$protein[,1:10])
    Y <- breast.TCGA$data.train$subtype
    
    d <- .minimal_train_test_subset(X, Y)
    
    test.keepX <- list(miRNA = c(6,9), 
                       mRNA = c(6,9), 
                       proteomics = c(6,9))
    
    expect_error(suppressWarnings(suppressMessages(tune.block.splsda(d$X.tr, d$Y.tr,
                                                           test.keepX = test.keepX, 
                                                           folds=3, design = "full",
                                                           measure = "random-value"))),
                 "measure")
})


test_that("(tune.block.splsda:error): catches invalid value of 'already.tested.X'", {
    
    data(breast.TCGA)
    X <- list(miRNA = breast.TCGA$data.train$mirna[,1:10], 
                mRNA = breast.TCGA$data.train$mrna[,1:10], 
                proteomics = breast.TCGA$data.train$protein[,1:10])
    Y <- breast.TCGA$data.train$subtype
    
    d <- .minimal_train_test_subset(X, Y)
    
    test.keepX <- list(miRNA = c(6,9), 
                       mRNA = c(6,9), 
                       proteomics = c(6,9))
    
    expect_error(suppressWarnings(suppressMessages(tune.block.splsda(d$X.tr, d$Y.tr,
                                                           test.keepX = test.keepX, 
                                                           folds=3, design = "full",
                                                           already.tested.X = "random-value"))),
                 "already.tested.X")
    
    expect_error(suppressWarnings(suppressMessages(tune.block.splsda(d$X.tr, d$Y.tr,
                                                           test.keepX = test.keepX, 
                                                           folds=3, design = "full",
                                                           already.tested.X = NULL))),
                 "already.tested.X")
    
    already.tested.X <- list(miRNA = c(6,6), 
                             mRNA = c(9,6), 
                             proteomics = c(6,6))
    
    expect_error(suppressWarnings(suppressMessages(tune.block.splsda(d$X.tr, d$Y.tr,
                                                           test.keepX = test.keepX, 
                                                           folds=3, design = "full",
                                                           already.tested.X = already.tested.X))),
                 "already.tested.X")
    
    already.tested.X <- list(miRNA = c(6), 
                             mRNA = c(9,6), 
                             proteomics = c(6,6))
    
    expect_error(suppressWarnings(suppressMessages(tune.block.splsda(d$X.tr, d$Y.tr,
                                                                     ncomp=3,
                                                           test.keepX = test.keepX, 
                                                           folds=3, design = "full",
                                                           already.tested.X = already.tested.X))),
                 "already.tested.X")
    
    already.tested.X <- list(miRNA = list(6,6), 
                             mRNA = c(6,6), 
                             proteomics = c(6,6))
    
    expect_error(suppressWarnings(suppressMessages(tune.block.splsda(d$X.tr, d$Y.tr,
                                                                     ncomp=3,
                                                           test.keepX = test.keepX, 
                                                           folds=3, design = "full",
                                                           already.tested.X = already.tested.X))),
                 "already.tested.X")
})


test_that("(tune.block.splsda:error): catches invalid value of 'test.keepX'", {
    
    data(breast.TCGA)
    X <- list(miRNA = breast.TCGA$data.train$mirna[,1:10], 
                mRNA = breast.TCGA$data.train$mrna[,1:10], 
                proteomics = breast.TCGA$data.train$protein[,1:10])
    Y <- breast.TCGA$data.train$subtype
    
    d <- .minimal_train_test_subset(X, Y)
    
    test.keepX <- list(miRNA = c(6,9), 
                       mRNA = c(6,9))
    
    expect_error(suppressWarnings(suppressMessages(tune.block.splsda(d$X.tr, d$Y.tr,
                                                           test.keepX = test.keepX, 
                                                           folds=3, design = "full"))),
                 "test.keepX")
})


test_that("(tune.block.splsda:error): catches edge cases when utilising 'near.zero.var'", {

    data(breast.TCGA)
    X <- list(miRNA = breast.TCGA$data.train$mirna[,1:10],
                mRNA = breast.TCGA$data.train$mrna[,1:10],
                proteomics = breast.TCGA$data.train$protein[,1:10])
    Y <- breast.TCGA$data.train$subtype
    
    X$miRNA[1:145,1:2] <- 0

    d <- .minimal_train_test_subset(X, Y)

    test.keepX <- list(miRNA = c(6,9),
                       mRNA = c(6,9),
                       proteomics = c(6,9))

    expect_warning(suppressMessages(tune.block.splsda(d$X.tr, d$Y.tr,
                                                       test.keepX = test.keepX,
                                                       folds = 3, design="full",
                                                       near.zero.var = T)),
                   "near-zero variance")
    
    X$miRNA[1:145,1:10] <- 0

    d <- .minimal_train_test_subset(X, Y)
    
    expect_error(suppressWarnings(suppressMessages(tune.block.splsda(d$X.tr, d$Y.tr,
                                                       test.keepX = test.keepX,
                                                       folds = 3, design="full",
                                                       near.zero.var = T))),
                   "No more variables")
})


test_that("(tune.block.splsda:error): catches invalid value of 'folds'", {
    
    data(breast.TCGA)
    X <- list(miRNA = breast.TCGA$data.train$mirna[,1:10], 
                mRNA = breast.TCGA$data.train$mrna[,1:10], 
                proteomics = breast.TCGA$data.train$protein[,1:10])
    Y <- breast.TCGA$data.train$subtype
    
    d <- .minimal_train_test_subset(X, Y)
    
    test.keepX <- list(miRNA = c(6,9),
                       mRNA = c(6,9),
                       proteomics = c(6,9))
    
    expect_error(suppressWarnings(suppressMessages(tune.block.splsda(d$X.tr, d$Y.tr,
                                                           test.keepX = test.keepX, 
                                                           folds=NULL, design = "full"))),
                 "folds")
    
    expect_error(suppressWarnings(suppressMessages(tune.block.splsda(d$X.tr, d$Y.tr,
                                                           test.keepX = test.keepX, 
                                                           folds="random-value", design = "full"))),
                 "folds")
    
    expect_error(suppressWarnings(suppressMessages(tune.block.splsda(d$X.tr, d$Y.tr,
                                                           test.keepX = test.keepX, 
                                                           folds=1, design = "full"))),
                 "folds")
    
    expect_error(suppressWarnings(suppressMessages(tune.block.splsda(d$X.tr, d$Y.tr,
                                                           test.keepX = test.keepX, 
                                                           folds=200, design = "full"))),
                 "folds")
})


###############################################################################
### =============================== WARNINGS ============================== ###
###############################################################################


test_that("(tune.block.splsda:warning): catches nrepeat>1 when validation='loo'", {
    
    data(breast.TCGA)
    X <- list(miRNA = breast.TCGA$data.train$mirna[,1:10], 
                mRNA = breast.TCGA$data.train$mrna[,1:10], 
                proteomics = breast.TCGA$data.train$protein[,1:10])
    Y <- breast.TCGA$data.train$subtype
    
    d <- .minimal_train_test_subset(X, Y)
    
    test.keepX <- list(miRNA = c(6,9), 
                       mRNA = c(6,9), 
                       proteomics = c(6,9))
    
    expect_warning(suppressMessages(tune.block.splsda(d$X.tr, d$Y.tr,
                                                           test.keepX = test.keepX, 
                                                           folds=3, design = "full",
                                                           validation = "loo",
                                                           nrepeat=2)),
                 "validation")
})
