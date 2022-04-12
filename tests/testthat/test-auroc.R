
###############################################################################
### ============================ MISSING TESTS ============================ ###
###############################################################################

# basic
## mint.block.plsda
## mint.block.splsda

###############################################################################
### ============================ GROUND TRUTHS ============================ ###
###############################################################################

Test.Data <- readRDS(system.file("testdata", "testdata-auroc.rda", package = "mixOmics"))
Testable.Components <- Test.Data$tc
Ground.Truths <- Test.Data$gt

###############################################################################
### ================================ BASIC ================================ ###
###############################################################################



test_that("(auroc:basic): plsda", {
    
    testable.components <- Testable.Components$basic.plsda
    GT <- Ground.Truths$basic.plsda
    
    data(breast.tumors)
    X <- breast.tumors$gene.exp
    Y <- breast.tumors$sample$treatment
    
    res.plsda <- plsda(X, Y, ncomp = 2)
    
    plsda.auroc = auroc(res.plsda, roc.comp = 1, print = FALSE)
    
    invisible(capture.output(TT <- dput(plsda.auroc[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(auroc:basic): splsda", {
    
    testable.components <- Testable.Components$basic.splsda
    GT <- Ground.Truths$basic.splsda
    
    data(breast.tumors)
    X <- breast.tumors$gene.exp
    Y <- breast.tumors$sample$treatment
    
    choice.keepX <- c(10, 10)
    
    res.splsda <- splsda(X, Y, ncomp = 2, keepX = choice.keepX)
    
    splsda.auroc = auroc(res.splsda, roc.comp = 1, print = FALSE)
    
    invisible(capture.output(TT <- dput(splsda.auroc[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(auroc:basic): mint.plsda", {
    
    testable.components <- Testable.Components$basic.mint.plsda
    GT <- Ground.Truths$basic.mint.plsda

    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    s <- stemcells$study

    res.mint.plsda <- mint.plsda(X, Y, ncomp = 2, study = s)

    mint.plsda.auroc = auroc(res.mint.plsda, roc.comp = 1, print = FALSE)
    
    invisible(capture.output(TT <- dput(mint.plsda.auroc[testable.components])))

    expect_equal(TT, GT)
})


test_that("(auroc:basic): mint.splsda", {

    testable.components <- Testable.Components$basic.mint.splsda
    GT <- Ground.Truths$basic.mint.splsda
    
    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    s <- stemcells$study

    choice.keepX <- c(10,10)

    res.mint.splsda <- mint.splsda(X, Y, ncomp = 2, study = s, keepX = choice.keepX)

    mint.splsda.auroc = auroc(res.mint.splsda, roc.comp = 1, print = FALSE)

    invisible(capture.output(TT <- dput(mint.splsda.auroc[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(auroc:basic): block.plsda", {

    testable.components <- Testable.Components$basic.block.plsda
    GT <- Ground.Truths$basic.block.plsda
        
    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.train$mirna,
             mRNA = breast.TCGA$data.train$mrna,
             proteomics = breast.TCGA$data.train$protein)
    Y = breast.TCGA$data.train$subtype

    res.block.plsda <- block.plsda(X, Y, design = "full")

    block.plsda.auroc = auroc(res.block.plsda, roc.comp = 1, print = FALSE)

    invisible(capture.output(TT <- dput(block.plsda.auroc[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(auroc:basic): block.splsda", {
    
    testable.components <- Testable.Components$basic.block.splsda
    GT <- Ground.Truths$basic.block.splsda

    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.train$mirna,
             mRNA = breast.TCGA$data.train$mrna,
             proteomics = breast.TCGA$data.train$protein)
    Y = breast.TCGA$data.train$subtype

    choice.keepX <- list(miRNA=c(10,10),
                         mRNA=c(10,10),
                         proteomics=c(10,10))

    res.block.splsda <- block.splsda(X, Y, design = "full", keepX = choice.keepX)

    block.splsda.auroc = auroc(res.block.splsda, roc.comp = 1, print = FALSE)

    invisible(capture.output(TT <- dput(block.splsda.auroc[testable.components])))
    
    expect_equal(TT, GT)
})


###############################################################################
### ================================ DATA ================================= ###
###############################################################################


test_that("(auroc:data): splsda, srbct", {
    
    testable.components <- Testable.Components$srbct.splsda
    GT <- Ground.Truths$srbct.splsda

    data(srbct)
    X <- srbct$gene
    Y <- srbct$class

    choice.keepX <- c(10, 10)

    res.splsda <- splsda(X, Y, ncomp = 2, keepX = choice.keepX)

    srbct.splsda.auroc = auroc(res.splsda, roc.comp = 1, print = FALSE)

    invisible(capture.output(TT <- dput(srbct.splsda.auroc[testable.components])))
    
    expect_equal(TT, GT)
})


###############################################################################
### ============================= PARAMETER =============================== ###
###############################################################################


test_that("(auroc:parameter): newdata/outcome.test", {

    testable.components <- Testable.Components$newdata.splsda
    GT <- Ground.Truths$newdata.splsda

    data(breast.tumors)
    set.seed(3)
    test=sample(1:47,5,replace=FALSE)

    X <- breast.tumors$gene.exp[-test,]
    Y <- breast.tumors$sample$treatment[-test]

    X.test<-breast.tumors$gene.exp[test,]
    Y.test<-breast.tumors$sample$treatment[test]

    res.plsda <- plsda(X, Y, ncomp = 2)

    newdata.outcome.test.auroc = auroc(res.plsda, print = FALSE, roc.comp = 1,
                                          newdata = X.test, outcome.test = Y.test)

    invisible(capture.output(TT <- dput(newdata.outcome.test.auroc[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(auroc:parameter): multilevel", {
    
    testable.components <- Testable.Components$multilevel.splsda
    GT <- Ground.Truths$multilevel.splsda

    data(diverse.16S)
    X <- diverse.16S$data.TSS
    Y <- diverse.16S$bodysite
    mL <- diverse.16S$sample

    res.plsda <- plsda(X, Y, ncomp = 2)

    multilevel.auroc <- auroc(res.plsda, print = FALSE, roc.comp = 2,
                                 multilevel = mL)

    invisible(capture.output(TT <- dput(multilevel.auroc[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(auroc:parameter): roc.comp", {

    testable.components <- Testable.Components$roc.comp.splsda
    GT <- Ground.Truths$roc.comp.splsda

    data(breast.tumors)
    X <- breast.tumors$gene.exp
    Y <- breast.tumors$sample$treatment

    res.plsda <- plsda(X, Y, ncomp = 2)

    roc.comp.auroc = auroc(res.plsda, roc.comp = 2, print = FALSE)

    invisible(capture.output(TT <- dput(roc.comp.auroc[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(auroc:parameter): roc.study", {
    
    testable.components <- Testable.Components$roc.study.mint.splsda
    GT <- Ground.Truths$roc.study.mint.splsda

    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    s <- stemcells$study

    res.mint.plsda <- mint.plsda(X, Y, ncomp = 2, study = s)

    roc.study.mint.splsda = auroc(res.mint.plsda, roc.study = 2, print = FALSE)

    invisible(capture.output(TT <- dput(roc.study.mint.splsda[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(auroc:parameter): roc.block", {
    
    testable.components <- Testable.Components$roc.block.splsda
    GT <- Ground.Truths$roc.block.splsda

    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.train$mirna,
             mRNA = breast.TCGA$data.train$mrna,
             proteomics = breast.TCGA$data.train$protein)
    Y = breast.TCGA$data.train$subtype

    res.block.plsda <- block.plsda(X, Y, ncomp = 2)

    roc.block.auroc <- auroc(res.block.plsda, roc.block = 2, print = FALSE)
    
    invisible(capture.output(TT <- dput(roc.block.auroc[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(auroc:parameter): study.test", {
    
    testable.components <- Testable.Components$study.test.mint.splsda
    GT <- Ground.Truths$study.test.mint.splsda

    data(stemcells)
    set.seed(3)
    test=sample(1:125,5,replace=FALSE)


    X.tr <- stemcells$gene[-test,]
    X.te <- stemcells$gene[test,]
    Y.tr <- stemcells$celltype[-test]
    Y.te <- stemcells$celltype[test]
    s.tr <- stemcells$study[-test]
    s.te <- stemcells$study[test]

    res.mint.plsda <- mint.plsda(X.tr, Y.tr, ncomp = 2, study = s.tr)

    roc.study.test.auroc = auroc(res.mint.plsda, print = FALSE,
                                    newdata = X.te, outcome.test = Y.te,
                                    study.test = s.te)
    
    invisible(capture.output(TT <- dput(roc.study.test.auroc[testable.components])))
    
    expect_equal(TT, GT)
})


###############################################################################
### ================================ ERROR ================================ ###
###############################################################################


test_that("(auroc:error): roc.block not greater than number of blocks", {

    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.train$mirna,
             mRNA = breast.TCGA$data.train$mrna,
             proteomics = breast.TCGA$data.train$protein)
    Y = breast.TCGA$data.train$subtype

    res.block.plsda <- block.plsda(X, Y, ncomp = 2)

    expect_error(auroc(res.block.plsda, roc.block = 4, print = FALSE),
                 "roc.block cannot be greater than 3",
                 fixed = TRUE)
})


test_that("(auroc:error): roc.study not an invalid value", {

    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    s <- stemcells$study

    res.mint.plsda <- mint.plsda(X, Y, study = s, ncomp = 2)

    expect_error(auroc(res.mint.plsda, roc.study = c(1, 2), print = FALSE),
                 "`roc.study' must be a single entry,
    either `global' or one of levels(object$study)",
                 fixed = TRUE)

    expect_error(auroc(res.mint.plsda, roc.study = "study1", print = FALSE),
                 "'roc.study' must be one of 'levels(object$study)'",
                 fixed = TRUE)
})


test_that("(auroc:error): prevent sgcca and rgcca objects being used", {

    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.train$mirna,
             mRNA = breast.TCGA$data.train$mrna,
             proteomics = breast.TCGA$data.train$protein)

    res.rgcca <-  wrapper.rgcca(X)


    expect_error(auroc(res.rgcca, print = FALSE),
                 "no applicable method for 'auroc' applied to an object of class \"c('sparse.rgcca', 'rgcca')",
                 fixed = TRUE)

    res.block.pls <- block.pls(X, indY=3)

    expect_error(auroc(res.block.pls, print = FALSE),
                 "no applicable method for 'auroc' applied to an object of class \"c('block.pls', 'sgcca')",
                 fixed = TRUE)
})


test_that("(auroc:error): ensure new.data contains all the same features", {

    data(breast.tumors)
    set.seed(3)
    test=sample(1:47,5,replace=FALSE)

    X <- breast.tumors$gene.exp[-test,]
    Y <- breast.tumors$sample$treatment[-test]

    X.test<-breast.tumors$gene.exp[test,1:500]
    Y.test<-breast.tumors$sample$treatment[test]

    res.plsda <- plsda(X, Y, ncomp = 2)

    expect_error(auroc(res.plsda, print = FALSE, roc.comp = 1,
                       newdata = X.test, outcome.test = Y.test),
                 "newdata' must include all the variables of 'object$X",
                 fixed = TRUE)
})



