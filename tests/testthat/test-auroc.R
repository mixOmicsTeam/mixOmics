
###############################################################################
### ============================ MISSING TESTS ============================ ###
###############################################################################

# basic
## mint.block.plsda
## mint.block.splsda - created a node within the Project workflow. Once addressed, add tests

###############################################################################
### ============================ GROUND TRUTHS ============================ ###
###############################################################################

Test.Data <- readRDS(system.file("testdata", "testdata-auroc.rda", package = "mixOmics"))
Testable.Components <- Test.Data$tc
Ground.Truths <- Test.Data$gt

###############################################################################
### ================================ BASIC ================================ ###
###############################################################################



test_that("(auroc:basic): (s)plsda", {
    
    testable.components <- Testable.Components$basic
    
    data(breast.tumors)
    X <- breast.tumors$gene.exp
    Y <- breast.tumors$sample$treatment
    
    # --- PLS-DA -- # 
    GT <- Ground.Truths$basic.plsda
    
    res.plsda <- plsda(X, Y, ncomp = 2)
    
    plsda.auroc = auroc(res.plsda, roc.comp = 1, print = FALSE)
    
    invisible(capture.output(TT <- dput(plsda.auroc[testable.components])))
    
    expect_equal(TT, GT)
    
    # --- sPLS-DA -- # 
    GT <- Ground.Truths$basic.splsda
    
    choice.keepX <- c(10, 10)
    
    res.splsda <- splsda(X, Y, ncomp = 2, keepX = choice.keepX)
    
    splsda.auroc = auroc(res.splsda, roc.comp = 1, print = FALSE)
    
    invisible(capture.output(TT <- dput(splsda.auroc[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(auroc:basic): mint.(s)plsda", {
    
    testable.components <- Testable.Components$mint

    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    s <- stemcells$study
    
    # --- MINT.PLS-DA -- # 
    GT <- Ground.Truths$basic.mint.plsda

    res.mint.plsda <- mint.plsda(X, Y, ncomp = 2, study = s)

    mint.plsda.auroc = auroc(res.mint.plsda, roc.comp = 1, print = FALSE)
    
    invisible(capture.output(TT <- dput(mint.plsda.auroc[testable.components])))

    expect_equal(TT, GT)
    
    # --- MINT.sPLS-DA -- # 
    GT <- Ground.Truths$basic.mint.splsda

    choice.keepX <- c(10,10)

    res.mint.splsda <- mint.splsda(X, Y, ncomp = 2, study = s, keepX = choice.keepX)

    mint.splsda.auroc = auroc(res.mint.splsda, roc.comp = 1, print = FALSE)

    invisible(capture.output(TT <- dput(mint.splsda.auroc[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(auroc:basic): block.(s)plsda", {

    testable.components <- Testable.Components$block
        
    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.train$mirna,
             mRNA = breast.TCGA$data.train$mrna,
             proteomics = breast.TCGA$data.train$protein)
    Y = breast.TCGA$data.train$subtype
    
    # --- BLOCK.PLS-DA -- # 
    GT <- Ground.Truths$basic.block.plsda

    res.block.plsda <- block.plsda(X, Y, design = "full")

    block.plsda.auroc = auroc(res.block.plsda, roc.comp = 1, print = FALSE)

    invisible(capture.output(TT <- dput(block.plsda.auroc[testable.components])))
    
    expect_equal(TT, GT)
    
    # --- BLOCK.sPLS-DA -- # 
    GT <- Ground.Truths$basic.block.splsda

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
    
    testable.components <- Testable.Components$basic
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

    testable.components <- Testable.Components$basic
    GT <- Ground.Truths$newdata.splsda

    data(breast.tumors)
    X <- breast.tumors$gene.exp
    Y <- breast.tumors$sample$treatment
    
    d <- .minimal_train_test_subset(X, as.factor(Y), n.tr = 15, n.te = 3)

    res.plsda <- plsda(d$X.tr, d$Y.tr, ncomp = 2)

    newdata.outcome.test.auroc = auroc(res.plsda, print = FALSE, roc.comp = 1,
                                          newdata = d$X.te, outcome.test = d$Y.te)

    invisible(capture.output(TT <- dput(newdata.outcome.test.auroc[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(auroc:parameter): multilevel", {
    
    testable.components <- Testable.Components$basic
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

    testable.components <- Testable.Components$basic
    GT <- Ground.Truths$basic.plsda

    data(breast.tumors)
    X <- breast.tumors$gene.exp
    Y <- breast.tumors$sample$treatment

    res.plsda <- plsda(X, Y, ncomp = 2)

    roc.comp.auroc = auroc(res.plsda, roc.comp = 2, print = FALSE)

    invisible(capture.output(TT <- dput(roc.comp.auroc[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(auroc:parameter): roc.study", {
    
    testable.components <- Testable.Components$study
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
    
    testable.components <- Testable.Components$block
    GT <- Ground.Truths$basic.block.plsda

    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.train$mirna,
             mRNA = breast.TCGA$data.train$mrna,
             proteomics = breast.TCGA$data.train$protein)
    Y = breast.TCGA$data.train$subtype

    res.block.plsda <- block.plsda(X, Y, ncomp = 2)

    roc.block.auroc <- auroc(res.block.plsda, roc.block = 2, print = FALSE)
    
    invisible(capture.output(TT <- dput(roc.block.auroc[testable.components])))
    
    expect_equal(TT, GT)
    
    roc.block.auroc <- auroc(res.block.plsda, roc.block = "mRNA", print = FALSE)
    
    invisible(capture.output(TT <- dput(roc.block.auroc[testable.components])))
    
    expect_equal(TT, GT)
})


test_that("(auroc:parameter): study.test", {
    
    testable.components <- Testable.Components$study
    GT <- Ground.Truths$study.test.mint.splsda

    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    S <- stemcells$study
    
    d <- .minimal_train_test_subset(X, Y, S=S)

    res.mint.plsda <- mint.plsda(d$X.tr, d$Y.tr, ncomp = 2, study = d$S.tr)

    study.test.auroc = auroc(res.mint.plsda, print = FALSE,
                                    newdata = d$X.te, outcome.test = d$Y.te,
                                    study.test = d$S.te)
    
    invisible(capture.output(TT <- dput(study.test.auroc[testable.components])))
    
    expect_equal(TT, GT)
})


###############################################################################
### ================================ ERROR ================================ ###
###############################################################################


test_that("(auroc:error): catches invalid `roc.comp` values", {

    data(breast.tumors)
    X <- breast.tumors$gene.exp
    Y <- breast.tumors$sample$treatment

    res.plsda <- plsda(X, Y, ncomp = 2)

    expect_error(auroc(res.plsda, roc.comp = c(1,2), print = FALSE),
                 "roc.comp")
    
    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    S <- stemcells$study

    res.mint.plsda <- mint.plsda(X, Y, ncomp = 2, study = S)

    expect_error(auroc(res.mint.plsda, roc.comp = c(1,2), print = FALSE),
                 "roc.comp")
})


test_that("(auroc:error): catches invalid `roc.block` values", {

    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.train$mirna,
             mRNA = breast.TCGA$data.train$mrna,
             proteomics = breast.TCGA$data.train$protein)
    Y = breast.TCGA$data.train$subtype

    res.block.plsda <- block.plsda(X, Y, ncomp = 2)

    expect_error(auroc(res.block.plsda, roc.block = 4, print = FALSE),
                 "roc.block")
    
    expect_error(auroc(res.block.plsda, roc.block = TRUE, print = FALSE),
                 "roc.block")
})


test_that("(auroc:error): catches invalid `roc.study` values", {

    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    s <- stemcells$study

    res.mint.plsda <- mint.plsda(X, Y, study = s, ncomp = 2)

    expect_error(auroc(res.mint.plsda, roc.study = c(1, 2), print = FALSE),
                 "roc.study")

    expect_error(auroc(res.mint.plsda, roc.study = "study1", print = FALSE),
                 "roc.study")
})


test_that("(auroc:error): prevent sgcca and rgcca objects being used", {

    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.train$mirna,
             mRNA = breast.TCGA$data.train$mrna,
             proteomics = breast.TCGA$data.train$protein)

    res.rgcca <-  wrapper.rgcca(X)


    expect_error(auroc(res.rgcca, print = FALSE),
                 "rgcca")

    res.block.pls <- block.pls(X, indY=3)

    expect_error(auroc(res.block.pls, print = FALSE),
                 "sgcca")
})


test_that("(auroc:error): confirm 'new.data', 'outcome.test' and 'study.test' are all of same length and width", {

    data(breast.tumors)
    X <- breast.tumors$gene.exp
    Y <- breast.tumors$sample$treatment
    
    d <- .minimal_train_test_subset(X, as.factor(Y), n.tr = 15, n.te = 3)

    res.plsda <- plsda(d$X.tr, d$Y.tr, ncomp = 2)

    expect_error(auroc(res.plsda, print = FALSE, roc.comp = 1,
                       newdata = d$X.te[, 1:500], outcome.test = d$Y.te),
                 "newdata")

    expect_error(auroc(res.plsda, print = FALSE, roc.comp = 1,
                       newdata = d$X.te[1:4,], outcome.test = d$Y.te),
                 "outcome.test")
    
    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    S <- stemcells$study
    
    d <- .minimal_train_test_subset(X, Y, S=S)

    res.mint.plsda <- mint.plsda(X, Y, ncomp = 2, study = S)

    expect_error(auroc(res.mint.plsda, print = FALSE, roc.comp = 1,
                       newdata = d$X.te[1:2,], outcome.test = d$Y.te,
                       study.test = d$S.te),
                 "outcome.test")

    expect_error(auroc(res.mint.plsda, print = FALSE, roc.comp = 1,
                       newdata = d$X.te, outcome.test = d$Y.te,
                       study.test = d$S.te[1:2]),
                 "study.test")
})



