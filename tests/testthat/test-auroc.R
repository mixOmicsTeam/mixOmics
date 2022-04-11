
###############################################################################
### ============================ MISSING TESTS ============================ ###
###############################################################################

# basic
## mint.block.plsda
## mint.block.splsda

###############################################################################
### ============================ GROUND TRUTHS ============================ ###
###############################################################################

Testable.Components <- list(basic.plsda = c("Comp1", "Comp2"),
                            basic.splsda = c("Comp1", "Comp2"),
                            basic.mint.plsda = c("Comp1"),
                            basic.mint.splsda = c("Comp1"),
                            basic.block.plsda = c("miRNA", "mRNA", "proteomics"),
                            basic.block.splsda = c("miRNA", "mRNA", "proteomics"),
                            srbct.splsda = c("Comp1", "Comp2"),
                            newdata.splsda = c("Comp1", "Comp2"),
                            multilevel.splsda = c("Comp1", "Comp2"),
                            roc.comp.splsda = c("Comp1", "Comp2"),
                            roc.study.mint.splsda = c("Comp1"),
                            roc.block.splsda = c("miRNA", "mRNA", "proteomics"),
                            study.test.mint.splsda = c("Comp2")
)

Ground.Truths <- list(basic.plsda = list(Comp1 = structure(c(0.863, 2.473e-05), 
                                                           .Dim = 1:2, 
                                                           .Dimnames = list("AF vs BE", 
                                                                            c("AUC", "p-value"))), 
                                         Comp2 = structure(c(0.9981, 7.124e-09), 
                                                           .Dim = 1:2, 
                                                           .Dimnames = list("AF vs BE", c("AUC", "p-value")))),
                      basic.splsda = list(Comp1 = structure(c(0.9704, 4.623e-08), 
                                                            .Dim = 1:2, 
                                                            .Dimnames = list("AF vs BE", c("AUC", "p-value"))), 
                                          Comp2 = structure(c(0.9981, 7.124e-09), 
                                                            .Dim = 1:2, 
                                                            .Dimnames = list("AF vs BE", c("AUC", "p-value")))),
                      basic.mint.plsda = list(Comp1 = structure(c(0.993, 0.6361, 0.8258, 4.441e-16, 0.01658, 3.671e-10), 
                                                                .Dim = 3:2, 
                                                                .Dimnames = list(c("Fibroblast vs Other(s)", 
                                                                                   "hESC vs Other(s)", 
                                                                                   "hiPSC vs Other(s)"), 
                                                                                 c("AUC", "p-value")))),
                      basic.mint.splsda = list(Comp1 = structure(c(1, 0.7844, 0.8075, 2.22e-16, 5.496e-07, 3.3e-09), 
                                                                 .Dim = 3:2, 
                                                                 .Dimnames = list(c("Fibroblast vs Other(s)", 
                                                                                    "hESC vs Other(s)",
                                                                                    "hiPSC vs Other(s)"), 
                                                                                  c("AUC", "p-value")))),
                      basic.block.plsda = list(miRNA = list(comp1 = structure(c(0.9365, 0.5861, 0.9218, 0, 0.1453, 0), 
                                                                              .Dim = 3:2, 
                                                                              .Dimnames = list(c("Basal vs Other(s)", 
                                                                                                 "Her2 vs Other(s)", 
                                                                                                 "LumA vs Other(s)"), 
                                                                                               c("AUC", "p-value"))), 
                                                            comp2 = structure(c(0.9494, 0.5856, 0.9257, 0, 0.1479, 0), 
                                                                              .Dim = 3:2, 
                                                                              .Dimnames = list(c("Basal vs Other(s)", 
                                                                                                 "Her2 vs Other(s)", 
                                                                                                 "LumA vs Other(s)"), 
                                                                                               c("AUC", "p-value")))), 
                                               mRNA = list(comp1 = structure(c(0.9937, 0.6206, 0.9918, 0, 0.04144, 0), 
                                                                             .Dim = 3:2, 
                                                                             .Dimnames = list(c("Basal vs Other(s)", 
                                                                                                "Her2 vs Other(s)", 
                                                                                                "LumA vs Other(s)"), 
                                                                                              c("AUC", "p-value"))), 
                                                           comp2 = structure(c(0.9962, 0.5725, 0.9893, 0, 0.2201, 0), 
                                                                             .Dim = 3:2, 
                                                                             .Dimnames = list(c("Basal vs Other(s)", 
                                                                                                "Her2 vs Other(s)", 
                                                                                                "LumA vs Other(s)"), 
                                                                                              c("AUC", "p-value")))), 
                                               proteomics = list(comp1 = structure(c(0.9564, 0.6578, 0.9844, 0, 0.007614, 0), 
                                                                                   .Dim = 3:2, 
                                                                                   .Dimnames = list(c("Basal vs Other(s)", 
                                                                                                      "Her2 vs Other(s)", 
                                                                                                      "LumA vs Other(s)"), 
                                                                                                    c("AUC", "p-value" ))), 
                                                                 comp2 = structure(c(0.9757, 0.6564, 0.9865, 0, 0.008164,0), 
                                                                                   .Dim = 3:2, 
                                                                                   .Dimnames = list(c("Basal vs Other(s)", 
                                                                                                      "Her2 vs Other(s)", 
                                                                                                      "LumA vs Other(s)"), 
                                                                                                    c("AUC", "p-value"))))),
                      basic.block.splsda = list(miRNA = list(comp1 = structure(c(0.9549, 0.5642, 0.9232, 0, 0.2778, 0), 
                                                                               .Dim = 3:2, 
                                                                               .Dimnames = list(c("Basal vs Other(s)",
                                                                                                  "Her2 vs Other(s)", 
                                                                                                  "LumA vs Other(s)"), 
                                                                                                c("AUC", "p-value"))), 
                                                             comp2 = structure(c(0.9619, 0.8644, 0.9595, 0, 7.078e-10, 0), 
                                                                               .Dim = 3:2, 
                                                                               .Dimnames = list(c("Basal vs Other(s)", 
                                                                                                  "Her2 vs Other(s)", 
                                                                                                  "LumA vs Other(s)"), 
                                                                                                c("AUC", "p-value")))), 
                                                mRNA = list(comp1 = structure(c(0.9968, 0.6083, 0.9867, 0, 0.06689, 0), 
                                                                              .Dim = 3:2, 
                                                                              .Dimnames = list(c("Basal vs Other(s)", 
                                                                                                 "Her2 vs Other(s)", 
                                                                                                 "LumA vs Other(s)"), 
                                                                                               c("AUC", "p-value"))), 
                                                            comp2 = structure(c(0.9983, 0.9631, 0.9954, 0, 4.885e-15, 0), 
                                                                              .Dim = 3:2, 
                                                                              .Dimnames = list(c("Basal vs Other(s)", 
                                                                                                 "Her2 vs Other(s)", 
                                                                                                 "LumA vs Other(s)"), 
                                                                                               c("AUC", "p-value")))), 
                                                proteomics = list(comp1 = structure(c(0.9729, 0.6378, 0.9854, 0, 0.01978, 0), 
                                                                                    .Dim = 3:2, 
                                                                                    .Dimnames = list(c("Basal vs Other(s)", 
                                                                                                       "Her2 vs Other(s)", 
                                                                                                       "LumA vs Other(s)"), 
                                                                                                     c("AUC", "p-value"))), 
                                                                  comp2 = structure(c(0.9886, 0.9383, 0.9956, 0, 1.223e-13, 0), 
                                                                                    .Dim = 3:2, 
                                                                                    .Dimnames = list(c("Basal vs Other(s)", 
                                                                                                       "Her2 vs Other(s)", 
                                                                                                       "LumA vs Other(s)"), 
                                                                                                     c("AUC", "p-value"))))),
                      srbct.splsda = list(Comp1 = structure(c(0.3891, 1, 0.8088, 0.6547, 0.1453, 5.586e-06, 0.0009391, 0.04955), 
                                                            .Dim = c(4L, 2L), 
                                                            .Dimnames = list(c("EWS vs Other(s)", 
                                                                               "BL vs Other(s)", 
                                                                               "NB vs Other(s)", 
                                                                               "RMS vs Other(s)"), 
                                                                             c("AUC", "p-value"))), 
                                          Comp2 = structure(c(1, 1, 0.8791, 0.8081, 5.135e-11, 5.586e-06, 4.89e-05, 9.12e-05), 
                                                            .Dim = c(4L, 2L), 
                                                            .Dimnames = list(c("EWS vs Other(s)", 
                                                                               "BL vs Other(s)", 
                                                                               "NB vs Other(s)", 
                                                                               "RMS vs Other(s)"), 
                                                                             c("AUC", "p-value")))),
                      newdata.splsda = list(Comp1 = structure(c(0.6667, 0.5637), 
                                                              .Dim = 1:2, 
                                                              .Dimnames = list("AF vs BE", 
                                                                               c("AUC", "p-value"))), 
                                            Comp2 = structure(c(0.6667, 0.5637), 
                                                              .Dim = 1:2, 
                                                              .Dimnames =
                                                                  list("AF vs BE",
                                                                       c("AUC", "p-value")))),
                      multilevel.splsda = list(Comp1 = structure(c(0.5024, 0.9976, 1, 0.9603, 0, 0), 
                                                                 .Dim = 3:2, 
                                                                 .Dimnames = list(c("Antecubital_fossa vs Other(s)", 
                                                                                    "Stool vs Other(s)", 
                                                                                    "Subgingival_plaque vs Other(s)"), 
                                                                                  c("AUC", "p-value"))), 
                                               Comp2 = structure(c(1, 1, 1, 0, 0, 0), 
                                                                 .Dim = 3:2, 
                                                                 .Dimnames = list(c("Antecubital_fossa vs Other(s)", 
                                                                                    "Stool vs Other(s)", 
                                                                                    "Subgingival_plaque vs Other(s)"), 
                                                                                  c("AUC", "p-value")))),
                      roc.comp.splsda = list(Comp1 = structure(c(0.863, 2.473e-05), 
                                                               .Dim = 1:2, 
                                                               .Dimnames = list("AF vs BE", 
                                                                                c("AUC", "p-value"))), 
                                             Comp2 = structure(c(0.9981, 7.124e-09), 
                                                               .Dim = 1:2, 
                                                               .Dimnames = list("AF vs BE", 
                                                                                c("AUC", "p-value")))),
                      roc.study.mint.splsda = list(Comp1 = structure(c(0.863, 2.473e-05), 
                                                                     .Dim = 1:2, 
                                                                     .Dimnames = list("AF vs BE", 
                                                                                      c("AUC", "p-value")))),
                      roc.block.splsda = list(miRNA = list(comp1 = structure(c(0.9365, 0.5861, 0.9218, 0, 0.1453, 0), 
                                                                             .Dim = 3:2, 
                                                                             .Dimnames = list(c("Basal vs Other(s)",
                                                                                                "Her2 vs Other(s)", 
                                                                                                "LumA vs Other(s)"), 
                                                                                              c("AUC", "p-value"))), 
                                                           comp2 = structure(c(0.9494, 0.5856, 0.9257, 0, 0.1479, 0), 
                                                                             .Dim = 3:2, 
                                                                             .Dimnames = list(c("Basal vs Other(s)", 
                                                                                                "Her2 vs Other(s)", 
                                                                                                "LumA vs Other(s)"), 
                                                                                              c("AUC", "p-value")))), 
                                              mRNA = list(comp1 = structure(c(0.9937, 0.6206, 0.9918, 0, 0.04144, 0), 
                                                                            .Dim = 3:2, 
                                                                            .Dimnames = list(c("Basal vs Other(s)", 
                                                                                               "Her2 vs Other(s)", 
                                                                                               "LumA vs Other(s)"), 
                                                                                             c("AUC", "p-value"))), 
                                                          comp2 = structure(c(0.9962, 0.5725, 0.9893, 0, 0.2201, 0), 
                                                                            .Dim = 3:2, 
                                                                            .Dimnames = list(c("Basal vs Other(s)", 
                                                                                               "Her2 vs Other(s)", 
                                                                                               "LumA vs Other(s)"), 
                                                                                             c("AUC", "p-value")))), 
                                              proteomics = list(comp1 = structure(c(0.9564, 0.6578, 0.9844, 0, 0.007614, 0), 
                                                                                  .Dim = 3:2, 
                                                                                  .Dimnames = list(c("Basal vs Other(s)", 
                                                                                                     "Her2 vs Other(s)", 
                                                                                                     "LumA vs Other(s)"), 
                                                                                                   c("AUC", "p-value"))), 
                                                                comp2 = structure(c(0.9757, 0.6564, 0.9865, 0, 0.008164, 0), 
                                                                                  .Dim = 3:2, 
                                                                                  .Dimnames = list(c("Basal vs Other(s)", 
                                                                                                     "Her2 vs Other(s)", 
                                                                                                     "LumA vs Other(s)"), 
                                                                                                   c("AUC", "p-value"))))),
                      study.test.mint.splsda = list(Comp2 = structure(c(1, 0.8333, 1, 0.1573, 0.2482, 0.08326), 
                                                                      .Dim = 3:2, 
                                                                      .Dimnames = list(c("Fibroblast vs Other(s)", 
                                                                                         "hESC vs Other(s)",
                                                                                         "hiPSC vs Other(s)"), 
                                                                                       c("AUC", "p-value"))))
                      )

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

    roc.study.auroc = auroc(res.mint.plsda, roc.study = 2, print = FALSE)

    invisible(capture.output(TT <- dput(roc.comp.auroc[testable.components])))
    
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



