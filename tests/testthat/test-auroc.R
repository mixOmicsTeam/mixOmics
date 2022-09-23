
###############################################################################
### ============================ MISSING TESTS ============================ ###
###############################################################################

# basic
## mint.block.plsda
## mint.block.splsda - created a node within the Project workflow. Once addressed, add tests

###############################################################################
### ================================ BASIC ================================ ###
###############################################################################


test_that("(auroc:basic): (s)plsda", {
    
    data(breast.tumors)
    X <- breast.tumors$gene.exp
    Y <- breast.tumors$sample$treatment
    
    # --- PLS-DA -- # 
    res.plsda <- plsda(X, Y, ncomp = 2)
    
    res.auroc = .quiet(auroc(res.plsda, print = TRUE))
    
    expect_equal(typeof(res.auroc), "list")
    expect_equal(matrix(res.auroc$Comp1), 
                 rbind(8.630e-01, 2.473e-05))
    expect_equal(matrix(res.auroc$Comp2), 
                 rbind(9.981e-01, 7.124e-09))
    
    # --- sPLS-DA -- # 
    choice.keepX <- c(10, 10)
    
    res.splsda <- splsda(X, Y, ncomp = 2, keepX = choice.keepX)
    
    res.auroc = auroc(res.splsda, print = FALSE)
    
    expect_equal(typeof(res.auroc), "list")
    expect_equal(matrix(res.auroc$Comp1), 
                 rbind(9.704e-01, 4.623e-08))
    expect_equal(matrix(res.auroc$Comp2), 
                 rbind(9.981e-01, 7.124e-09))
})


test_that("(auroc:basic): mint.(s)plsda", {
    
    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    s <- stemcells$study
    
    # --- MINT.PLS-DA -- # 
    res.mint.plsda <- mint.plsda(X, Y, ncomp = 2, study = s)
    
    res.auroc = .quiet(auroc(res.mint.plsda, roc.comp = 1, print = TRUE))
    
    expect_equal(typeof(res.auroc), "list")
    expect_equal(matrix(res.auroc$Comp1)[c(1,6),1], 
                 c(9.930e-01, 3.671e-10))
    
    # --- MINT.sPLS-DA -- # 
    choice.keepX <- c(10,10)
    
    res.mint.splsda <- mint.splsda(X, Y, ncomp = 2, study = s, keepX = choice.keepX)
    
    res.auroc = auroc(res.mint.splsda, roc.comp = 1, print = FALSE)
    
    expect_equal(typeof(res.auroc), "list")
    expect_equal(matrix(res.auroc$Comp1)[c(1,6),1], 
                 c(1.0000, 3.300e-09))
})


test_that("(auroc:basic): block.(s)plsda", {
    
    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.train$mirna,
             mRNA = breast.TCGA$data.train$mrna,
             proteomics = breast.TCGA$data.train$protein)
    Y = breast.TCGA$data.train$subtype
    
    # --- BLOCK.PLS-DA -- # 
    res.block.plsda <- block.plsda(X, Y, design = "full")
    
    res.auroc = auroc(res.block.plsda, roc.comp = 1, print = FALSE)
    
    expect_equal(typeof(res.auroc), "list")
    expect_equal(matrix(res.auroc$miRNA$comp1)[c(1,6),1], 
                 c(0.9365, 0.0000))
    expect_equal(matrix(res.auroc$proteomics$comp2)[c(1,6),1], 
                 c(0.9757, 0.000000))
    
    # --- BLOCK.sPLS-DA -- # 
    choice.keepX <- list(miRNA=c(10,10),
                         mRNA=c(10,10),
                         proteomics=c(10,10))
    
    res.block.splsda <- block.splsda(X, Y, design = "full", keepX = choice.keepX)
    
    res.auroc = .quiet(auroc(res.block.splsda, roc.comp = 1, print = TRUE))
    
    expect_equal(typeof(res.auroc), "list")
    expect_equal(matrix(res.auroc$miRNA$comp1)[c(1,6),1], 
                 c(0.9549, 0.0000))
    expect_equal(matrix(res.auroc$proteomics$comp2)[c(1,6),1], 
                 c(0.9886, 0.0000))
})


###############################################################################
### ================================ DATA ================================= ###
###############################################################################


test_that("(auroc:data): splsda, srbct", {
    
    data(srbct)
    X <- srbct$gene
    Y <- srbct$class
    
    choice.keepX <- c(10, 10)
    
    res.splsda <- splsda(X, Y, ncomp = 2, keepX = choice.keepX)
    
    res.auroc = auroc(res.splsda, roc.comp = 1, print = FALSE)
    
    expect_equal(typeof(res.auroc), "list")
    expect_equal(matrix(res.auroc$Comp1)[c(1,8),1], 
                 c(0.38910, 0.04955))
    expect_equal(matrix(res.auroc$Comp2)[c(1,8),1], 
                 c(1, 9.12e-05))
})


###############################################################################
### ============================= PARAMETER =============================== ###
###############################################################################


test_that("(auroc:parameter): newdata/outcome.test", {
    
    data(breast.tumors)
    X <- breast.tumors$gene.exp
    Y <- breast.tumors$sample$treatment
    
    d <- .minimal_train_test_subset(X, as.factor(Y), n.tr = 15, n.te = 3)
    
    res.plsda <- plsda(d$X.tr, d$Y.tr, ncomp = 2)
    
    res.auroc = auroc(res.plsda, print = FALSE, roc.comp = 1,
                      newdata = d$X.te, outcome.test = d$Y.te)
    
    expect_equal(typeof(res.auroc), "list")
    expect_equal(matrix(res.auroc$Comp1), 
                 rbind(0.7778, 0.2752))
    expect_equal(matrix(res.auroc$Comp2),
                 rbind(0.5556, 0.8273))
})


test_that("(auroc:parameter): multilevel", {
    
    data(diverse.16S)
    X <- diverse.16S$data.TSS
    Y <- diverse.16S$bodysite
    mL <- diverse.16S$sample
    
    res.plsda <- plsda(X, Y, ncomp = 2)
    
    res.auroc <- auroc(res.plsda, print = FALSE, roc.comp = 2,
                       multilevel = mL)
    
    expect_equal(typeof(res.auroc), "list")
    expect_equal(matrix(res.auroc$Comp1)[c(1,6),1], 
                 c(0.5024, 0.0000))
    expect_equal(matrix(res.auroc$Comp2)[c(1,6),1],
                 c(1, 0))
})


test_that("(auroc:parameter): roc.comp", {
    
    data(breast.tumors)
    X <- breast.tumors$gene.exp
    Y <- breast.tumors$sample$treatment
    
    res.plsda <- plsda(X, Y, ncomp = 2)
    
    res.auroc = auroc(res.plsda, roc.comp = 2, print = FALSE)
    
    expect_equal(typeof(res.auroc), "list")
    expect_equal(matrix(res.auroc$Comp1), 
                 rbind(0.863, 2.473e-05))
    expect_equal(matrix(res.auroc$Comp2),
                 rbind(0.9981, 7.124e-09))
})


test_that("(auroc:parameter): roc.study", {
    
    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    s <- stemcells$study
    
    res.mint.plsda <- mint.plsda(X, Y, ncomp = 2, study = s)
    
    res.auroc = auroc(res.mint.plsda, roc.study = 2, print = FALSE)
    
    expect_equal(typeof(res.auroc), "list")
    expect_equal(matrix(res.auroc$Comp2)[c(1,6),1],
                 c(1, 3.751e-09))
})


test_that("(auroc:parameter): roc.block", {
    
    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.train$mirna,
             mRNA = breast.TCGA$data.train$mrna,
             proteomics = breast.TCGA$data.train$protein)
    Y = breast.TCGA$data.train$subtype
    
    res.block.plsda <- block.plsda(X, Y, ncomp = 2)
    
    # --- as integer --- #
    res.auroc <- auroc(res.block.plsda, roc.block = 2, print = FALSE)
    
    expect_equal(typeof(res.auroc), "list")
    expect_equal(matrix(res.auroc$miRNA$comp1)[c(1,6),1], 
                 c(0.9365, 0.0000))
    expect_equal(matrix(res.auroc$proteomics$comp2)[c(1,6),1], 
                 c(0.9757, 0.000000))
    
    # --- as character --- #
    res.auroc <- auroc(res.block.plsda, roc.block = "mRNA", print = FALSE)
    
    expect_equal(typeof(res.auroc), "list")
    expect_equal(matrix(res.auroc$miRNA$comp1)[c(1,6),1], 
                 c(0.9365, 0.0000))
    expect_equal(matrix(res.auroc$proteomics$comp2)[c(1,6),1], 
                 c(0.9757, 0.000000))
})


test_that("(auroc:parameter): study.test", {
    
    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    S <- stemcells$study
    
    d <- .minimal_train_test_subset(X, Y, S=S)
    
    res.mint.plsda <- mint.plsda(d$X.tr, d$Y.tr, ncomp = 2, study = d$S.tr)
    
    res.auroc = auroc(res.mint.plsda, print = FALSE,
                      newdata = d$X.te, outcome.test = d$Y.te,
                      study.test = d$S.te)
    
    expect_equal(typeof(res.auroc), "list")
    expect_equal(matrix(res.auroc$Comp2)[c(1,6),1], 
                 c(1, 0.1797))
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
