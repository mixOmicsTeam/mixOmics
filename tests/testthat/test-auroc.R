# basic
## mint.block.plsda
## mint.block.splsda

# data


# parameter
##  newdata
##  multilevel
##  plot ??? can we test this?
##  roc.comp --|
##  roc.study -|-- maybe one test?
##  roc.block -|
##  study.test

# error
## auroc.mint.block.plsda
#### roc.block cannot be greater than %s AND roc.block should be integer or character 

## auroc.mint.plsda
#### factor outcome.test must be a factor with ... elements
#### roc.study' must be a single entry, either `global' or one of levels(object$study)

## auroc.mixo_plsda
#### Factor outcome.test must be a factor with ... elements.

## predict.mixo_spls
#### no prediction for RGCCA methods
#### ERROR : choose one of the four following modes: 'all', 'max.dist', 'centroids.dist' or 'mahalanobis.dist'
#### 'newdata' must include all the variables of 'object$X
#### samples should have a unique identifier/rowname
#### Some blocks are missing in 'newdata'; the prediction is based on the following blocks only



test_that("(basic): plsda auroc", {
    
    GH.URL <- "https://github.com/Max-Bladen/mixOmicsWork/blob/main/Test%20Ground%20Truths/auroc/basic.plsda.RData?raw=true"
    .load_url(GH.URL)
    
    data(breast.tumors)
    X <- breast.tumors$gene.exp
    Y <- breast.tumors$sample$treatment
    
    res.plsda <- plsda(X, Y, ncomp = 2)
    
    plsda.auroc = auroc(res.plsda, plot = TRUE, roc.comp = 1, print = FALSE)
    
    expect_equal(plsda.auroc, GT.plsda.auroc)
 })


test_that("(basic): splsda auroc", {
    
    GH.URL <- "https://github.com/Max-Bladen/mixOmicsWork/blob/main/Test%20Ground%20Truths/auroc/basic.splsda.RData?raw=true"
    .load_url(GH.URL)
    
    data(breast.tumors)
    X <- breast.tumors$gene.exp
    Y <- breast.tumors$sample$treatment
    
    choice.keepX <- c(10, 10)
    
    res.splsda <- splsda(X, Y, ncomp = 2, keepX = choice.keepX)
    
    splsda.auroc = auroc(res.splsda, plot = TRUE, roc.comp = 1, print = FALSE)
    
    expect_equal(splsda.auroc, GT.splsda.auroc)
})


test_that("(basic): mint.plsda auroc", {
    
    GH.URL <- "https://github.com/Max-Bladen/mixOmicsWork/blob/main/Test%20Ground%20Truths/auroc/basic.mint.plsda.RData?raw=true"
    .load_url(GH.URL)
   
    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    s <- stemcells$study
    
    res.mint.plsda <- mint.plsda(X, Y, ncomp = 2, study = s)
    
    mint.plsda.auroc = auroc(res.mint.plsda, plot = TRUE, roc.comp = 1, print = FALSE)
    
    expect_equal(mint.plsda.auroc, GT.mint.plsda.auroc)
})


test_that("(basic): mint.splsda auroc", {
    
    GH.URL <- "https://github.com/Max-Bladen/mixOmicsWork/blob/main/Test%20Ground%20Truths/auroc/basic.mint.splsda.RData?raw=true"
    .load_url(GH.URL)
    
    data(stemcells)
    X <- stemcells$gene
    Y <- stemcells$celltype
    s <- stemcells$study
    
    choice.keepX <- c(10,10)
    
    res.mint.splsda <- mint.splsda(X, Y, ncomp = 2, study = s, keepX = choice.keepX)
    
    mint.splsda.auroc = auroc(res.mint.splsda, plot = TRUE, roc.comp = 1, print = FALSE)
    
    expect_equal(mint.splsda.auroc, GT.mint.splsda.auroc)
})


test_that("(basic): block.plsda auroc", {
    
    GH.URL <- "https://github.com/Max-Bladen/mixOmicsWork/blob/main/Test%20Ground%20Truths/auroc/basic.block.plsda.RData?raw=true"
    .load_url(GH.URL)
    
    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.train$mirna,
             mRNA = breast.TCGA$data.train$mrna,
             proteomics = breast.TCGA$data.train$protein)
    Y = breast.TCGA$data.train$subtype
    
    res.block.plsda <- block.plsda(X, Y, design = "full")
    
    block.plsda.auroc = auroc(res.block.plsda, plot = TRUE, roc.comp = 1, print = FALSE)
    
    expect_equal(block.plsda.auroc, GT.block.plsda.auroc)
})


test_that("(basic): block.splsda auroc", {
    
    GH.URL <- "https://github.com/Max-Bladen/mixOmicsWork/blob/main/Test%20Ground%20Truths/auroc/basic.block.splsda.RData?raw=true"
    .load_url(GH.URL)
    
    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.train$mirna,
             mRNA = breast.TCGA$data.train$mrna,
             proteomics = breast.TCGA$data.train$protein)
    Y = breast.TCGA$data.train$subtype
    
    choice.keepX <- list(miRNA=c(10,10),
                         mRNA=c(10,10), 
                         proteomics=c(10,10))
    
    res.block.splsda <- block.splsda(X, Y, design = "full", keepX = choice.keepX)
    
    block.splsda.auroc = auroc(res.block.splsda, plot = TRUE, roc.comp = 1, print = FALSE)
    
    expect_equal(block.splsda.auroc, GT.block.splsda.auroc)
})


test_that("(data): splsda auroc, srbct", {
    
    GH.URL <- "https://github.com/Max-Bladen/mixOmicsWork/blob/main/Test%20Ground%20Truths/auroc/srbct.plsda.RData?raw=true"
    .load_url(GH.URL)
    
    data(srbct)
    X <- srbct$gene
    Y <- srbct$class
    
    choice.keepX <- c(10, 10)
    
    res.splsda <- splsda(X, Y, ncomp = 2, keepX = choice.keepX)
    
    srbct.splsda.auroc = auroc(res.splsda, plot = TRUE, roc.comp = 1, print = FALSE)
    
    expect_equal(srbct.splsda.auroc, GT.srbct.splsda.auroc)
})






# save(GT.plsda.auroc, file = "C:/Users/Work/Desktop/mO Work/Test Ground Truths/auroc/basic.plsda.RData")