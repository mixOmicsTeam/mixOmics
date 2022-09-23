###############################################################################
### ============================ MISSING TESTS ============================ ###
###############################################################################

# basic
## mint.(s)plsda

###############################################################################
### ================================ BASIC ================================ ###
###############################################################################


test_that("(background.predict:basic): (s)plsda", {
    
    data(breast.tumors)
    X <- breast.tumors$gene.exp
    Y <- breast.tumors$sample$treatment
    
    # --- PLS-DA --- # 
    res.plsda <- plsda(X, Y, ncomp = 2)
    
    res.bgp = background.predict(res.plsda, comp.predicted = 2, resolution = 10)
    
    expect_equal(c(sapply(res.bgp, class)),
                 rep(c("matrix", "array"),2))
    .expect_numerically_close(res.bgp$AF[c(1,length(res.bgp$AF))],
                              c(21.4039434, 0.3059848))
    .expect_numerically_close(res.bgp$BE[c(1,length(res.bgp$BE))],
                              c(-18.08271, -14.92983))
    
    # --- sPLS-DA --- # 
    choice.keepX <- c(10, 10)
    
    res.splsda <- splsda(X, Y, ncomp = 2, keepX = choice.keepX)
    
    res.bgp = background.predict(res.splsda, comp.predicted = 2, resolution = 10)
    
    expect_equal(c(sapply(res.bgp, class)),
                 rep(c("matrix", "array"),2))
    .expect_numerically_close(res.bgp$AF[c(1,length(res.bgp$AF))],
                              c(4.36332, -3.15242))
    .expect_numerically_close(res.bgp$BE[c(1,length(res.bgp$BE))],
                              c(-3.575758, -3.973891))
})


###############################################################################
### ================================ DATA ================================= ###
###############################################################################


test_that("(background.predict:data): liver.toxicity", {
    
    data(liver.toxicity)
    X = liver.toxicity$gene
    Y = as.factor(liver.toxicity$treatment[, 4])
    
    res.plsda <- plsda(X, Y, ncomp = 2)
    res.bgp = background.predict(res.plsda, comp.predicted = 2, resolution = 10)
    
    expect_equal(c(sapply(res.bgp, class)),
                 rep(c("matrix", "array"),4))
    .expect_numerically_close(res.bgp$`6`[c(1,length(res.bgp$`6`))],
                              c(22.53952, -20.37578))
    .expect_numerically_close(res.bgp$`48`[c(1,length(res.bgp$`48`))],
                              c(-63.46130, -28.74109))
})


test_that("(background.predict:data): srbct", {
    
    data(srbct)
    X <- srbct$gene
    Y <- srbct$class
    
    res.plsda <- plsda(X, Y, ncomp = 2)
    res.bgp = background.predict(res.plsda, comp.predicted = 2, resolution = 10)
    
    expect_equal(c(sapply(res.bgp, class)),
                 rep(c("matrix", "array"),4))
    .expect_numerically_close(res.bgp$EWS[c(1,length(res.bgp$EWS))],
                              c(-28.23491, -18.97067))
    .expect_numerically_close(res.bgp$RMS[c(1,length(res.bgp$RMS))],
                              c(-10.488557,6.814814))
})


###############################################################################
### ============================= PARAMETER =============================== ###
###############################################################################


test_that("(background.predict:parameter): comp.predicted", {
    
    data(breast.tumors)
    X <- breast.tumors$gene.exp
    Y <- breast.tumors$sample$treatment
    
    res.plsda <- plsda(X, Y, ncomp = 2)
    
    res.bgp = background.predict(res.plsda, comp.predicted = 1, resolution = 10)
    
    expect_equal(c(sapply(res.bgp, class)),
                 rep(c("matrix", "array"),2))
    .expect_numerically_close(res.bgp$AF[c(1,length(res.bgp$AF))],
                              c(5.60928, -14.92983))
    .expect_numerically_close(res.bgp$BE[c(1,length(res.bgp$BE))],
                              c(-18.08271, -14.92983))
})


test_that("(background.predict:parameter): dist", {
    
    data(breast.tumors)
    X <- breast.tumors$gene.exp
    Y <- breast.tumors$sample$treatment
    
    res.plsda <- plsda(X, Y, ncomp = 2)
    
    # --- max.dist --- #
    res.bgp = background.predict(res.plsda, comp.predicted = 2,
                                 dist = "max.dist", resolution = 10)
    
    expect_equal(c(sapply(res.bgp, class)),
                 rep(c("matrix", "array"),2))
    .expect_numerically_close(res.bgp$AF[c(1,length(res.bgp$AF))],
                              c(21.4039434, 0.3059848))
    
    # --- centroids.dist --- #
    res.bgp = background.predict(res.plsda, comp.predicted = 2,
                                 dist = "centroids.dist", resolution = 10)
    
    expect_equal(c(sapply(res.bgp, class)),
                 rep(c("matrix", "array"),2))
    .expect_numerically_close(res.bgp$AF[c(1,length(res.bgp$AF))],
                              c(13.506612, -8.835507))
    
    # --- mahalanobis.dist --- #
    res.bgp = background.predict(res.plsda, comp.predicted = 2,
                                 dist = "mahalanobis.dist", resolution = 10)
    
    expect_equal(c(sapply(res.bgp, class)),
                 rep(c("matrix", "array"),2))
    .expect_numerically_close(res.bgp$AF[c(1,length(res.bgp$AF))],
                              c(21.403943, -5.788343))
})


test_that("(background.predict:parameter): resolution", {
    
    data(breast.tumors)
    X <- breast.tumors$gene.exp
    Y <- breast.tumors$sample$treatment
    
    res.plsda <- plsda(X, Y, ncomp = 2)
    
    # --- medium resolution --- #
    res.bgp = background.predict(res.plsda, comp.predicted = 2,
                                 resolution = 30)
    
    expect_equal(c(sapply(res.bgp, class)),
                 rep(c("matrix", "array"),2))
    expect_equal(c(sapply(res.bgp, dim)),
                 c(102, 2, 88, 2))
    
    # --- high resolution --- #
    res.bgp = background.predict(res.plsda, comp.predicted = 2,
                                 resolution = 100)
    
    expect_equal(c(sapply(res.bgp, class)),
                 rep(c("matrix", "array"),2))
    expect_equal(c(sapply(res.bgp, dim)),
                 c(301, 2, 299, 2))
})


###############################################################################
### ================================ ERROR ================================ ###
###############################################################################



test_that("(background.predict:error): catches misuse of 'block.(s)plsda' objects", {
    
    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.train$mirna,
             mRNA = breast.TCGA$data.train$mrna,
             proteomics = breast.TCGA$data.train$protein)
    Y = breast.TCGA$data.train$subtype
    
    res.block.plsda <- block.plsda(X, Y, design = "full")
    
    expect_error(background.predict(res.block.plsda),
                 "'plsda'")
})


test_that("(background.predict:error): catches invalid 'dist' values", {
    
    data(breast.tumors)
    X <- breast.tumors$gene.exp
    Y <- breast.tumors$sample$treatment
    
    res.plsda <- plsda(X, Y, ncomp = 2)
    
    expect_error(background.predict(res.plsda, dist = "random.dist"),
                 "distances")
    
})


test_that("(background.predict:error): catches invalid 'comp.predicted' values", {
    
    data(breast.tumors)
    X <- breast.tumors$gene.exp
    Y <- breast.tumors$sample$treatment
    
    res.plsda <- plsda(X, Y, ncomp = 2)
    
    expect_error(background.predict(res.plsda, comp.predicted = 3),
                 "1 or 2 components")
})


test_that("(background.predict:error): catches invalid 'xlim' values", {
    
    data(breast.tumors)
    X <- breast.tumors$gene.exp
    Y <- breast.tumors$sample$treatment
    
    res.plsda <- plsda(X, Y, ncomp = 2)
    
    expect_error(background.predict(res.plsda, comp.predicted = 2, xlim = c(1)),
                 "xlim")
})


test_that("(background.predict:error): catches invalid 'ylim' values", {
    
    data(breast.tumors)
    X <- breast.tumors$gene.exp
    Y <- breast.tumors$sample$treatment
    
    res.plsda <- plsda(X, Y, ncomp = 2)
    
    expect_error(background.predict(res.plsda, comp.predicted = 2, ylim = c(1)),
                 "ylim")
})


test_that("(background.predict:error): catches invalid 'resolution' values", {
    
    data(breast.tumors)
    X <- breast.tumors$gene.exp
    Y <- breast.tumors$sample$treatment
    
    res.plsda <- plsda(X, Y, ncomp = 2)
    
    expect_error(background.predict(res.plsda, comp.predicted = 2, resolution = -1),
                 "resolution")
})