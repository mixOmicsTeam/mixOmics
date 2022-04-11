
###############################################################################
### ============================ MISSING TESTS ============================ ###
###############################################################################

# basic
## mint.plsda
## mint.splsda

###############################################################################
### ============================ GROUND TRUTHS ============================ ###
###############################################################################

source("ground.truths/background.predict.R")

###############################################################################
### ================================ BASIC ================================ ###
###############################################################################


test_that("(background.predict:basic): plsda", {
    
    testable.components <- Testable.Components$basic.plsda
    GT <- Ground.Truths$basic.plsda
    
    data(breast.tumors)
    X <- breast.tumors$gene.exp
    Y <- breast.tumors$sample$treatment
    
    res.plsda <- plsda(X, Y, ncomp = 2)
    
    plsda.bgp = background.predict(res.plsda, comp.predicted = 2)
    
    invisible(capture.output(TT <- dput(plsda.bgp[testable.components])))
    
    expect_equal(TT, GT)
})


# test_that("(background.predict:basic): splsda", {
#     
#     GH.URL <- "https://github.com/Max-Bladen/mixOmicsWork/blob/main/Test%20Ground%20Truths/background.predict/basic.splsda.RData?raw=true"
#     .load_url(GH.URL)
#     
#     data(breast.tumors)
#     X <- breast.tumors$gene.exp
#     Y <- breast.tumors$sample$treatment
#     
#     choice.keepX <- c(10, 10)
#     
#     res.splsda <- splsda(X, Y, ncomp = 2, keepX = choice.keepX)
#     
#     splsda.bgp = background.predict(res.splsda, comp.predicted = 2)
#     
#     expect_equal(splsda.bgp, GT.splsda.bgp)
# })
# 
# 
# ###############################################################################
# ### ================================ DATA ================================= ###
# ###############################################################################
# 
# 
# test_that("(background.predict:data): liver.toxicity", {
#     
#     GH.URL <- "https://github.com/Max-Bladen/mixOmicsWork/blob/main/Test%20Ground%20Truths/background.predict/liver.toxicity.plsda.RData?raw=true"
#     .load_url(GH.URL)
#     
#     data(liver.toxicity)
#     X = liver.toxicity$gene
#     Y = as.factor(liver.toxicity$treatment[, 4])
#     
#     res.plsda <- plsda(X, Y, ncomp = 2)
#     liver.toxicity.plsda.bgp = background.predict(res.plsda, comp.predicted = 2)
#     
#     expect_equal(liver.toxicity.plsda.bgp, GT.liver.toxicity.plsda.bgp)
# })
# 
# 
# test_that("(background.predict:data): srbct", {
#     
#     GH.URL <- "https://github.com/Max-Bladen/mixOmicsWork/blob/main/Test%20Ground%20Truths/background.predict/srbct.plsda.RData?raw=true"
#     .load_url(GH.URL)
#     
#     data(srbct)
#     X <- srbct$gene
#     Y <- srbct$class
#     
#     res.plsda <- plsda(X, Y, ncomp = 2)
#     srbct.plsda.bgp = background.predict(res.plsda, comp.predicted = 2)
#     
#     expect_equal(srbct.plsda.bgp, GT.srbct.plsda.bgp)
# })
# 
# 
# ###############################################################################
# ### ============================= PARAMETER =============================== ###
# ###############################################################################
# 
# 
# test_that("(background.predict:parameter): comp.predicted", {
#     
#     GH.URL <- "https://github.com/Max-Bladen/mixOmicsWork/blob/main/Test%20Ground%20Truths/background.predict/comp.predicted.plsda.RData?raw=true"
#     .load_url(GH.URL)
#     
#     data(breast.tumors)
#     X <- breast.tumors$gene.exp
#     Y <- breast.tumors$sample$treatment
#     
#     res.plsda <- plsda(X, Y, ncomp = 2)
#     
#     comp.predicted.bgp = background.predict(res.plsda, comp.predicted = 1)
#     
#     expect_equal(comp.predicted.bgp, GT.comp.predicted.bgp)
# })
# 
# 
# test_that("(background.predict:parameter): dist", {
#     
#     GH.URLs <- list("https://github.com/Max-Bladen/mixOmicsWork/blob/main/Test%20Ground%20Truths/background.predict/max.dist.plsda.RData?raw=true",
#                    "https://github.com/Max-Bladen/mixOmicsWork/blob/main/Test%20Ground%20Truths/background.predict/centroids.dist.plsda.RData?raw=true",
#                    "https://github.com/Max-Bladen/mixOmicsWork/blob/main/Test%20Ground%20Truths/background.predict/mahalanobis.dist.plsda.RData?raw=true")
#     invisible(lapply(GH.URLs, .load_url))
#     
#     data(breast.tumors)
#     X <- breast.tumors$gene.exp
#     Y <- breast.tumors$sample$treatment
#     
#     res.plsda <- plsda(X, Y, ncomp = 2)
#     
#     max.dist.bgp = background.predict(res.plsda, comp.predicted = 2,
#                                             dist = "max.dist")
#     
#     centroids.dist.bgp = background.predict(res.plsda, comp.predicted = 2,
#                                          dist = "centroids.dist")
#     
#     mahalanobis.dist.bgp = background.predict(res.plsda, comp.predicted = 2,
#                                          dist = "mahalanobis.dist")
#     
#     expect_equal(max.dist.bgp, GT.max.dist.bgp)
#     expect_equal(centroids.dist.bgp, GT.centroids.dist.bgp)
#     expect_equal(mahalanobis.dist.bgp, GT.mahalanobis.dist.bgp)
# })
# 
# 
# test_that("(background.predict:parameter): resolution", {
#     
#     GH.URLs <- list("https://github.com/Max-Bladen/mixOmicsWork/blob/main/Test%20Ground%20Truths/background.predict/resolution.20.plsda.RData?raw=true",
#                     "https://github.com/Max-Bladen/mixOmicsWork/blob/main/Test%20Ground%20Truths/background.predict/resolution.50.plsda.RData?raw=true")
#     invisible(lapply(GH.URLs, .load_url))
#     
#     data(breast.tumors)
#     X <- breast.tumors$gene.exp
#     Y <- breast.tumors$sample$treatment
#     
#     res.plsda <- plsda(X, Y, ncomp = 2)
#     
#     res.20.bgp = background.predict(res.plsda, comp.predicted = 2,
#                                        resolution = 20)
#     
#     res.50.bgp = background.predict(res.plsda, comp.predicted = 2,
#                                        resolution = 50)
#     
#     expect_equal(res.20.bgp, GT.res.20.bgp)
#     expect_equal(res.50.bgp, GT.res.50.bgp)
# })
# 
# 
# ###############################################################################
# ### ================================ ERROR ================================ ###
# ###############################################################################
# 
# 
# 
# test_that("(background.predict:error): cannot use block.(s)plsda objects", {
#     
#     data(breast.TCGA)
#     X = list(miRNA = breast.TCGA$data.train$mirna,
#              mRNA = breast.TCGA$data.train$mrna,
#              proteomics = breast.TCGA$data.train$protein)
#     Y = breast.TCGA$data.train$subtype
#     
#     res.block.plsda <- block.plsda(X, Y, design = "full")
#     
#     expect_error(background.predict(res.block.plsda),
#                  "'background.predict' can only be calculated for 'plsda'
#         and 'splsda' objects",
#                  fixed = TRUE)
#     
#     choice.keepX <- list(miRNA=c(10,10),
#                          mRNA=c(10,10), 
#                          proteomics=c(10,10))
#     
#     res.block.splsda <- block.splsda(X, Y, design = "full", keepX = choice.keepX)
#     
#     expect_error(background.predict(res.block.splsda),
#                  "'background.predict' can only be calculated for 'plsda'
#         and 'splsda' objects",
#                  fixed = TRUE)
# })
# 
# 
# test_that("(background.predict:error): ensure dist has valid value", {
#     
#     data(breast.tumors)
#     X <- breast.tumors$gene.exp
#     Y <- breast.tumors$sample$treatment
#     
#     res.plsda <- plsda(X, Y, ncomp = 2)
#     
#     expect_error(background.predict(res.plsda, dist = "incorrect.dist"),
#                  "Choose one of the three following distances: 'max.dist',
#         'centroids.dist' or 'mahalanobis.dist'",
#                  fixed = TRUE)
# 
# })
# 
# 
# test_that("(background.predict:error): ensure dist has valid value", {
#     
#     data(breast.tumors)
#     X <- breast.tumors$gene.exp
#     Y <- breast.tumors$sample$treatment
#     
#     res.plsda <- plsda(X, Y, ncomp = 2)
#     
#     expect_error(background.predict(res.plsda, comp.predicted = 3),
#                  "Can only show predicted background for 1 or 2 components",
#                  fixed = TRUE)
# })


