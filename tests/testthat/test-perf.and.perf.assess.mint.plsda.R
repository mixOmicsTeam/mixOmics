context("perf.assess.mint.plsda")
library(BiocParallel)

## ------------------------------------------------------------------------ ##
## Test perf.mint.splsda()

test_that("perf.mint.splsda works", code = {
    data(stemcells)
    res = mint.splsda(
        X = stemcells$gene,
        Y = stemcells$celltype,
        ncomp = 3,
        keepX = c(5, 10, 15),
        study = stemcells$study
    )
    
    out.perf = perf(res, auc = FALSE)
    expect_is(out.perf, "perf")
    expect_true(all(out.perf$choice.ncomp == 1))
    
    out.perf.assess <- perf.assess(res, auc = FALSE)
    expect_is(out.perf.assess, "perf")
    expect_equal(out.perf.assess$choice.ncomp, NULL)
    
    expect_equal(out.perf$study.specific.error$`1`$BER[3,], out.perf.assess$study.specific.error$`1`$BER[1,])
    
})

test_that("perf.mint.splsda works with custom alpha", code = {
    data(stemcells)
    res = mint.splsda(
        X = stemcells$gene,
        Y = stemcells$celltype,
        ncomp = 3,
        keepX = c(5, 10, 15),
        study = stemcells$study
    )
    
    out = perf(res, auc = FALSE, alpha = 0.1)
    expect_is(out, "perf")
    expect_true(all(out$choice.ncomp == 1))
    
})

## ------------------------------------------------------------------------ ##
## Test perf.assess.mixo_plsda() give informative error message when one sample in one class
# can't run this because get error message already when try to build model like this

# test_that("perf.assess.mixo_plsda error when one sample in one class", code = {
#   
#   data(stemcells)
#   res = mint.splsda(
#     X = stemcells$gene[1:33,],
#     Y = stemcells$celltype[1:33],
#     ncomp = 3,
#     keepX = c(5, 10, 15),
#     study = stemcells$study[1:33]
#   )
#   
#   Error in internal_wrapper.mint(X = X, Y = Y.mat, ncomp = ncomp, near.zero.var = near.zero.var,  : 
#   At least one study has only one sample, please consider removing before calling the function again
# })

# test_that("perf.mint.splsda and tune.mint.splsda yield the same error rates for all measures", code = {
#     data(stemcells)
#     X = stemcells$gene
#     Y = stemcells$celltype
#     study <- stemcells$study
#     
#     metricsC1 <- matrix(0, nrow = 2, ncol = 3)
#     colnames(metricsC1) <- c('max.dist', 'centroids.dist', 'mahalanobis.dist')
#     rownames(metricsC1) <- c("overall", "BER")
#     
#     metricsC2 <- metricsC1
#     
#     for (dist in c('max.dist', 'centroids.dist', 'mahalanobis.dist')) {
#         for (measure in c("overall", "BER")) {
#             tune.mint = tune.mint.splsda(X = X, Y = Y, study = study, ncomp = 2, test.keepX = seq(1, 51, 5),
#                                          dist = dist, progressBar = FALSE, measure = measure)
#             
#             mint.splsda.res = mint.splsda(X = X, Y = Y, study = study, ncomp = 2,
#                                           keepX = tune.mint$choice.keepX)
#             
#             perf.mint = perf(mint.splsda.res, progressBar = FALSE, dist = dist)
#             
#             metricsC1[measure, dist] <- tune.mint$error.rate[which(rownames(tune.mint$error.rate) == tune.mint$choice.keepX[1]), "comp1"]
#             metricsC2[measure, dist] <- tune.mint$error.rate[which(rownames(tune.mint$error.rate) == tune.mint$choice.keepX[2]), "comp2"]
#             
#             expect_equal(round(perf.mint$global.error[[measure]]["comp1",], 4), round(metricsC1[[measure, dist]], 4))
#             expect_equal(round(perf.mint$global.error[[measure]]["comp2",], 4), round(metricsC2[[measure, dist]], 4))
#         }
#     }
# })
