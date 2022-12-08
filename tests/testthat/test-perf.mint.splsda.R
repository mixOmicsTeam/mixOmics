context("tune.mint.splsda")


test_that("perf.mint.splsda works", code = {
    data(stemcells)
    res = mint.splsda(
        X = stemcells$gene,
        Y = stemcells$celltype,
        ncomp = 3,
        keepX = c(5, 10, 15),
        study = stemcells$study
    )
    
    out = perf(res, auc = FALSE)
    expect_is(out, "perf")
    expect_true(all(out$choice.ncomp == 1))
    
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

test_that("perf.mint.splsda and tune.mint.splsda yield the same error rates for all measures", code = {
    data(stemcells)
    X = stemcells$gene
    Y = stemcells$celltype
    study <- stemcells$study
    
    metricsC1 <- matrix(0, nrow = 2, ncol = 3)
    colnames(metricsC1) <- c('max.dist', 'centroids.dist', 'mahalanobis.dist')
    rownames(metricsC1) <- c("overall", "BER")
    
    metricsC2 <- metricsC1
    
    for (dist in c('max.dist', 'centroids.dist', 'mahalanobis.dist')) {
        for (measure in c("overall", "BER")) {
            tune.mint = tune.mint.splsda(X = X, Y = Y, study = study, ncomp = 2, test.keepX = seq(1, 51, 5),
                                         dist = dist, progressBar = FALSE, measure = measure)
            
            mint.splsda.res = mint.splsda(X = X, Y = Y, study = study, ncomp = 2,
                                          keepX = tune.mint$choice.keepX)
            
            perf.mint = perf(mint.splsda.res, progressBar = FALSE, dist = dist)
            
            metricsC1[measure, dist] <- tune.mint$error.rate[which(rownames(tune.mint$error.rate) == tune.mint$choice.keepX[1]), "comp1"]
            metricsC2[measure, dist] <- tune.mint$error.rate[which(rownames(tune.mint$error.rate) == tune.mint$choice.keepX[2]), "comp2"]
            
            expect_equal(round(perf.mint$global.error[[measure]]["comp1",], 4), round(metricsC1[[measure, dist]], 4))
            expect_equal(round(perf.mint$global.error[[measure]]["comp2",], 4), round(metricsC2[[measure, dist]], 4))
        }
    }
})
