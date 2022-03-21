context("perf.mixo_splsda")

test_that("perf.mixo_splsda functions", code = {
    data(liver.toxicity)
    X <- liver.toxicity$gene
    Y <- liver.toxicity$treatment$Dose.Group
    # create a class with one sample only
    #Y[c(1,2)] <- 'FOO'
    
    res <- plsda(X, Y, ncomp = 2)
    set.seed(12)
    out <- perf(res, validation = "Mfold", folds = 3, nrepeat = 3)
    expect_is(out, "perf")
    
    ground.ncomp <- matrix(c(2,1,2,2,1,2), ncol = 3, byrow=T,
                           dimnames = list(c("overall", "BER"),
                                           c("max.dist", "centroids.dist", "mahalanobis.dist")))
    
    expect_true(all(out$choice.ncomp == ground.ncomp))
})

test_that("perf.mixo_splsda does not allow for class with 1 associated sample", code = {
    data(liver.toxicity)
    X <- liver.toxicity$gene
    Y <- liver.toxicity$treatment$Dose.Group
    # create a class with one sample only
    Y[c(1)] <- 'asdf'
    
    res <- plsda(X, Y, ncomp = 2)
    set.seed(12)
    expect_error(perf(res, validation = "Mfold", folds = 3, nrepeat = 3),
                 paste("Cannot evaluate performance when a class level ('", 
                       names(table(res$Y))[which(table(res$Y) == 1)], 
                       "') has only a single assocaited sample.", 
                       sep = ""),
                 fixed = T)
})
