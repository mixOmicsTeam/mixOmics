test_that("perf.mixo_splsda functions", code = {
    data(liver.toxicity)
    X <- liver.toxicity$gene
    Y <- liver.toxicity$treatment$Dose.Group

    set.seed(12)
    res <- plsda(X, Y, ncomp = 2)
    out <- perf(res, validation = "Mfold", folds = 3, nrepeat = 3)

    ground.ncomp <- matrix(c(2,1,2,2,1,2), ncol = 3, byrow=T,
                           dimnames = list(c("overall", "BER"),
                                           c("max.dist", "centroids.dist", "mahalanobis.dist")))

    expect_equal(out$choice.ncomp, ground.ncomp)
})

test_that("does not allow for class with 1 associated sample", code = {
    data(liver.toxicity)
    X <- liver.toxicity$gene
    Y <- liver.toxicity$treatment$Dose.Group
    # create a class with one sample only
    Y[c(1)] <- 'random.class'

    res <- plsda(X, Y, ncomp = 2)

    expect_error(perf(res, validation = "Mfold", folds = 3, nrepeat = 3),
                 "single assocaited")
})