context("tune.mint.splsda")

test_that("tune.mint.splsda works", code = {
    data(stemcells)
    data = stemcells$gene
    type.id = stemcells$celltype
    exp = stemcells$study
    
    out = tune.mint.splsda(
        X = data,
        Y = type.id,
        ncomp = 2,
        near.zero.var = FALSE,
        study = exp,
        test.keepX = seq(1, 5, 2)
    )
    
    out$choice.ncomp
    out$choice.keepX
    expect_is(out, "tune.mint.splsda")
    expect_equal(out$choice.ncomp$ncomp, 1)
    
})

test_that("tune.mint.splsda works with custom alpha", code = {
    data(stemcells)
    data = stemcells$gene
    type.id = stemcells$celltype
    exp = stemcells$study
    
    out = tune.mint.splsda(
        X = data,
        Y = type.id,
        ncomp = 2,
        near.zero.var = FALSE,
        study = exp,
        test.keepX = seq(1, 5, 2),
        signif.threshold = 0.05
    )
    
    out$choice.ncomp
    out$choice.keepX
    expect_is(out, "tune.mint.splsda")
    expect_equal(out$choice.ncomp$ncomp, 1)
    
})
