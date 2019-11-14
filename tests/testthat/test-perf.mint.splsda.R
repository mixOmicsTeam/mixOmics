context("test-tune.mint.splsda")


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
        signif.threshold = 0.1
    )
    
    out$choice.ncomp
    out$choice.keepX
    expect_is(out, "tune.mint.splsda")
    expect_equal(out$choice.ncomp$ncomp, 1)
    
})