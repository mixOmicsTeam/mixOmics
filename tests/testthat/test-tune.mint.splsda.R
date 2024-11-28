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
    
    # check can plot
    pdf(NULL)
    on.exit(dev.off())
    expect_error(plot(out), NA) # makes note about not being able to plot SD bars
    
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
    
    # check can plot
    pdf(NULL)
    on.exit(dev.off())
    expect_error(plot(out), NA) # makes note about not being able to plot SD bars
    
})

test_that("tune.mint.splsda works when test.keepX = NULL and gives same result as perf()", code = {
  
  # set up data
  data(stemcells)
  X = stemcells$gene
  Y = stemcells$celltype
  study = stemcells$study
  
  # tune on components only
  tune_res <- suppressWarnings(
    tune.mint.splsda(X, Y, study = study, ncomp = 2,
                test.keepX = NULL)
  )
  
  # run perf
  mint.splsda_res <- mint.splsda(X, Y, study = study, ncomp = 2)
  perf_res <- suppressWarnings(
    perf(mint.splsda_res, ncomp = 2, dist = "max.dist")
  )
  
  # check results match
  expect_equal(tune_res$global.error$BER[1,1], perf_res$global.error$BER[1,1])
})
