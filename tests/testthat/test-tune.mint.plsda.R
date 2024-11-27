context("tune.mint.plsda")
library(BiocParallel)

test_that("tune.mint.plsda works and is the same perf alone and in tune wrapper", code = {
  
  # set up data
  data(stemcells)
  X = stemcells$gene
  Y = stemcells$celltype
  study = stemcells$study
  
  # run alone
  tune.res.1 = suppressWarnings(
    tune.mint.plsda(X, Y, study = study, ncomp = 2)
  )
  
  # run in wrapper
  tune.res.2 = suppressWarnings(
    tune(method = "mint.plsda", X, Y, study = study, ncomp = 2)
  )
  
  # run perf
  model <- mint.plsda(X, Y, ncomp = 2, study = study)
  tune.res.3 = suppressWarnings(
    perf(model, ncomp = 2)
  )
  
  # check outputs format
  expect_equal(class(tune.res.1)[1], "perf")
  expect_equal(class(tune.res.2)[1], "perf")
  expect_equal(class(tune.res.3)[1], "perf")
  
  # check outputs values
  .expect_numerically_close(tune.res.1$global.error$BER[1,1], 0.3803556)
  .expect_numerically_close(tune.res.2$global.error$BER[1,1], 0.3803556)
  .expect_numerically_close(tune.res.3$global.error$BER[1,1], 0.3803556)
  
  # check can plot outputs
  expect_silent(plot(tune.res.1))
  expect_silent(plot(tune.res.2))
  expect_silent(plot(tune.res.3))
  
  
})