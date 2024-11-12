context("perf.assess.mixo_splsda")
library(BiocParallel)

## ------------------------------------------------------------------------ ##
## Test perf.assess.mixo_plsda()

test_that("perf.assess.mixo_plsda in serial and parallel", code = {
  data(liver.toxicity)
  X <- liver.toxicity$gene
  Y <- liver.toxicity$treatment$Dose.Group
  res <- plsda(X, Y, ncomp = 2)
  
  # Serial execution
  out <- suppressWarnings(
    perf.assess(res, validation = "Mfold", folds = 2, nrepeat = 1, BPPARAM = SerialParam(), seed = 12)
  )
  
  # Parallel execution
  out.parallel <- suppressWarnings(
    perf.assess(res, validation = "Mfold", folds = 2, nrepeat = 1, BPPARAM = SnowParam(workers = 2), seed = 12)
  )
  
  # Check result
  expect_equal(out$choice.ncomp, NULL)
  
  # Expect the same result
  expect_equal(out$error.rate$overall, out.parallel$error.rate$overall)
})