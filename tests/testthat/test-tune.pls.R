context("tune.pls")
library(BiocParallel)

test_that("tune.pls works and is the same as perf alone and in tune wrapper", code = {
  
  # set up data
  data(liver.toxicity)
  X <- liver.toxicity$gene
  Y <- liver.toxicity$clinic
  
  # independently
  tune.res.1 <- tune.pls( X, Y, ncomp = 2, measure = "cor", logratio = "none",
                        folds = 5, nrepeat = 1, progressBar = FALSE, seed = 100)
  
  # in tune wrapper
  tune.res.2 <- tune( X, Y, ncomp = 2, measure = "cor", logratio = "none",
                          folds = 5, nrepeat = 1, method = "pls", seed = 100)
  
  # with perf
  model <- pls(X, Y, ncomp = 10, logratio = "none")
  tune.res.3 <- perf(model, folds = 5, nrepeat = 1, seed = 100)
  
  # check outputs format
  expect_equal(class(tune.res.1)[1], "perf.pls.mthd")
  expect_equal(class(tune.res.2)[1], "perf.pls.mthd")
  expect_equal(class(tune.res.3)[1], "perf.pls.mthd")
  
  # check output values
  .expect_numerically_close(tune.res.1$measures$Q2$summary[1,3], -0.1304012)
  .expect_numerically_close(tune.res.2$measures$Q2$summary[1,3], -0.1304012)
  .expect_numerically_close(tune.res.3$measures$Q2$summary[1,3], -0.1304012)
  
  # check can plot outputs
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot(tune.res.1))
  expect_silent(plot(tune.res.2))
  expect_silent(plot(tune.res.3))
  
})
  