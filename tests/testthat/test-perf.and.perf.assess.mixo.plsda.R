context("perf.mixo.plsda")
library(BiocParallel)

## ------------------------------------------------------------------------ ##
## Test perf.mixo_plsda() and perf.assess.mixo_plsda()

test_that("perf.assess.mixo_plsda works same as perf for same components in serial and in parallel", code = {
  
  # set up data and model
  data(liver.toxicity)
  X <- liver.toxicity$gene
  Y <- liver.toxicity$treatment$Dose.Group
  res <- plsda(X, Y, ncomp = 2)
  
  # Execution using old perf() function in serial
  out.1 <- suppressWarnings(
    perf(res, validation = "Mfold", folds = 2, nrepeat = 1, BPPARAM = SerialParam(), seed = 12)
  )
  
  # Execution using old perf() function in parallel
  out.2 <- suppressWarnings(
    perf(res, validation = "Mfold", folds = 2, nrepeat = 1, BPPARAM = SnowParam(workers = 2), seed = 12)
  )
  
  # New perf.assess in series
  out.3 <- suppressWarnings(
    perf.assess(res, validation = "Mfold", folds = 2, nrepeat = 1, BPPARAM = SerialParam(), seed = 12)
  )
  
  # New perf.assess in parallel
  out.4 <- suppressWarnings(
    perf.assess(res, validation = "Mfold", folds = 2, nrepeat = 1, BPPARAM = SnowParam(workers = 2), seed = 12)
  )
  
  # Check result of ncomp - only for perf not for perf.assess
  ground.ncomp <- matrix(c(2,2,2,2,2,2), ncol = 3, byrow=T,
                         dimnames = list(c("overall", "BER"),
                                         c("max.dist", "centroids.dist", "mahalanobis.dist")))
  expect_equal(out.1$choice.ncomp, ground.ncomp)
  expect_equal(out.2$choice.ncomp, ground.ncomp)
  expect_equal(out.3$choice.ncomp, NULL)
  expect_equal(out.4$choice.ncomp, NULL)
  
  # Expect the same errors
  .expect_numerically_close(out.1$error.rate.all$overall$max.dist[2,1], 0.390625)
  expect_equal(out.1$error.rate$overall, out.2$error.rate$overall)
  expect_equal(out.3$error.rate$overall, out.4$error.rate$overall)
  expect_equal(as.vector(out.1$error.rate$overall[2,]), as.vector(out.3$error.rate$overall))
  
})

# ------------------------------------------------------------------------ ##
# Test perf.mixo_splsda() and perf.assess.mixo_splsda()

test_that("perf.assess.mixo_splsda works same as perf for same components in serial and in parallel", code = {

  # set up data and model
  data(srbct)
  X <- srbct$gene
  Y <- srbct$class
  srbct.splsda <- splsda(X, Y, ncomp = 2, keepX = rep(10, 2))
  
  # Original perf() in serial execution
  out.1 <- suppressWarnings(
    perf(srbct.splsda, validation = "Mfold", folds = 2, nrepeat = 1,
                BPPARAM = SerialParam(RNGseed = 1000), seed = 12)
  )

  # Original perf() in parallel execution
  out.2 <- suppressWarnings(
    perf(srbct.splsda, validation = "Mfold", folds = 2, nrepeat = 1,
              BPPARAM = SnowParam(workers = 2), seed = 12)
  )

  # New perf.assess() in serial execution
  out.3 <- suppressWarnings(
    perf.assess(srbct.splsda, validation = "Mfold", folds = 2, nrepeat = 1,
         BPPARAM = SerialParam(RNGseed = 1000), seed = 12)
  )
  
  # New perf.assess() in parallel execution
  out.4 <- suppressWarnings(
    perf.assess(srbct.splsda, validation = "Mfold", folds = 2, nrepeat = 1,
                BPPARAM = SnowParam(workers = 2), seed = 12)
  )

  # Check result of ncomp - only for perf not for perf.assess
  ground.ncomp <- matrix(c(2,2,2,2,2,2), ncol = 3, byrow=T,
                         dimnames = list(c("overall", "BER"),
                                         c("max.dist", "centroids.dist", "mahalanobis.dist")))
  expect_equal(out.1$choice.ncomp, ground.ncomp)
  expect_equal(out.2$choice.ncomp, ground.ncomp)
  expect_equal(out.3$choice.ncomp, NULL)
  expect_equal(out.4$choice.ncomp, NULL)

  # Expect the same errors
  .expect_numerically_close(out.1$error.rate.all$overall$max.dist[2,1], 0.2380952)
  expect_equal(out.1$error.rate$overall, out.2$error.rate$overall)
  expect_equal(out.3$error.rate$overall, out.4$error.rate$overall)
  expect_equal(as.vector(out.1$error.rate$overall[2,]), as.vector(out.3$error.rate$overall))

})
