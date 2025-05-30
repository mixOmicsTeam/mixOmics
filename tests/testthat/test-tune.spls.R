context("tune.spls")
library(BiocParallel)

test_that("tune.spls works and is the same in parallel and when run in tune wrapper", code = {
  
  # set up data
  data("nutrimouse")
  X <- nutrimouse$gene
  Y <- nutrimouse$lipid
  
  # run in serial
  tune.spls.res.1 = suppressWarnings(tune.spls(X, Y, ncomp = 2,
                                             test.keepX = seq(5, 10, 5), 
                                             test.keepY = seq(3, 6, 3), measure = "cor",
                                             logratio = "none", mode = "regression", 
                                             folds = 2, nrepeat = 1, progressBar = FALSE,
                                             BPPARAM = SerialParam(RNGseed = 5), # RNGseed is ignored
                                             seed = 5212))
  
  # run in parallel
  tune.spls.res.2 = suppressWarnings(tune.spls(X, Y, ncomp = 2,
                                               test.keepX = seq(5, 10, 5), 
                                               test.keepY = seq(3, 6, 3), measure = "cor",
                                               logratio = "none", mode = "regression", 
                                               folds = 2, nrepeat = 1, progressBar = FALSE,
                                               BPPARAM = SnowParam(workers = 2),
                                               seed = 5212))
  
  # in tune wrapper in serial
  tune.spls.res.3 = suppressWarnings(tune(X, Y, ncomp = 2,
                                               test.keepX = seq(5, 10, 5), 
                                               test.keepY = seq(3, 6, 3), measure = "cor",
                                               logratio = "none", mode = "regression", 
                                               folds = 2, nrepeat = 1, progressBar = FALSE,
                                               BPPARAM = SerialParam(RNGseed = NULL),
                                          seed = 5212,
                                          method = "spls"),
                                     )

  # in tune wrapper in parallel
  tune.spls.res.4 = suppressWarnings(tune(X, Y, ncomp = 2,
                                          test.keepX = seq(5, 10, 5), 
                                          test.keepY = seq(3, 6, 3), measure = "cor",
                                          logratio = "none", mode = "regression", 
                                          folds = 2, nrepeat = 1, progressBar = FALSE,
                                          BPPARAM = SnowParam(workers = 2),
                                          method = "spls",
                                          seed = 5212),
  )
  
  
  # check outputs
  expect_equal(class(tune.spls.res.1), class(tune.spls.res.2), class(tune.spls.res.3), class(tune.spls.res.4), "tune.spls")
  expect_equal(unname(tune.spls.res.1$choice.keepX), c(10,5))
  expect_equal(unname(tune.spls.res.2$choice.keepX), c(10,5))
  expect_equal(unname(tune.spls.res.3$choice.keepX), c(10,5))
  expect_equal(unname(tune.spls.res.4$choice.keepX), c(10,5))
  expect_equal(unname(tune.spls.res.1$choice.keepY), c(3,6))
  expect_equal(unname(tune.spls.res.2$choice.keepY), c(3,6))
  expect_equal(unname(tune.spls.res.3$choice.keepY), c(3,6))
  expect_equal(unname(tune.spls.res.4$choice.keepY), c(3,6))
  
  # check outputs exactly the same regardless of how the function was run
  .expect_numerically_close(tune.spls.res.1$measure.pred$mean, tune.spls.res.2$measure.pred$mean)
  .expect_numerically_close(tune.spls.res.1$measure.pred$mean, tune.spls.res.3$measure.pred$mean)
  .expect_numerically_close(tune.spls.res.1$measure.pred$mean, tune.spls.res.4$measure.pred$mean)
  
  # check can plot outputs
  expect_silent(plot(tune.spls.res.1))
  expect_silent(plot(tune.spls.res.2))
  expect_silent(plot(tune.spls.res.3))
  expect_silent(plot(tune.spls.res.4))
})

## If ncol(Y) == 1 tune.spls calls tune.spls1
## checking these cases where spls1 method is called
# test setup copied from vignette 04-spls1-tuning
# for tune.spls1 need to have seed set globally and in BPPARAM

test_that("tune.spls.1 works and is the same in parallel and when run in tune wrapper", code = {
  
  # set up data
  data(liver.toxicity)
  X <- liver.toxicity$gene
  y <- liver.toxicity$clinic[, "ALB.g.dL."]
  list.keepX <- c(5:10, seq(15, 50, 5))     
  
  # run in serial
  tune.spls1.MAE.1 <- tune.spls(X, y, ncomp = 2, test.keepX = list.keepX, 
                                validation = 'Mfold', folds = 2, nrepeat = 1, 
                                progressBar = FALSE, measure = 'MAE',
                                BPPARAM = SerialParam(RNGseed = NULL), seed = 33)
  
  # run in parallel
  tune.spls1.MAE.2 <- tune.spls(X, y, ncomp = 2, test.keepX = list.keepX, 
                           validation = 'Mfold', folds = 2, nrepeat = 1, 
                           progressBar = FALSE, measure = 'MAE',
                           BPPARAM = SnowParam(workers = 2), seed = 33)
  
  # in tune wrapper in serial
  tune.spls1.MAE.3 <- tune(X, y, ncomp = 2, test.keepX = list.keepX, 
                           validation = 'Mfold', folds = 2, nrepeat = 1, 
                           progressBar = FALSE, measure = 'MAE',
                           BPPARAM = SerialParam(),
                           method = "spls", seed = 33)
  
  # in tune wrapper in parallel
  tune.spls1.MAE.4 <- tune(X, y, ncomp = 2, test.keepX = list.keepX, 
                           validation = 'Mfold', folds = 2, nrepeat = 1, 
                           progressBar = FALSE, measure = 'MAE',
                           BPPARAM = SnowParam(workers = 2),
                           method = "spls", seed = 33)
  
  
  # check outputs
  expect_equal(class(tune.spls1.MAE.1), class(tune.spls1.MAE.2), class(tune.spls1.MAE.3), class(tune.spls1.MAE.4), "tune.spls1")
  expect_equal(unname(tune.spls1.MAE.1$choice.keepX), c(25, 20))
  expect_equal(unname(tune.spls1.MAE.2$choice.keepX), c(25, 20))
  expect_equal(unname(tune.spls1.MAE.3$choice.keepX), c(25, 20))
  expect_equal(unname(tune.spls1.MAE.4$choice.keepX), c(25, 20))
  
  # check outputs exactly the same regardless of how the function was run
  expect_equal(tune.spls1.MAE.1$measure.pred$mean, tune.spls1.MAE.2$measure.pred$mean)
  expect_equal(tune.spls1.MAE.1$measure.pred$mean, tune.spls1.MAE.3$measure.pred$mean)
  expect_equal(tune.spls1.MAE.1$measure.pred$mean, tune.spls1.MAE.4$measure.pred$mean)
})

test_that("tune.spls works when test.keepX and test.keepY = NULL and gives same result as perf()", code = {
  
  # set up data
  data(liver.toxicity)
  X <- liver.toxicity$gene
  Y <- liver.toxicity$clinic
  
  # tune on components only
  tune_res <- suppressWarnings(
    tune.spls(X, Y, ncomp = 2, logratio = "none",
                nrepeat = 1, folds = 3, measure = "cor",
                test.keepX = NULL, test.keepY = NULL, seed = 20)
  )
  
  # run perf
  spls_res <- spls(X, Y, ncomp = 2, logratio = "none")
  perf_res <- suppressWarnings(
    perf(spls_res, ncomp = 2, nrepeat = 1, folds = 3, seed = 20)
  )
  
  # check results match
  expect_equal(tune_res$measures$MSEP$summary[1,3], perf_res$measures$MSEP$summary[1,3])
})

test_that("tune.spls works in different test cases", code = {
  
  # set up data
  data(liver.toxicity)
  X <- liver.toxicity$gene
  Y <- liver.toxicity$clinic
  list.keepX <- c(5:10, seq(15, 50, 5))     
  
  # near.zero.var = TRUE
  expect_no_error(
    tune.spls(X, Y, ncomp = 2, logratio = "none", near.zero.var = TRUE,
              nrepeat = 1, folds = 3, measure = "cor",
              test.keepX = NULL, test.keepY = NULL, seed = 20)
  )
  expect_no_error(
    tune.spls(X, Y, ncomp = 2, logratio = "none", near.zero.var = TRUE,
              nrepeat = 1, folds = 3, measure = "cor",
              test.keepX = list.keepX, test.keepY = NULL, seed = 20)
  )
  
  # canonical mode
  expect_no_error(
    tune.spls(X, Y, ncomp = 2, logratio = "none", near.zero.var = TRUE,
              nrepeat = 1, folds = 3, measure = "cor", mode = "canonical",
              test.keepX = NULL, test.keepY = NULL, seed = 20)
  )
  expect_no_error(
    tune.spls(X, Y, ncomp = 2, logratio = "none", near.zero.var = TRUE,
              nrepeat = 1, folds = 3, measure = "cor", mode = "canonical",
              test.keepX = list.keepX, test.keepY = NULL, seed = 20)
  )
})


