context("tune.splsda")

test_that("tune.splsda works ", code = {
  data(breast.tumors)
  X = breast.tumors$gene.exp
  Y = as.factor(breast.tumors$sample$treatment)
  
  set.seed(42)
  tune.splsda.res = tune.splsda(X, Y, ncomp = 2, nrepeat = 3, logratio = "none",
                                test.keepX = c(5, 10, 15), folds = 3, dist = "max.dist",
                                BPPARAM = SerialParam(RNGseed = 100))
  
  expect_equal(class(tune.splsda.res), "tune.splsda")
  expect_equal(unname(tune.splsda.res$choice.keepX), c(5,15))
  expect_equal(tune.splsda.res$choice.ncomp$ncomp, 1)
  .expect_numerically_close(tune.splsda.res$error.rate[1,1], 0.162037)
})

test_that("tune.splsda works in serial and in parallel", {
    data(breast.tumors)
    X = breast.tumors$gene.exp
    Y = as.factor(breast.tumors$sample$treatment)
    RNGversion(.mixo_rng()) ## in case RNG changes!
    set.seed(100)
    tune = tune.splsda(X, Y, ncomp = 2, nrepeat = 3, logratio = "none",
                       test.keepX = c(5, 10, 15), folds = 3, dist = "max.dist",
                       BPPARAM = SerialParam(RNGseed = 100))
    expect_equal(tune$choice.ncomp$ncomp, 1L)
    expect_equal(tune$choice.keepX, c(comp1 = 5, comp2 = 15))
    
    set.seed(100)
    tune2 = tune.splsda(X, Y, ncomp = 2, nrepeat = 3, logratio = "none",
                       test.keepX = c(5, 10, 15), folds = 3, dist = "max.dist", 
                       BPPARAM = MulticoreParam(workers = 2, RNGseed = 100))
    expect_equal(tune$choice.keepX, tune2$choice.keepX)
    expect_equal(tune$choice.keepX, tune2$choice.keepX)
    
})

test_that("tune.splsda works faster in parallel", {
  data(liver.toxicity)
  X <- liver.toxicity$gene
  Y <- liver.toxicity$clinic
  set.seed(42)
  # Serial execution
  serial_time <- system.time(
    tune_res_serial <- tune.spls(X, Y, ncomp = 3, 
                                 test.keepX = c(5, 10, 15), 
                                 test.keepY = c(3, 6, 8), 
                                 folds = 5,
                                 BPPARAM = SerialParam())
  )
  # Serial execution - 2 cores
  parallel_time_2_cores <- system.time(
    tune_res_serial <- tune.spls(X, Y, ncomp = 3, 
                                 test.keepX = c(5, 10, 15), 
                                 test.keepY = c(3, 6, 8), 
                                 folds = 5,
                                 BPPARAM = MulticoreParam(workers = 2))
  )
  expect_true(serial_time[3] > parallel_time_2_cores[3])
  
})

test_that("tune.splsda and tune(method='splsda') are equivalent", {
  data(breast.tumors)
  X = breast.tumors$gene.exp
  Y = as.factor(breast.tumors$sample$treatment)
  set.seed(42)
  tune.splsda.res.1 = tune.splsda(X, Y, ncomp = 2, nrepeat = 3, logratio = "none",
                                test.keepX = c(5, 10, 15), folds = 3, dist = "max.dist",
                                BPPARAM = SerialParam(RNGseed = 100))
  tune.splda.res.2 = suppressWarnings(tune(method = "splsda", 
                                             X, Y, ncomp = 2, nrepeat = 3, logratio = "none",
                                             test.keepX = c(5, 10, 15), folds = 3, dist = "max.dist",
                                             BPPARAM = SerialParam(RNGseed = 100)))
  expect_equal(tune.splsda.res.1$measure.pred, tune.splda.res.2$measure.pred)
})
  