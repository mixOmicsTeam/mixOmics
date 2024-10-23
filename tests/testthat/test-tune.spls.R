context("tune.spls")
library(BiocParallel)

test_that("tune.spls works ", code = {
    data("nutrimouse")
    X <- nutrimouse$gene
    Y <- nutrimouse$lipid
    
    set.seed(42)
    tune.spls.res = suppressWarnings(tune.spls(X, Y, ncomp = 3,
                              test.keepX = seq(5, 10, 5),
                              test.keepY = seq(3, 6, 3), measure = "cor",
                              folds = 5, nrepeat = 3, progressBar = FALSE,
                              BPPARAM = SerialParam(RNGseed = 5212)))
    
    expect_equal(class(tune.spls.res), c("tune.pls", "tune.spls"))
    expect_equal(unname(tune.spls.res$choice.keepX), c(5,5,5))
    expect_equal(unname(tune.spls.res$choice.keepY), c(3,3,6))
})

test_that("tune.spls works in parallel", code = {

    data("nutrimouse")
    X <- nutrimouse$gene
    Y <- nutrimouse$lipid
    
    ## added this to avoid errors where num_workers exceeded limits set by devtools::check()
    chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
    if (nzchar(chk) && chk == "TRUE") {
      # use 2 cores in CRAN/Travis/AppVeyor
      num_workers <- 2L
    } else {
      # use all cores in devtools::test()
      num_workers <- parallel::detectCores()
    }

    set.seed(42)
    tune.spls.res = suppressWarnings(tune.spls(X, Y, ncomp = 3,
                                      test.keepX = seq(1, 5, 1),
                                      test.keepY = seq(3, 6, 3), measure = "cor",
                                      folds = 5, nrepeat = 3, progressBar = FALSE,
                                      BPPARAM = BiocParallel::SnowParam(RNGseed = 5212, workers = num_workers)))
    
    expect_equal(class(tune.spls.res), c("tune.pls", "tune.spls"))
    expect_equal(unname(tune.spls.res$choice.keepX), c(1,1,1))
    expect_equal(unname(tune.spls.res$choice.keepY), c(3,3,3))
})

test_that("tune.spls works faster in parallel", {
  data("nutrimouse")
  X <- nutrimouse$gene
  Y <- nutrimouse$lipid
  set.seed(42)
  test.keepX <- c(1, 2, 3)
  test.keepY <- c(1, 2, 3)
  # Serial execution
  serial_time <- system.time(
    tune.spls.res.serial <- suppressWarnings(tune.spls(X, Y, ncomp = 3, test.keepX = test.keepX,
                                               test.keepY = test.keepY, measure = "cor",
                                               folds = 5, nrepeat = 20, progressBar = FALSE,
                                               BPPARAM = BiocParallel::SerialParam(RNGseed = 5212)))
  )
  # Serial execution - 2 cores
  parallel_time_2_cores <- system.time(
    tune.spls.res.parallel <- suppressWarnings(tune.spls(X, Y, ncomp = 3,test.keepX = test.keepX,
                                                       test.keepY = test.keepY, measure = "cor",
                                                       folds = 5, nrepeat = 20, progressBar = FALSE,
                                                       BPPARAM = BiocParallel::SnowParam(workers = 2, RNGseed = 5212)))
  )
  expect_true(serial_time[3] > parallel_time_2_cores[3])
})

test_that("tune.spls and tune(method='spls') are equivalent", {
    
    data("nutrimouse")
    X <- nutrimouse$gene
    Y <- nutrimouse$lipid
    
    set.seed(42)
    tune.spls.res.1 = suppressWarnings(tune.spls(X, Y, ncomp = 3,
                                      test.keepX = seq(1, 2, 1),
                                      test.keepY = seq(3, 6, 3),
                                      folds = 5, nrepeat = 3, progressBar = FALSE,
                                      BPPARAM = SerialParam(RNGseed = 5212)))
    
    tune.spls.res.2 = suppressWarnings(tune(method = "spls", X, Y, ncomp = 3,
                                      test.keepX = seq(1, 2, 1),
                                      test.keepY = seq(3, 6, 3),
                                      folds = 5, nrepeat = 3, progressBar = FALSE,
                                      BPPARAM = SerialParam(RNGseed = 5212)))
    
    expect_equal(tune.spls.res.1$measure.pred, tune.spls.res.2$measure.pred)
})

## If ncol(Y) == 1 tune.spls calls tune.spls1
## checking these cases where spls1 method is called
# test setup copied from vignette 04-spls1-tuning
# for tune.spls1 need to have seed set globally and in BPPARAM

test_that("tune.spls.1 works and is reproducible", code = {
  data(liver.toxicity)
  X <- liver.toxicity$gene
  y <- liver.toxicity$clinic[, "ALB.g.dL."]
  list.keepX <- c(5:10, seq(15, 50, 5))     
  
  set.seed(33)
  tune.spls1.MAE.1 <- tune.spls(X, y, ncomp= 2, 
                              test.keepX = list.keepX, 
                              validation = 'Mfold', 
                              folds = 7,
                              nrepeat = 7, 
                              progressBar = FALSE, 
                              measure = 'MAE',
                              BPPARAM = SerialParam(RNGseed = 33))
  set.seed(33)
  tune.spls1.MAE.2 <- tune.spls(X, y, ncomp= 2, 
                              test.keepX = list.keepX, 
                              validation = 'Mfold', 
                              folds = 7,
                              nrepeat = 7, 
                              progressBar = FALSE, 
                              measure = 'MAE',
                              BPPARAM = SerialParam(RNGseed = 33))
  
  expect_equal(tune.spls1.MAE.1$choice.keepX, tune.spls1.MAE.2$choice.keepX)
  expect_equal(class(tune.spls1.MAE.1), "tune.spls1")
  expect_equal(unname(tune.spls1.MAE.1$choice.keepX), c(20, 45))
})

test_that("tune.spls.1 and tune(method='spls') are equivalent", code = {
  data(liver.toxicity)
  X <- liver.toxicity$gene
  y <- liver.toxicity$clinic[, "ALB.g.dL."]
  list.keepX <- c(5:10, seq(15, 50, 5))     
  
  set.seed(33)
  tune.spls1.MAE.1 <- tune.spls(X, y, ncomp= 2, 
                              test.keepX = list.keepX, 
                              validation = 'Mfold', 
                              folds = 6,
                              nrepeat = 3, 
                              progressBar = FALSE, 
                              measure = 'MAE',
                              BPPARAM = SerialParam(RNGseed = 33))
  
  expect_equal(class(tune.spls1.MAE.1), "tune.spls1")
  expect_equal(unname(tune.spls1.MAE.1$choice.keepX), c(20, 25))
  
  set.seed(33)
  tune.spls1.MAE.2 <- tune(method = "spls",
                           X, y, ncomp= 2, 
                           test.keepX = list.keepX, 
                           validation = 'Mfold', 
                           folds = 6,
                           nrepeat = 3, 
                           progressBar = FALSE, 
                           measure = 'MAE',
                           BPPARAM = SerialParam(RNGseed = 33))
  expect_equal(class(tune.spls1.MAE.2), "tune.spls1")
  expect_equal(unname(tune.spls1.MAE.2$choice.keepX), c(20, 25))
})

test_that("tune.spls.1 works in serial and in parallel", code = {
  data(liver.toxicity)
  X <- liver.toxicity$gene
  y <- liver.toxicity$clinic[, "ALB.g.dL."]
  list.keepX <- c(5:10, seq(15, 50, 5))     
  
  set.seed(22)
  tune.spls1.MAE.1 <- tune.spls(X, y, ncomp= 2, 
                                test.keepX = list.keepX, 
                                validation = 'Mfold', 
                                folds = 10,
                                nrepeat = 10, 
                                progressBar = FALSE, 
                                measure = 'MAE',
                                BPPARAM = SerialParam(RNGseed = 22))
  set.seed(22)
  # avoids multiple warnings of 'package:stats' may not be available when loading' 
  tune.spls1.MAE.2 <- suppressWarnings(
    tune.spls(X, y, ncomp= 2, test.keepX = list.keepX, validation = 'Mfold', 
              folds = 10, nrepeat = 10, progressBar = FALSE, measure = 'MAE',
              BPPARAM = SnowParam(workers = 2, RNGseed = 2, exportglobals = TRUE))
  )
  
  expect_equal(tune.spls1.MAE.1$choice.keepX, tune.spls1.MAE.2$choice.keepX)
  expect_equal(class(tune.spls1.MAE.1), "tune.spls1")
  expect_equal(unname(tune.spls1.MAE.1$choice.keepX), c(20, 40))
})


