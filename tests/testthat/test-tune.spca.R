
test_that("tune.spca works", {
  
  
  data(srbct)
  X <- srbct$gene[1:20, 1:200]
  
  grid.keepX <- seq(5, 35, 10)
  
  set.seed(5212)
  object <- tune.spca(X,ncomp = 2, 
                      folds = 5, 
                      test.keepX = grid.keepX, nrepeat = 3,
                      BPPARAM = SerialParam(RNGseed = 5212))
  
  expect_equal(object$choice.keepX[[1]], 35)
  expect_equal(object$choice.keepX[[2]], 5)
})


test_that("tune.spca works with NA input", {
  
  data(srbct)
  X <- srbct$gene[1:20, 1:200]
  
  set.seed(5212)
  na.feats <- sample(1:ncol(X), 20)
  
  for (c in na.feats) {
    na.samples <- sample(1:nrow(X),1) #  sample.int(3, 1)
    
    X[na.samples, c] <- NA
  }
  
  grid.keepX <- seq(5, 35, 10)
  
  expect_warning({object <- tune.spca(X,ncomp = 2, 
                                      folds = 5, 
                                      test.keepX = grid.keepX, nrepeat = 3,
                                      BPPARAM = SerialParam(RNGseed = 5212))},
                 "NAs present")
  
  expect_equal(object$choice.keepX[[1]], 15)
  expect_equal(object$choice.keepX[[2]], 5)
})