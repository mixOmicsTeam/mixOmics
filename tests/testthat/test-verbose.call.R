test_that("verbose.call works on pca", code ={
  data("nutrimouse")
  data1 <- nutrimouse$gene
  
  optimal.ncomp <- 5
  
  obj <- suppressWarnings(pca(X=data1, ncomp = optimal.ncomp, verbose.call=T))
  expect_equal(5, obj$call$ncomp)
})

test_that("verbose.call works on spca", code ={
  data("nutrimouse")
  data1 <- nutrimouse$gene
  
  optimal.ncomp <- 5
  obj <- suppressWarnings(spca(X=data1, ncomp = optimal.ncomp, verbose.call=T))
  expect_equal(5, obj$call$ncomp)
})

test_that("verbose.call works on pls", code ={
  data("nutrimouse")
  data1 <- nutrimouse$gene
  data2 <- nutrimouse$lipid
  
  optimal.ncomp <- 5
  
  obj <- suppressWarnings(pls(X=data1, Y=data2, ncomp = optimal.ncomp, verbose.call=T))
  expect_equal(5, obj$call$ncomp)
})

test_that("verbose.call works on spls", code ={
  data("nutrimouse")
  data1 <- nutrimouse$gene
  data2 <- nutrimouse$lipid
  
  optimal.ncomp <- 5
  
  obj <- suppressWarnings(spls(X=data1, Y=data2, ncomp = optimal.ncomp, verbose.call=T))
  expect_equal(5, obj$call$ncomp)
})

test_that("verbose.call works on rcca", code ={
  data("nutrimouse")
  data1 <- nutrimouse$gene
  data2 <- nutrimouse$lipid
  
  optimal.ncomp <- 5
  
  obj <- suppressWarnings(rcc(X=data1, Y=data2, ncomp = optimal.ncomp, 
             verbose.call=T, method = "shrinkage"))
  expect_equal(5, obj$call$ncomp)
})


test_that("verbose.call works on block.pls", code ={
  data("nutrimouse")
  data1 <- nutrimouse$gene
  data2 <- nutrimouse$lipid
  
  optimal.ncomp <- 5
  
  obj <- suppressWarnings(block.pls(X=list(data1 = data1,
                          data2 = data2), indY=2, ncomp = optimal.ncomp, 
                   verbose.call=T))
  expect_equal(5, obj$call$ncomp)
})


test_that("verbose.call works on block.spls", code ={
  data("nutrimouse")
  data1 <- nutrimouse$gene
  data2 <- nutrimouse$lipid
  
  optimal.ncomp <- 5
  
  obj <- suppressWarnings(block.spls(X=list(data1 = data1,
                           data2 = data2), indY=2, ncomp = optimal.ncomp, 
                    verbose.call=T))
  expect_equal(5, obj$call$ncomp)
})

test_that("verbose.call works on block.(s)plsda", code ={
  data("nutrimouse")
  data1 <- nutrimouse$gene
  data2 <- nutrimouse$lipid
  
  optimal.ncomp <- 5
  
  obj <- obj <- block.plsda(X=list(data1 = data1,
                                   data2 = data2), Y = c(rep(1,20), rep(2,20)), 
                            ncomp = optimal.ncomp, verbose.call=T)
  expect_equal(5, obj$call$ncomp)
  
  obj <- obj <- block.plsda(X=list(data1 = data1,
                                   data2 = data2), Y = c(rep(1,20), rep(2,20)), 
                            ncomp = optimal.ncomp, verbose.call=T)
  expect_equal(5, obj$call$ncomp)
})



test_that("verbose.call works on mint.pca", code ={
  data(stemcells)
  
  optimal.ncomp <- 5
  
  obj <- mint.pca(
    X = stemcells$gene,
    ncomp = optimal.ncomp,
    study = stemcells$study, 
    verbose.call=T
  )
  expect_equal(5, obj$call$ncomp)
})


test_that("verbose.call works on mint.(s)pls", code ={
  data(stemcells)
  
  optimal.ncomp <- 5
  
  obj <- mint.pls(
    X = stemcells$gene,
    Y = stemcells$gene,
    ncomp = optimal.ncomp,
    study = stemcells$study, 
    verbose.call=T
  )
  expect_equal(5, obj$call$ncomp)
  
  obj <- mint.spls(
    X = stemcells$gene,
    Y = stemcells$gene,
    ncomp = optimal.ncomp ,
    keepX = c(5, 10, 15),
    study = stemcells$study, 
    verbose.call=T
  )
  expect_equal(5, obj$call$ncomp)
})

test_that("verbose.call works on mint.(s)plsda", code ={
  data(stemcells)
  
  optimal.ncomp <- 5
  
  obj <- mint.plsda(
    X = stemcells$gene,
    Y = stemcells$celltype,
    ncomp = optimal.ncomp,
    study = stemcells$study, 
    verbose.call=T
  )
  expect_equal(5, obj$call$ncomp)
  
  obj <- mint.splsda(
    X = stemcells$gene,
    Y = stemcells$celltype,
    ncomp = optimal.ncomp ,
    keepX = c(5, 10, 15),
    study = stemcells$study, 
    verbose.call=T
  )
  expect_equal(5, obj$call$ncomp)
})