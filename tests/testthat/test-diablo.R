
###############################################################################
### ============================ GROUND TRUTHS ============================ ###
###############################################################################

Test.Data <- readRDS(system.file("testdata", "testdata-diablo.rda", package = "mixOmics"))
Testable.Components <- Test.Data$tc
Ground.Truths <- Test.Data$gt

###############################################################################
### ================================ BASIC ================================ ###
###############################################################################


test_that("(diablo:basic): block.splsda", {
    
  testable.components <- Testable.Components$basic.block.splsda
  GT <- Ground.Truths$basic.block.splsda
  
  set.seed(16)
  samples <- sample(1:150, 5)
  
  data(breast.TCGA)
  X = list(miRNA = breast.TCGA$data.train$mirna[samples,1:10],
           mRNA = breast.TCGA$data.train$mrna[samples,1:10],
           proteomics = breast.TCGA$data.train$protein[samples,1:10])
  Y = breast.TCGA$data.train$subtype[samples]
  
  choice.keepX = list(miRNA = c(3,3),
                      mRNA = c(3,3),
                      proteomics = c(3,3))
  
  res.block.splsda <- block.splsda(X, Y, keepX = choice.keepX)
    
  invisible(capture.output(TT <- dput(res.block.splsda[testable.components])))
    
  expect_equal(TT, GT)
})


test_that("(diablo:basic): wrapper.sgccda", {
    
  testable.components <- Testable.Components$basic.wrapper.sgccda
  GT <- Ground.Truths$basic.wrapper.sgccda

  set.seed(16)
  samples <- sample(1:150, 5)
  
  data(breast.TCGA)
  X = list(miRNA = breast.TCGA$data.train$mirna[samples,1:10],
           mRNA = breast.TCGA$data.train$mrna[samples,1:10],
           proteomics = breast.TCGA$data.train$protein[samples,1:10])
  Y = breast.TCGA$data.train$subtype[samples]
  
  choice.keepX = list(miRNA = c(3,3),
                      mRNA = c(3,3),
                      proteomics = c(3,3))

  res.sgccda <- wrapper.sgccda(X, Y, keepX = choice.keepX)
    
  invisible(capture.output(TT <- dput(res.sgccda[testable.components])))
    
  expect_equal(TT, GT)
})


###############################################################################
### ================================ DATA ================================= ###
###############################################################################


test_that("(diablo:data): breast.test", {
    
  testable.components <- Testable.Components$breast.test.block.splsda
  GT <- Ground.Truths$breast.test.block.splsda
  
  set.seed(16)
  samples <- sample(1:70, 5)
  
  data(breast.TCGA)
  X = list(miRNA = breast.TCGA$data.test$mirna[samples,1:10],
           mRNA = breast.TCGA$data.test$mrna[samples,1:10])
  Y = breast.TCGA$data.test$subtype[samples]
  
  choice.keepX = list(miRNA = c(3,3),
                      mRNA = c(3,3))

  breast.test.res.block.splsda <- block.splsda(X, Y, keepX = choice.keepX)
    
  invisible(capture.output(TT <- dput(breast.test.res.block.splsda[testable.components])))
    
  expect_equal(TT, GT)
})


###############################################################################
### ============================== PARAMETER ============================== ###
###############################################################################


test_that("(diablo:parameter): indY", {
    
  testable.components <- Testable.Components$indY.block.splsda
  GT <- Ground.Truths$indY.block.splsda
  
  set.seed(16)
  samples <- sample(1:70, 5)

  data(breast.TCGA)
  X = list(miRNA = breast.TCGA$data.train$mirna[samples,1:10],
           mRNA = breast.TCGA$data.train$mrna[samples,1:10],
           subtype = breast.TCGA$data.train$subtype[samples])

  choice.keepX = list(miRNA = c(3,3),
                      mRNA = c(3,3))

  res.indY.block.splsda <- block.splsda(X, indY=3, keepX = choice.keepX)
    
  invisible(capture.output(TT <- dput(res.indY.block.splsda[testable.components])))
    
  expect_equal(TT, GT)
})


test_that("(diablo:parameter): ncomp", {
    
  testable.components <- Testable.Components$ncomp.block.splsda
  GT <- Ground.Truths$ncomp.block.splsda

  set.seed(16)
  samples <- sample(1:150, 5)
  
  data(breast.TCGA)
  X = list(miRNA = breast.TCGA$data.train$mirna[samples,1:10],
           mRNA = breast.TCGA$data.train$mrna[samples,1:10],
           proteomics = breast.TCGA$data.train$protein[samples,1:10])
  Y = breast.TCGA$data.train$subtype[samples]
  
  choice.keepX = list(miRNA = c(3,3),
                      mRNA = c(3,3),
                      proteomics = c(3,3))

  ncomp.res.block.splsda <- block.splsda(X, Y, keepX = choice.keepX,
                                            ncomp = 3)
    
  invisible(capture.output(TT <- dput(ncomp.res.block.splsda[testable.components])))
    
  expect_equal(TT, GT)
})


test_that("(diablo:parameter): design", {
    
  testable.components <- Testable.Components$design.block.splsda

  set.seed(16)
  samples <- sample(1:150, 5)
  
  data(breast.TCGA)
  X = list(miRNA = breast.TCGA$data.train$mirna[samples,1:10],
           mRNA = breast.TCGA$data.train$mrna[samples,1:10],
           proteomics = breast.TCGA$data.train$protein[samples,1:10])
  Y = breast.TCGA$data.train$subtype[samples]
  
  choice.keepX = list(miRNA = c(3,3),
                      mRNA = c(3,3),
                      proteomics = c(3,3))

  design.0.res.block.splsda <- block.splsda(X, Y, keepX = choice.keepX,
                                            design = 0)
  design.null.res.block.splsda <- block.splsda(X, Y, keepX = choice.keepX,
                                               design = "null")

  design.0.5.res.block.splsda <- block.splsda(X, Y, keepX = choice.keepX,
                                            design = 0.5)
  
  design.1.res.block.splsda <- block.splsda(X, Y, keepX = choice.keepX,
                                            design = 1)
  design.full.res.block.splsda <- block.splsda(X, Y, keepX = choice.keepX,
                                            design = "full")
    
  invisible(capture.output(TT <- dput(design.0.res.block.splsda[testable.components])))
  GT <- Ground.Truths$design.0.block.splsda
  expect_equal(TT, GT)
  invisible(capture.output(expect_equal(TT, dput(design.null.res.block.splsda[testable.components]))))
  
  invisible(capture.output(TT <- dput(design.0.5.res.block.splsda[testable.components])))
  GT <- Ground.Truths$design.0.5.block.splsda
  expect_equal(TT, GT)
    
  invisible(capture.output(TT <- dput(design.1.res.block.splsda[testable.components])))
  GT <- Ground.Truths$design.1.block.splsda
  expect_equal(TT, GT)
  invisible(capture.output(expect_equal(TT, dput(design.full.res.block.splsda[testable.components]))))
})


test_that("(diablo:parameter): scheme", {
    
  testable.components <- Testable.Components$scheme.block.splsda
  GT <- Ground.Truths$scheme.block.splsda

  set.seed(16)
  samples <- sample(1:150, 5)
  
  data(breast.TCGA)
  X = list(miRNA = breast.TCGA$data.train$mirna[samples,1:10],
           mRNA = breast.TCGA$data.train$mrna[samples,1:10],
           proteomics = breast.TCGA$data.train$protein[samples,1:10])
  Y = breast.TCGA$data.train$subtype[samples]
  
  choice.keepX = list(miRNA = c(3,3),
                      mRNA = c(3,3),
                      proteomics = c(3,3))

  scheme.res.block.splsda <- block.splsda(X, Y, keepX = choice.keepX,
                                         scheme = "factorial")
    
  invisible(capture.output(TT <- dput(scheme.res.block.splsda[testable.components])))
    
  expect_equal(TT, GT)
})


test_that("(diablo:parameter): near.zero.var", {
    
  testable.components <- Testable.Components$nzv.block.splsda
  GT <- Ground.Truths$nzv.block.splsda

  create.many.zeroes <- function(x) {
    for (i in 1:length(x)) {
      set.seed(16)
      x[[i]][, 1:7] <- 0
    }
    return(x)
  }
    
  set.seed(16)
  samples <- sample(1:150, 10)
  
  data(breast.TCGA)
  X = list(miRNA = breast.TCGA$data.train$mirna[samples,1:10],
           mRNA = breast.TCGA$data.train$mrna[samples,1:10])
  Y = breast.TCGA$data.train$subtype[samples]
  
  choice.keepX = list(miRNA = c(3,3),
                      mRNA = c(3,3))

  X <- create.many.zeroes(X)

  nzv.res.block.splsda <- suppressWarnings(block.splsda(X, Y, keepX = choice.keepX,
                                                        near.zero.var = T))
    
  invisible(capture.output(TT <- dput(nzv.res.block.splsda[testable.components])))
    
  expect_equal(TT, GT)
})


test_that("(diablo:parameter): all.outputs", {

  data(breast.TCGA)
  X = list(miRNA = breast.TCGA$data.train$mirna,
           mRNA = breast.TCGA$data.train$mrna,
           proteomics = breast.TCGA$data.train$protein)
  Y = breast.TCGA$data.train$subtype

  choice.keepX = list(miRNA = c(10,10),
                      mRNA = c(10,10),
                      proteomics = c(10,10))

  all.outputs.res.block.splsda <- block.splsda(X, Y, keepX = choice.keepX,
                                          all.outputs = F)
    
  expect_equal(all.outputs.res.block.splsda$AVE, NULL)
  expect_equal(all.outputs.res.block.splsda$prop_expl_var, NULL)
})


###############################################################################
### ================================ ERROR ================================ ###
###############################################################################


test_that("(diablo:error): ensure row names are the same", {
  
  set.seed(16)
  samples <- sample(1:150, 5)
  
  data(breast.TCGA)
  X = list(miRNA = breast.TCGA$data.train$mirna[samples,1:10],
           mRNA = breast.TCGA$data.train$mrna[samples[c(5,4,3,2,1)],1:10],
           proteomics = breast.TCGA$data.train$protein[samples,1:10])
  Y = breast.TCGA$data.train$subtype[samples]
  
  choice.keepX = list(miRNA = c(3,3),
                      mRNA = c(3,3),
                      proteomics = c(3,3))
  
  expect_error(block.splsda(X, Y, keepX = choice.keepX),
               "Please check the rownames of the data, there seems to be some
    discrepancies",
    fixed=T)
})


test_that("(diablo:error): ensure Y is a class vector/factor", {
  
  set.seed(16)
  samples <- sample(1:150, 5)
  
  data(breast.TCGA)
  X = list(miRNA = breast.TCGA$data.train$mirna[samples,1:10],
           mRNA = breast.TCGA$data.train$mrna[samples,1:10],
           proteomics = breast.TCGA$data.train$protein[samples,1:10])
  Y = breast.TCGA$data.train$mrna[samples, ]
  
  expect_error(block.splsda(X, Y),
               "'Y' should be a factor or a class vector.",
    fixed=T)
})


test_that("(diablo:error): ensure each block has unique name", {
  
  set.seed(16)
  samples <- sample(1:150, 5)
  
  data(breast.TCGA)
  X = list(miRNA = breast.TCGA$data.train$mirna[samples,1:10],
           miRNA = breast.TCGA$data.train$mrna[samples,1:10],
           proteomics = breast.TCGA$data.train$protein[samples,1:10])
  Y = breast.TCGA$data.train$subtype[samples]
  
  expect_error(block.splsda(X, Y),
               "Each block of 'X' must have a unique name.",
    fixed=T)
})


test_that("(diablo:error): ensure each block has same number of rows", {
  
  set.seed(16)
  samples <- sample(1:150, 5)
  
  data(breast.TCGA)
  X = list(miRNA = breast.TCGA$data.train$mirna[samples,1:10],
           mRNA = breast.TCGA$data.train$mrna[samples[c(1,2,3)],1:10],
           proteomics = breast.TCGA$data.train$protein[samples,1:10])
  Y = breast.TCGA$data.train$subtype[samples]
  
  expect_error(block.splsda(X, Y),
               "Unequal number of rows among the blocks of 'X'",
    fixed=T)
})


test_that("(diablo:error): ensure design is in the right format", {
  
  set.seed(16)
  samples <- sample(1:150, 5)
  
  data(breast.TCGA)
  X = list(miRNA = breast.TCGA$data.train$mirna[samples,1:10],
           mRNA = breast.TCGA$data.train$mrna[samples,1:10],
           proteomics = breast.TCGA$data.train$protein[samples,1:10])
  Y = breast.TCGA$data.train$subtype[samples]
  
  d <- matrix(c(0,1,1,0), nrow = 2)
  
  expect_error(block.splsda(X, Y, design=d),
               "'design' must be a square matrix with 3columns.",
    fixed=T)
})


###############################################################################
### ============================== EDGE CASES ============================= ###
###############################################################################


test_that("(diablo:error): notify user of ignoring indY", {
  
  set.seed(16)
  samples <- sample(1:150, 5)
  
  data(breast.TCGA)
  X = list(miRNA = breast.TCGA$data.train$mirna[samples,1:10],
           mRNA = breast.TCGA$data.train$mrna[samples,1:10],
           proteomics = breast.TCGA$data.train$protein[samples,1:10])
  Y = breast.TCGA$data.train$subtype[samples]
  
  expect_warning(block.splsda(X, Y, indY=3),
               "'Y' and 'indY' are provided, 'Y' is used.",
    fixed=T)
})


test_that("(diablo:error): notify user of automatic lowering of ncomp", {
  
  set.seed(16)
  samples <- sample(1:150, 5)
  
  data(breast.TCGA)
  X = list(miRNA = breast.TCGA$data.train$mirna[samples,1:10],
           mRNA = breast.TCGA$data.train$mrna[samples,1:10],
           proteomics = breast.TCGA$data.train$protein[samples,1:10])
  Y = breast.TCGA$data.train$subtype[samples]
  
  expect_warning(block.splsda(X, Y, ncomp = 11),
               "Reset maximum number of variates 'ncomp[1]'
            to ncol(X[[1]])= 10.",
    fixed=T)
})
