###############################################################################
### ============================ MISSING TESTS ============================ ###
###############################################################################

# parameter - no parameter type tests as all only affect visualisation

###############################################################################
### ============================ GROUND TRUTHS ============================ ###
###############################################################################

Test.Data <- readRDS(system.file("testdata", "testdata-circosPlot.rda", package = "mixOmics"))
Ground.Truths <- Test.Data$gt


###############################################################################
### ================================ BASIC ================================ ###
###############################################################################


test_that("(circosPlot:basic): block.spls", {
    
  GT <- Ground.Truths$basic.block.spls

  data(breast.TCGA)
  X <- list(miRNA = breast.TCGA$data.train$mirna[, 1:10],
           mRNA = breast.TCGA$data.train$mrna[, 1:10])
  Y <- breast.TCGA$data.train$protein[, 1:5]

  choice.keepX <- list(miRNA = c(3,3),
                       mRNA = c(3,3))

  res.block.spls <- block.spls(X, Y,
                               keepX = choice.keepX)

  block.spls.circos <- circosPlot(res.block.spls, group = breast.TCGA$data.train$subtype,
                                     cutoff=0.0)
    
  invisible(capture.output(TT <- dput(block.spls.circos)))
    
  expect_equal(TT, GT)
})


test_that("(circosPlot:basic): block.splsda", {
    
  GT <- Ground.Truths$basic.block.splsda

  data(breast.TCGA)
  X <- list(miRNA = breast.TCGA$data.train$mirna[, 1:10],
           mRNA = breast.TCGA$data.train$mrna[, 1:10])
  Y <- breast.TCGA$data.train$subtype

  choice.keepX <- list(miRNA = c(3,3),
                       mRNA = c(3,3))

  res.block.splsda <- block.splsda(X, Y,
                                   keepX = choice.keepX)

  block.splsda.circos <- circosPlot(res.block.splsda,
                                    cutoff=0.0)
    
  invisible(capture.output(TT <- dput(block.splsda.circos)))
    
  expect_equal(TT, GT)
})


###############################################################################
### ================================ DATA ================================= ###
###############################################################################


test_that("(circosPlot:data): nutrimouse", {
    
  GT <- Ground.Truths$nutrimouse.block.splsda

  data(nutrimouse)
  Y = nutrimouse$diet
  X = list(gene = nutrimouse$gene[, 1:10],
           lipid = nutrimouse$lipid[, 1:10])

  choice.keepX <- list(gene = c(3,3),
                       lipid = c(3,3))

  res.block.splsda <- block.splsda(X, Y, keepX = choice.keepX)

  nutrimouse.circos <- circosPlot(res.block.splsda, cutoff = 0.0)
    
  invisible(capture.output(TT <- dput(nutrimouse.circos)))
    
  expect_equal(TT, GT)
})


###############################################################################
### ================================ ERROR ================================ ###
###############################################################################


test_that("(circosPlot:error): group parameter has appropriate value", {
  
  data(breast.TCGA)
  X = list(miRNA = breast.TCGA$data.train$mirna,
           mRNA = breast.TCGA$data.train$mrna)
  
  choice.keepX <- list(miRNA = c(10, 10),
                       mRNA = c(10, 10))
  
  res.block.spls <- block.spls(X, Y <- breast.TCGA$data.train$protein, 
                               keepX = choice.keepX)
  
  err <- "group must be a factor of length: nrow(object$X$Y) = 150\n"
  
  expect_error(circosPlot(res.block.spls, group = breast.TCGA$data.test$subtype,
                          cutoff=0.7),
               err,
               fixed = T)
  
  expect_error(circosPlot(res.block.spls, cutoff=0.7),
               err,
               fixed = T)
})


test_that("(circosPlot:error): ensure there is minimum of three blocks", {
  
  data(breast.TCGA)
  X = list(miRNA = breast.TCGA$data.train$mirna)
  
  choice.keepX <- list(miRNA = c(10, 10))
  
  res.block.splsda <- block.splsda(X, Y = breast.TCGA$data.train$subtype, 
                                   keepX = choice.keepX)
  
  expect_error(circosPlot(res.block.splsda,
                         cutoff=0.7),
              "This function is only available when there are more than 3 blocks
    (2 in object$X + an outcome object$Y)",
              fixed=T)
})


test_that("(circosPlot:error): ensure cutoff parameter is specified", {
  
  data(breast.TCGA)
  X = list(miRNA = breast.TCGA$data.train$mirna,
           mRNA = breast.TCGA$data.train$mrna)
  
  choice.keepX <- list(miRNA = c(10, 10),
                       mRNA = c(10, 10))
  
  res.block.splsda <- block.splsda(X, Y = breast.TCGA$data.train$subtype, 
                                   keepX = choice.keepX)
  
  expect_error(circosPlot(res.block.splsda),
               "'cutoff' is missing",
               fixed = T)
               
})


test_that("(circosPlot:error): cannot take block.pls objects", {
  
  data(breast.TCGA)
  X = list(miRNA = breast.TCGA$data.train$mirna,
           mRNA = breast.TCGA$data.train$mrna)
  
  res.block.pls <- block.pls(X, Y = breast.TCGA$data.train$protein)
  
  expect_error(circosPlot(res.block.pls, group = breast.TCGA$data.train$subtype,
                          cutoff=0.7),
               "no applicable method for 'circosPlot' applied to an object of class \"c('block.pls', 'sgcca')\"",
               fixed = T)
})


test_that("(circosPlot:error): cannot take block.plsda objects", {
  
  data(breast.TCGA)
  X = list(miRNA = breast.TCGA$data.train$mirna,
           mRNA = breast.TCGA$data.train$mrna)
  
  res.block.plsda <- block.plsda(X, Y = breast.TCGA$data.train$subtype)
  
  expect_error(circosPlot(res.block.plsda, cutoff=0.7),
               "no applicable method for 'circosPlot' applied to an object of class \"c('block.plsda', 'block.pls', 'sgccda', 'sgcca', 'DA')\"",
               fixed=T)
})


###############################################################################
### ============================== EDGE CASES ============================= ###
###############################################################################


test_that("(circosPlot:edge.case): works with similar feature names in different blocks", {
  
  create_similar_feature_names <- function(data_list)
      {
          lapply(data_list, function(x){
              colnames(x) <- paste0('feature_', seq_len(ncol(x)))
              x
          })
      }
      
  data("breast.TCGA")
  X = list(miRNA = breast.TCGA$data.train$mirna,
           mRNA = breast.TCGA$data.train$mrna,
           proteomics = breast.TCGA$data.train$protein)
  
  X <- create_similar_feature_names(X)
  
  choice.keepX <- list(miRNA = c(10, 10),
                       mRNA = c(10, 10),
                       proteomics = c(10,10))
  
  res.block.splsda = block.splsda(X = X, Y = breast.TCGA$data.train$subtype,
                                   keepX = choice.keepX)
  
  expect_output(circosPlot(res.block.splsda, cutoff = 0.7),
                "adding block name to feature names in the output similarity matrix as there are similar feature names across blocks.",
                fixed=T)
})


test_that("(circosPlot:edge.case): warning for high cutoff value", {
  
  data(breast.TCGA)
  X = list(miRNA = breast.TCGA$data.train$mirna,
           mRNA = breast.TCGA$data.train$mrna)
  
  choice.keepX <- list(miRNA = c(10, 10),
                       mRNA = c(10, 10))
  
  res.block.splsda <- block.splsda(X, Y = breast.TCGA$data.train$subtype, 
                                   keepX = choice.keepX)
  
  expect_warning(circosPlot(res.block.splsda, cutoff=0.99),
                 "Choose a lower correlation threshold to highlight
    links between datasets",
                 fixed=T)
})

