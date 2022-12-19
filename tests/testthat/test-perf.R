context("perf")

test_that("perf works when nrepeat < 3 and validation != 'loo'", code = {
  
  data(srbct) # extract the small round bull cell tumour data
  X <- srbct$gene # use the gene expression data as the X matrix
  Y <- srbct$class # use the class data as the Y matrix
  
  initial.plsda <- plsda(X, Y, ncomp = 5)
  
  set.seed(12)
  plsda.perf <- suppressWarnings(perf(initial.plsda, progressBar = FALSE, auc = FALSE, 
                     folds = 3, nrepeat = 1))
  
  trueVals <- matrix(5, ncol = 3, nrow = 2)
  colnames(trueVals) <- c("max.dist", "centroids.dist", "mahalanobis.dist")
  rownames(trueVals) <- c("overall", "BER")
  
  expect_equal(plsda.perf$choice.ncomp, trueVals)
})