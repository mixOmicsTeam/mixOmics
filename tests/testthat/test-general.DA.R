context("general.DA")

test_that("ALL DA functions raise specific error when Y contains NAs", {
  
  data(srbct)
  X <- srbct$gene
  Y <- srbct$class
  Y[c(1,2,3)] <- NA
  
  data("breast.TCGA")
  X.b <- list(mirna = breast.TCGA$data.train$mirna,
              mrna = breast.TCGA$data.train$mrna)
  Y.b <- breast.TCGA$data.train$subtype
  Y.b[c(1,2,3)] <- NA
  Y.b.2 <- breast.TCGA$data.train$protein
  
  data("stemcells")
  X.m <- stemcells$gene
  Y.m <- stemcells$celltype
  Y.m[c(1,2,3)] <- NA
  S.m <- stemcells$study
  
  # ------------------------------------------------------------------------- #
  
  # plsda
  expect_error(plsda(X, Y), 
               "Unmapped Y contains samples with no associated class. May be caused by NAs in input Y vector", 
               fixed = TRUE)
  
  # splsda
  expect_error(plsda(X, Y), 
               "Unmapped Y contains samples with no associated class. May be caused by NAs in input Y vector", 
               fixed = TRUE)
  
  # block.plsda
  expect_error(block.plsda(X.b, Y.b), 
               "Unmapped Y contains samples with no associated class. May be caused by NAs in input Y vector", 
               fixed = TRUE)
  
  # block.splsda
  expect_error(block.splsda(X.b, Y.b), 
               "Unmapped Y contains samples with no associated class. May be caused by NAs in input Y vector", 
               fixed = TRUE)
  
  # mint.plsda
  expect_error(mint.plsda(X.m, Y.m, study = S.m), 
               "Unmapped Y contains samples with no associated class. May be caused by NAs in input Y vector", 
               fixed = TRUE)
  
  # mint.splsda
  expect_error(mint.splsda(X.m, Y.m, study = S.m), 
               "Unmapped Y contains samples with no associated class. May be caused by NAs in input Y vector", 
               fixed = TRUE)
  
  # block.pls - ensure no error is raised for regression type problem
  expect_is(block.pls(X.b, Y.b.2), "block.pls")
  
  # block.spls - ensure no error is raised for regression type problem
  expect_is(block.spls(X.b, Y.b.2), "block.spls")
})