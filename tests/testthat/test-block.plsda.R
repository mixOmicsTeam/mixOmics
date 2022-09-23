

###############################################################################
### ================================ BASIC ================================ ###
###############################################################################


test_that("(block.plsda:basic): default", {
  
  data(breast.TCGA)
  X <- list(miRNA = breast.TCGA$data.train$mirna,
            mRNA = breast.TCGA$data.train$mrna,
            proteomics = breast.TCGA$data.train$protein)
  Y <- breast.TCGA$data.train$subtype
  
  res.block.plsda <- block.plsda(X, Y)
  
  expect_equal(class(res.block.plsda)[1], "block.plsda")
  
  .expect_numerically_close(res.block.plsda$variates$mRNA[1,1], 7.50)
  .expect_numerically_close(res.block.plsda$variates$Y[50,2], -0.83)
  
  .expect_numerically_close(res.block.plsda$loadings$miRNA[1,1], -0.0672)
  .expect_numerically_close(res.block.plsda$loadings$Y[3,2], -0.711)
})


###############################################################################
### ================================ DATA ================================= ###
###############################################################################





###############################################################################
### ============================== PARAMETER ============================== ###
###############################################################################


test_that("(block.plsda:parameter): indY", {
  
  data(breast.TCGA)
  X <- list(miRNA = breast.TCGA$data.train$mirna,
            mRNA = breast.TCGA$data.train$mrna,
            proteomics = breast.TCGA$data.train$protein,
            response = breast.TCGA$data.train$subtype)
  
  res.block.plsda <- block.plsda(X, indY = 4)
  
  expect_equal(class(res.block.plsda)[1], "block.plsda")
  
  .expect_numerically_close(res.block.plsda$variates$mRNA[1,1], 7.50)
  .expect_numerically_close(res.block.plsda$variates$response[50,2], -0.83)
  
  .expect_numerically_close(res.block.plsda$loadings$miRNA[1,1], -0.0672)
  .expect_numerically_close(res.block.plsda$loadings$response[3,2], -0.711)
})


###############################################################################
### ================================ ERROR ================================ ###
###############################################################################





###############################################################################
### ============================== EDGE CASES ============================= ###
###############################################################################





###############################################################################

# if the method is a graphical function, include the following:
#dev.off()
#unlink(list.files(pattern = "*.pdf"))