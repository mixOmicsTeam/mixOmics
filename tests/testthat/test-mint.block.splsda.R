

###############################################################################
### ================================ BASIC ================================ ###
###############################################################################

test_that("mint.block.plsda and mint.block.splsda works", {
  
  data(breast.TCGA)
  mrna <- rbind(breast.TCGA$data.train$mrna, breast.TCGA$data.test$mrna)
  mirna <- rbind(breast.TCGA$data.train$mirna, breast.TCGA$data.test$mirna)
  X <- list(mrna = mrna, mirna = mirna)
  Y <- c(breast.TCGA$data.train$subtype, breast.TCGA$data.test$subtype)
  
  study <- c(rep("study1",150), rep("study2",70))
  
  res.mint.block.plsda <- mint.block.plsda(X, Y, study=study, design="full")
  
  expect_true("mint.block.plsda" %in% class(res.mint.block.plsda))
  .expect_numerically_close(res.mint.block.plsda$variates$mrna[1,1], -7.503)
  .expect_numerically_close(res.mint.block.plsda$variates$Y[220,2], -0.796)
  
  .expect_numerically_close(res.mint.block.plsda$loadings$mrna[1,1], 0.100)
  .expect_numerically_close(res.mint.block.plsda$loadings$mirna[100,2], -0.019)
  
  
  res.mint.block.splsda <- mint.block.splsda(X, Y, study=study, design="full", keepX = list(mrna=c(2,2), mirna=c(3,3)))
  
  expect_true("mint.block.splsda" %in% class(res.mint.block.splsda))
  .expect_numerically_close(res.mint.block.splsda$variates$mrna[1,1], -1.48)
  .expect_numerically_close(res.mint.block.splsda$variates$Y[220,2], -0.5758)
  
})


###############################################################################
### ================================ DATA ================================= ###
###############################################################################





###############################################################################
### ============================== PARAMETER ============================== ###
###############################################################################





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