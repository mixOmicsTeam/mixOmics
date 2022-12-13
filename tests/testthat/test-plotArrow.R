test_that("plotArrow does not function on (mint).(s).plsda objects", code = {
  
  data("stemcells")
  
  X <- stemcells$gene
  Y <- stemcells$celltype
  S <- stemcells$study
  
  optimal.ncomp <- 4
  optimal.keepX <- c(24, 45, 20, 30)
  
  splsda.stemcells <- splsda(X = X, Y = Y, 
                             ncomp = optimal.ncomp, 
                             keepX = optimal.keepX)
  expect_error(plotArrow(splsda.stemcells), "'plotArrow' not implemented for (s)PLSDA or MINT sPLSDA", fixed = TRUE)
})

test_that("plotArrow functions on DIABLO objects", code = {
  
  data(breast.TCGA)
  
  X <- list(miRNA = breast.TCGA$data.train$mirna, 
            mRNA = breast.TCGA$data.train$mrna,
            proteomics = breast.TCGA$data.train$protein)
  Y <- breast.TCGA$data.train$subtype
  
  optimal.ncomp <- 2
  optimal.keepX <-  list(miRNA = c(10, 5), 
                         mRNA = c(25,16),
                         proteomics = c(8,5))
  tcga.diablo <- block.splsda(X, Y,
                              ncomp = optimal.ncomp,
                              keepX = optimal.keepX)
  
  pA_res <- plotArrow(tcga.diablo)
  
  expect_is(pA_res, "ggplot")
})