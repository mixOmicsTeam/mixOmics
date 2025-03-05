context("mint.block.spls")

# Test
test_that("mint.block.spls works", {
  data(breast.TCGA)
  study = c(rep("study1", 150), rep("study2", 70))
  mrna = rbind(breast.TCGA$data.train$mrna, breast.TCGA$data.test$mrna)
  mirna = rbind(breast.TCGA$data.train$mirna, breast.TCGA$data.test$mirna)
  Y = mrna[, 1]
  mrna = mrna[, -1]
  data = list(mrna = mrna, mirna = mirna)
  res = mint.block.splsda(data, Y, study=study, ncomp=2,
                          keepX = list(mrna=c(10,10), mirna=c(20,20)))
  expect_equal(res$loadings$mrna[5], 0)
  expect_equal(res$mode, "regression")
})
  
               