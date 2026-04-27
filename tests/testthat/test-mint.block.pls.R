context("mint.block.pls")

# Test
test_that("mint.block.pls works", {
  data(breast.TCGA)
  study = c(rep("study1", 150), rep("study2", 70))
  mrna = rbind(breast.TCGA$data.train$mrna, breast.TCGA$data.test$mrna)
  mirna = rbind(breast.TCGA$data.train$mirna, breast.TCGA$data.test$mirna)
  Y = mrna[, 1]
  mrna = mrna[, -1]
  data = list(mrna = mrna, mirna = mirna)
  res = mint.block.plsda(data, Y, study=study, ncomp=2)
  expect_equal(res$loadings$mrna[5], -0.08031856, tolerance = 1e-5)
  expect_equal(res$mode, "regression")
})
  
               