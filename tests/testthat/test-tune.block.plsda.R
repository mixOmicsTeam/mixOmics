context("tune.block.plsda")
library(BiocParallel)

test_that("tune.block.plsda works and is the same as perf alone and in tune wrapper", code = {
  
  data("breast.TCGA")
  X <- list(mrna = breast.TCGA$data.train$mrna, mirna = breast.TCGA$data.train$mirna, protein = breast.TCGA$data.train$protein)
  Y <- breast.TCGA$data.train$subtype
  set.seed(100)
  subset <- mixOmics:::stratified.subsampling(breast.TCGA$data.train$subtype, folds = 10)[[1]][[1]]
  X <- lapply(X, function(omic) omic[subset,])
  Y <- Y[subset]
  design = matrix(1, ncol = length(X), nrow = length(X),
                  dimnames = list(names(X), names(X)))
  diag(design) =  0
  design
  
  # run tune.block.plsda
  tune.res.1 <- suppressWarnings(
    tune.block.plsda(X, Y, design = design, ncomp = 2, nrepeat = 1, folds = 3,
                               seed = 13, dist = c("all"))
  )
  # run tune wrapper
  tune.res.2 <- suppressWarnings(
    tune(X, Y, design = design, ncomp = 2, nrepeat = 1, folds = 3,
                     seed = 13, dist = c("all"), method = "block.plsda")
  )
  # run perf
  model <- block.plsda(X, Y, design = design, ncomp = 2)
  tune.res.3 <- suppressWarnings(
    perf(model, nrepeat = 1, seed = 13, folds = 3)
  )
  
  # check outputs format
  expect_equal(class(tune.res.1)[1], "perf.sgccda.mthd")
  expect_equal(class(tune.res.2)[1], "perf.sgccda.mthd")
  expect_equal(class(tune.res.3)[1], "perf.sgccda.mthd")
  
  # check output values
  .expect_numerically_close(tune.res.1$error.rate.per.class$mrna$max.dist[3,2], 0.25)
  .expect_numerically_close(tune.res.2$error.rate.per.class$mrna$max.dist[3,2], 0.25)
  .expect_numerically_close(tune.res.3$error.rate.per.class$mrna$max.dist[3,2], 0.25)
  
  # check can plot outputs
  expect_silent(plot(tune.res.1))
  expect_silent(plot(tune.res.2))
  expect_silent(plot(tune.res.3))
  
})
  
  