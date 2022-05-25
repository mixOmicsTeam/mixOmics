
###############################################################################
### ============================ MISSING TESTS ============================ ###
###############################################################################


###############################################################################
### ============================ GROUND TRUTHS ============================ ###
###############################################################################

Test.Data <- readRDS(system.file("testdata", "testdata-network.rda", package = "mixOmics"))
Testable.Components <- Test.Data$tc
Ground.Truths <- Test.Data$gt

###############################################################################
### ================================ BASIC ================================ ###
###############################################################################


test_that("(basic:network): matrix", {

  testable.components <- Testable.Components$basic.matrix
  GT <- Ground.Truths$basic.matrix

  data(liver.toxicity)
  X <- liver.toxicity$gene[, 1:10]
  Y <- liver.toxicity$clinic[, 1:10]

  res.cor <- cor(X, Y)

  matrix.network <- network(res.cor)

  invisible(capture.output(TT <- dput(matrix.network[testable.components])))

  expect_equal(TT, GT)
})


test_that("(basic:network): pls", {

  testable.components <- Testable.Components$basic.pls
  GT <- Ground.Truths$basic.pls

  data(liver.toxicity)
  X <- liver.toxicity$gene[, 1:10]
  Y <- liver.toxicity$clinic[, 1:10]

  res.pls <- pls(X, Y)
  
  pls.network <- network(res.pls)

  invisible(capture.output(TT <- dput(pls.network[testable.components])))

  expect_equal(TT, GT)
})


test_that("(basic:network): spls", {

  testable.components <- Testable.Components$basic.spls
  GT <- Ground.Truths$basic.spls

  data(liver.toxicity)
  X <- liver.toxicity$gene[, 1:10]
  Y <- liver.toxicity$clinic[, 1:10]

  choice.keepX <- c(3,3)
  choice.keepY <- c(3,3)

  res.spls <- spls(X, Y, keepX = choice.keepX, keepY = choice.keepY)

  spls.network <- network(res.spls)

  invisible(capture.output(TT <- dput(spls.network[testable.components])))

  expect_equal(TT, GT)
})


test_that("(basic:network): plsda", {

  testable.components <- Testable.Components$basic.plsda
  GT <- Ground.Truths$basic.plsda

  data(srbct)
  X <- srbct$gene[, 1:10]
  Y <- srbct$class

  res.plsda <- plsda(X, Y)

  plsda.network <- network(res.plsda)

  invisible(capture.output(TT <- dput(plsda.network[testable.components])))

  expect_equal(TT, GT)
})


test_that("(basic:network): splsda", {

  testable.components <- Testable.Components$basic.splsda
  GT <- Ground.Truths$basic.splsda

  data(srbct)
  X <- srbct$gene[, 1:10]
  Y <- srbct$class

  choice.keepX <- c(3,3)

  res.splsda <- splsda(X, Y, keepX = choice.keepX)

  splsda.network <- network(res.splsda)

  invisible(capture.output(TT <- dput(splsda.network[testable.components])))

  expect_equal(TT, GT)
})


test_that("(basic:network): rcc", {

  testable.components <- Testable.Components$basic.rcc
  GT <- Ground.Truths$basic.rcc

  data(nutrimouse)
  X <- nutrimouse$lipid[, 1:10]
  Y <- nutrimouse$gene[, 1:10]

  res.rcc <- rcc(X, Y, method = "shrinkage")

  rcc.network <- network(res.rcc)

  invisible(capture.output(TT <- dput(rcc.network[testable.components])))

  expect_equal(TT, GT)
})


test_that("(basic:network): sgcca", {

  testable.components <- Testable.Components$basic.sgcca
  GT <- Ground.Truths$basic.sgcca

  data(breast.TCGA)
  X = list(miRNA = breast.TCGA$data.train$mirna[, 1:10],
           mRNA = breast.TCGA$data.train$mrna[, 1:10],
           proteomics = breast.TCGA$data.train$protein[, 1:10])

  choice.keepX <- list(miRNA = c(3,3),
                       mRNA = c(3,3),
                       proteomics = c(3,3))

  res.sgcca <- wrapper.sgcca(X, keepX = choice.keepX, ncomp = 3)

  sgcca.network <- network(res.sgcca)

  invisible(capture.output(TT <- dput(sgcca.network[testable.components])))

  expect_equal(TT, GT)
})


test_that("(basic:network): block.spls", {

  testable.components <- Testable.Components$basic.block.spls
  GT <- Ground.Truths$basic.block.spls

  data(breast.TCGA)
  X = list(miRNA = breast.TCGA$data.train$mirna[, 1:10],
           mRNA = breast.TCGA$data.train$mrna[, 1:10],
           proteomics = breast.TCGA$data.train$protein[, 1:10])

  choice.keepX <- list(miRNA = c(3,3),
                       mRNA = c(3,3),
                       proteomics = c(3,3))

  res.block.spls <- block.spls(X, indY=3, keepX = choice.keepX, ncomp = 3)

  block.spls.network <- network(res.block.spls)

  invisible(capture.output(TT <- dput(block.spls.network[testable.components])))

  expect_equal(TT, GT)
})


test_that("(basic:network): sgccda", {

  testable.components <- Testable.Components$basic.sgccda
  GT <- Ground.Truths$basic.sgccda

  data(breast.TCGA)
  X = list(miRNA = breast.TCGA$data.train$mirna[, 1:10],
           mRNA = breast.TCGA$data.train$mrna[, 1:10],
           proteomics = breast.TCGA$data.train$protein[, 1:10])
  Y = breast.TCGA$data.train$subtype

  choice.keepX <- list(miRNA = c(3,3),
                       mRNA = c(3,3),
                       proteomics = c(3,3))

  res.sgccda <- wrapper.sgccda(X, Y, keepX = choice.keepX)

  sgccda.network <- network(res.sgccda)

  invisible(capture.output(TT <- dput(sgccda.network[testable.components])))

  expect_equal(TT, GT)
})


test_that("(basic:network): block.splsda", {

  testable.components <- Testable.Components$basic.sgccda
  GT <- Ground.Truths$basic.sgccda

  data(breast.TCGA)
  X = list(miRNA = breast.TCGA$data.train$mirna[, 1:10],
           mRNA = breast.TCGA$data.train$mrna[, 1:10],
           proteomics = breast.TCGA$data.train$protein[, 1:10])
  Y = breast.TCGA$data.train$subtype

  choice.keepX <- list(miRNA = c(3,3),
                       mRNA = c(3,3),
                       proteomics = c(3,3))

  res.block.splsda <- block.splsda(X, Y, keepX = choice.keepX)

  block.splsda.network <- network(res.block.splsda)

  invisible(capture.output(TT <- dput(block.splsda.network[testable.components])))

  expect_equal(TT, GT)
})


###############################################################################
### ============================== PARAMETER ============================== ###
###############################################################################


test_that("(parameter:network): cutoff", {

  testable.components <- Testable.Components$cutoff.pls
  GT <- Ground.Truths$cutoff.pls

  data(liver.toxicity)
  X <- liver.toxicity$gene[, 1:10]
  Y <- liver.toxicity$clinic[, 1:10]

  res.pls <- pls(X, Y)

  cutoff.network <- network(res.pls, cutoff = 0.4)

  invisible(capture.output(TT <- dput(cutoff.network[testable.components])))

  expect_equal(TT, GT)
})


test_that("(parameter:network): blocks", {

  testable.components <- Testable.Components$blocks.block.spls
  GT <- Ground.Truths$blocks.block.spls

  data(breast.TCGA)
  X = list(miRNA = breast.TCGA$data.train$mirna[, 1:10],
           mRNA = breast.TCGA$data.train$mrna[, 1:10],
           proteomics = breast.TCGA$data.train$protein[, 1:10])

  choice.keepX <- list(miRNA = c(3,3),
                       mRNA = c(3,3),
                       proteomics = c(3,3))

  res.block.spls <- block.spls(X, indY=3, keepX = choice.keepX, ncomp = 3)

  blocks.network <- network(res.block.spls, blocks = c(1,3))

  invisible(capture.output(TT <- dput(blocks.network[testable.components])))

  expect_equal(TT, GT)
})


###############################################################################
### ================================ ERROR ================================ ###
###############################################################################


test_that("(error:network): catches unaccepted file extensions", {
  
  data(liver.toxicity)
  X <- liver.toxicity$gene[, 1:10]
  Y <- liver.toxicity$clinic[, 1:10]

  choice.keepX <- c(3,3)
  choice.keepY <- c(3,3)

  res.spls <- spls(X, Y, keepX = choice.keepX, keepY = choice.keepY)

  expect_error(network(res.spls, name.save = "network", save = "raw"),
               "'save' must be one of 'jpeg', 'png', 'tiff' or 'pdf'.",
               fixed = T)
})


test_that("(error:network): catches unaccepted input objects", {

  data(stemcells)
  X <- stemcells$gene

  res.pca <- pca(X)

  expect_error(network(res.pca),
               "'network' is only implemented for the following objects: matrix, pls, plsda, spls, splsda, rcc, sgcca, rgcca, sgccda",
               fixed = T)
})


test_that("(error:network): catches unaccepted row.names/col.names values", {
  
  data(liver.toxicity)
  X <- liver.toxicity$gene[, 1:10]
  Y <- liver.toxicity$clinic[, 1:10]

  choice.keepX <- c(3,3)
  choice.keepY <- c(3,3)

  res.spls <- spls(X, Y, keepX = choice.keepX, keepY = choice.keepY)

  expect_error(network(res.spls, row.names = "random.name"),
               "'row.names' must be a character vector of 10 unique entries.",
               fixed = T)

  expect_error(network(res.spls, col.names = "random.name"),
               "'col.names' must be a character vector of 10 unique entries.",
               fixed = T)
})


test_that("(error:network): catches unaccepted blocks values", {
  
  data(breast.TCGA)
  X = list(miRNA = breast.TCGA$data.train$mirna[, 1:10],
           mRNA = breast.TCGA$data.train$mrna[, 1:10],
           proteomics = breast.TCGA$data.train$protein[, 1:10])

  choice.keepX <- list(miRNA = 3,
                       mRNA = 3,
                       proteomics = 3)

  res.block.spls <- block.spls(X, ncomp=1, indY=3, keepX = choice.keepX)

  expect_error(network(res.block.spls, blocks = c("random.block", "random.block")),
               "One element of 'blocks' does not match with the names of the blocks",
               fixed = T)
  
  expect_error(network(res.block.spls, blocks = c(3,4)),
               "Incorrect value for 'blocks",
               fixed = T)
  
})


test_that("(error:network): catches cutoff values which are too high", {
  
  data(liver.toxicity)
  X <- liver.toxicity$gene[, 1:10]
  Y <- liver.toxicity$clinic[, 1:10]

  res.spls <- spls(X, Y)

  expect_error(network(res.spls, cutoff = 0.99),
               "You have chosen a high cutoff value of 0.99 which is greaer than the max value in the similarity matrix which is 0.43",
               fixed=T)
})



dev.off()
unlink(list.files(pattern = "*.pdf"))
