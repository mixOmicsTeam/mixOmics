context("plotArrow")

## ------------------------------------------------------------------------ ##

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

## ------------------------------------------------------------------------ ##
## vdiffr testing - "ggplot2"
library(vdiffr)

## spls model
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
spls.obj <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50),
                 keepY = c(10, 10, 10))

test_that("plotArrow works for spls objects", {
  skip_on_ci() # only run the vdiffr tests locally
  
  # simple plot showing sample names
  invisible(capture.output(
    expect_doppelganger(
      title = "Arrow plot spls model", 
      fig = plotArrow(spls.obj))
  ))
  # samples coloured by primary groups, sample names, default colours
  invisible(capture.output(
    expect_doppelganger(
      title = "Arrow plot spls with customised colours",
      fig = plotArrow(spls.obj,  group = as.factor(liver.toxicity$treatment[, 'Time.Group']),
                      col = c("red", "blue", "purple", "darkgreen"), 
                      ind.names  = liver.toxicity$treatment[, 'Dose.Group'],
                      legend = TRUE, position.names = 'start', legend.title = 'Time.Group'))
  ))
})


## diablo model
data(breast.TCGA)
idx = seq(1, length(breast.TCGA$data.train$subtype), 10)
X <- list(mRNA = breast.TCGA$data.train$mrna[idx,],
          miRNA = breast.TCGA$data.train$mirna[idx,],
          protein = breast.TCGA$data.train$protein[idx,])
Y <- breast.TCGA$data.train$subtype[idx] # set the response variable

diablo.obj <- block.splsda(X, Y, ncomp = 2) # undergo multiblock sPLS-DA

test_that("plotArrow works for diablo objects", {
  skip_on_ci() # only run the vdiffr tests locally
  
  # simple plot showing sample names
  invisible(capture.output(
    expect_doppelganger(
      title = "Arrow plot diablo model", 
      fig = plotArrow(diablo.obj, 
                      ind.names = FALSE,
                      legend = TRUE,
                      title = 'TCGA, DIABLO comp 1 - 2'))
  ))
  # samples coloured by primary groups, sample names, default colours
  pchs <- c(3, 2, 1)
  names(pchs) <- c("miRNA", "mRNA", "protein")
  invisible(capture.output(
    expect_doppelganger(
      title = "Arrow plot daiblo with customised colours",
      fig = plotArrow(diablo.obj, 
                      ind.names = FALSE,
                      legend = TRUE,
                      title = 'TCGA, DIABLO comp 1 - 2',
                      pch = pchs) )
  ))
})

