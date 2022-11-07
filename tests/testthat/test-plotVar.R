context("plotVar")
## ------------------------------------------------------------------------ ##
test_that("plotVar works for pls with var.names", {
  data(nutrimouse)
  x <- nutrimouse$gene
  y <- nutrimouse$lipid
  ## custom var.names
  var.names <- list(x = seq_along(x), y = seq_along(y))
  
  pls.res <- pls(x, y) 
  df <- plotVar(pls.res , var.names = var.names, plot = FALSE)
  
  var.names.char.vec <- unname(unlist(lapply(var.names, as.character)))
  
  expect_true(all(df$names == var.names.char.vec))
  
  ## ------------- spls
  spls.res <- spls(x, y , keepX = c(10, 10))
  df <- plotVar(spls.res , var.names = var.names, plot = FALSE)
  expect_true(is(df, 'data.frame'))
  expect_true(all(df$names %in% as.character(unlist(var.names))))
  
  ## ------------- spca
  var.names = list(seq_along(x))
  spca.res <- spca(x, keepX = c(10, 10))
  df <- plotVar(spca.res, var.names = var.names, plot = FALSE)
  expect_true(all(df$names %in% as.character(unlist(var.names))))
  
})

test_that("plotVar works in block.(s)PLS1 cases", {
  
  data(breast.TCGA)
  X <- list(miRNA = breast.TCGA$data.train$mirna[,1:10],
            mRNA = breast.TCGA$data.train$mrna[,1:10])
  
  Y <- matrix(breast.TCGA$data.train$protein[,4], ncol=1)
  rownames(Y) <- rownames(X$miRNA)
  colnames(Y) <- "response"
  
  block.pls.result <- block.spls(X, Y, design = "full",
                                 keepX = list(miRNA=c(3,3),
                                              mRNA=c(3,3)))
  
  plotVar.result <- plotVar(block.pls.result, plot = FALSE)
  
  expect_equal(as.character(unique(plotVar.result$Block)), c("miRNA", "mRNA"))
})