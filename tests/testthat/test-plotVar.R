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
