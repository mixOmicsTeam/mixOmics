context("plotVar")
## ------------------------------------------------------------------------ ##
test_that("plotVar works for pls with var.names", {
  data(nutrimouse)
  x <- nutrimouse$gene
  y <- nutrimouse$lipid
  ## custom var.names
  var.names <- list(x = seq_along(x), y = seq_along(y))
  
  pls.res <- pls(x, y) 
  df <- plotVar(pls.res , var.names = var.names)
  
  var.names.char.vec <- unname(unlist(lapply(var.names, as.character)))
  
  expect_true(all(df$names == var.names.char.vec))
})
