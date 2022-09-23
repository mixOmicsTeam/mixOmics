
###############################################################################
### ================================ BASIC ================================ ###
###############################################################################


test_that("(bin.color:basic): defaults", {
  
  data(nutrimouse)
  X <- nutrimouse$lipid[, 1:10]
  Y <- nutrimouse$gene[, 1:10]
  
  mat <- cor(X, Y)
  cols <- colorRampPalette(c("blue", "red"))(100) 
  
  out <- bin.color(mat=mat,
                   cutoff=0,
                   breaks=NULL,
                   col=cols,
                   symkey=TRUE)
  
  .expect_numerically_close(out$breaks[[1]], -0.761)
  .expect_numerically_close(out$breaks[[100]], 0.75)
  expect_equal(out$col[[51]], "#80007E")
})


###############################################################################
### ============================== PARAMETER ============================== ###
###############################################################################


test_that("(bin.color:basic): cutoff=0 & breaks is not NULL", {
  
  data(nutrimouse)
  X <- nutrimouse$lipid[, 1:10]
  Y <- nutrimouse$gene[, 1:10]
  
  mat <- cor(X, Y)
  cols <- colorRampPalette(c("blue", "red"))(100) 
  
  out <- bin.color(mat=mat,
                   cutoff=0,
                   breaks=seq(min(mat), max(mat), length.out=101),
                   col=cols,
                   symkey=TRUE)
  
  .expect_numerically_close(out$breaks[[1]], -0.578)
  .expect_numerically_close(out$breaks[[100]], 0.747)
  expect_equal(out$col[[51]], "#80007E")
})


###############################################################################
### ================================ ERROR ================================ ###
###############################################################################


test_that("(bin.color:error): catches when 'breaks' has an odd number", {
  
  data(nutrimouse)
  X <- nutrimouse$lipid[, 1:10]
  Y <- nutrimouse$gene[, 1:10]
  
  mat <- cor(X, Y)
  cols <- colorRampPalette(c("blue", "red"))(100) 
  expect_error(bin.color(mat=mat,
                         cutoff=0,
                         breaks=1,
                         col=cols,
                         symkey=TRUE),
               "even number")
})


test_that("(bin.color:basic): symkey = FALSE", {
  
  data(nutrimouse)
  X <- nutrimouse$lipid[, 1:10]
  Y <- nutrimouse$gene[, 1:10]
  
  mat <- cor(X, Y)
  cols <- colorRampPalette(c("blue", "red"))(100)
  
  expect_error(bin.color(mat=mat,
                         cutoff=0.3,
                         breaks=c(1),
                         col=cols,
                         symkey=FALSE),
               "one more")
  
  expect_error(bin.color(mat=mat,
                         cutoff=0,
                         breaks=c(1),
                         col=cols,
                         symkey=FALSE),
               "one more")
})
