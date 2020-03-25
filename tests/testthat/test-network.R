context("network")

test_that("network works for rcc", {
  
  ## network representation for objects of class 'rcc'
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
  
  ## create a tmp file
  tmp.file <- tempfile("network", fileext = ".jpeg")
  network.rcc.res <- network(nutri.res, comp = 1:3, cutoff = 0.6, save = "jpeg", 
                             name.save = tmp.file)
  expect_equal(names(network.rcc.res), c("gR", "M", "cutoff"))
  .expect_numerically_close(sum(network.rcc.res$M), 10.8786, digits = 3)
  unlink(tmp.file)
})

test_that("network works for spls", {
  data(liver.toxicity)
  X <- liver.toxicity$gene
  Y <- liver.toxicity$clinic
  toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50),
                        keepY = c(10, 10, 10))
  ## create a tmp file
  tmp.file <- tempfile("network", fileext = ".jpeg")
  network.spls.res <- network(toxicity.spls, comp = 1:3, cutoff = 0.8,
          color.node = c("mistyrose", "lightcyan"),
          shape.node = c("rectangle", "circle"),
          color.edge = color.spectral(100),
          lty.edge = "solid", lwd.edge =  1,
          show.edge.labels = FALSE, interactive = FALSE, save = "jpeg", 
          name.save = tmp.file)
  
  expect_equal(names(network.spls.res), c("gR", "M", "cutoff"))
 .expect_numerically_close(sum(network.spls.res$M), 45.6061, digits = 3)
 unlink(tmp.file)
})

unlink(list.files(pattern = "*.pdf"))
