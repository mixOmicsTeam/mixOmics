context("test-diabolo")

test_that("block.splsda works", {
  data(nutrimouse)
  Y = nutrimouse$diet

  data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
  design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE,
                  dimnames = list(c("gene", "lipid", "Y"),
                                  c("gene", "lipid", "Y")))


  nutrimouse.sgccda <- block.splsda(X = data,
                                    Y = Y,
                                    design = design,
                                    keepX = list(gene = c(10,10),
                                                 lipid = c(15,15)),
                                    ncomp = 2,
                                    scheme = "centroid",
                                    tol = 1e-30)
  expect_length(nutrimouse.sgccda, 24L)
  expect_equal(names(nutrimouse.sgccda),
               c("call", "X", "Y", "ind.mat", "ncomp", "mode", "keepX", "keepY",
                 "variates", "loadings", "crit", "AVE", "names", "init", "tol",
                 "iter", "max.iter", "nzv", "scale", "design", "scheme", "indY",
                 "weights", "explained_variance"))
  expect_is(nutrimouse.sgccda$X, "list")
  expect_is(nutrimouse.sgccda$design, "matrix")
  expect_equal(nutrimouse.sgccda$design, design)
})
