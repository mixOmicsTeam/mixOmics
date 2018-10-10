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
  expect_equal(nutrimouse.sgccda$X$gene, data$gene)
  expect_equal(nutrimouse.sgccda$X$lipid, data$lipid)

  expect_is(nutrimouse.sgccda$ind.mat, "matrix")

  expect_equal(dim(nutrimouse.sgccda$ind.mat), c(40L, 5L))

  expect_equal(nutrimouse.sgccda$Y, Y)

  expect_equal(nutrimouse.sgccda$ncomp, c("gene" = 2L, "lipid" = 2L, "Y" = 2L))

  expect_equal(nutrimouse.sgccda$mode, "regression")

  expect_equal(dim(nutrimouse.sgccda$loadings$gene), c(120L, 2L))
  expect_equal(dim(nutrimouse.sgccda$loadings$lipid), c(21L, 2L))
  expect_equal(dim(nutrimouse.sgccda$loadings$Y), c(5L, 2L))

  expect_null(nutrimouse.sgccda$nzv)

  expect_true(nutrimouse.sgccda$scale)

  expect_equal(nutrimouse.sgccda$scheme, "centroid")
  expect_equal(nutrimouse.sgccda$scheme, "centroid")

  expect_equal(nutrimouse.sgccda$indY, 3L)

  expect_equal(nutrimouse.sgccda$weights,
               c(gene = 0.694506104274723, lipid = 0.915845972615744))

  expect_length(nutrimouse.sgccda$explained_variance, 3L)
  expect_is(nutrimouse.sgccda$explained_variance, "list")
  expect_equal(names(nutrimouse.sgccda$explained_variance), colnames(design))

  expect_length(nutrimouse.sgccda$AVE, 3L)
  expect_equal(names(nutrimouse.sgccda$AVE),
               c("AVE_X", "AVE_outer", "AVE_inner"))
  expect_equal(nutrimouse.sgccda$AVE$AVE_outer[1], 0.217938372815004)
  expect_equal(nutrimouse.sgccda$AVE$AVE_inner[1], 0.663209598406049)
  expect_equal(nutrimouse.sgccda$AVE$AVE_X$Y[1], c(`comp 1` = 0.25))

})
