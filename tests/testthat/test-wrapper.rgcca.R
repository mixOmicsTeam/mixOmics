context("wrapper.rgcca")

# Test 1: Check if wrapper.rgcca runs without errors
test_that("wrapper.rgcca runs without error", {
  data(nutrimouse)
  Y = unmap(nutrimouse$diet)
  data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)
  design = matrix(c(0, 1, 1,
                    1, 0, 1,
                    1, 1, 0), ncol = 3, nrow = 3, byrow = TRUE)
  res = wrapper.rgcca(X = data, design = design, tau = c(1, 1, 0), ncomp = 2)
  expect_s3_class(res, "sparse.rgcca")
})

# Test 2: Check if the number of components in the result is correct
test_that("wrapper.rgcca handles ncomp correctly", {
  data(nutrimouse)
  Y = unmap(nutrimouse$diet)
  data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)
  design = matrix(c(0, 1, 1,
                    1, 0, 1,
                    1, 1, 0), ncol = 3, nrow = 3, byrow = TRUE)
  res = wrapper.rgcca(X = data, design = design, tau = c(1, 1, 0), ncomp = 2)
  expect_equal(length(res$ncomp), 3) # 2 ncomp and also Y
})

# Test 3: Check if the 'X' in the result matches the input
test_that("wrapper.rgcca returns correct X data", {
  data(nutrimouse)
  Y = unmap(nutrimouse$diet)
  data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)
  design = matrix(c(0, 1, 1,
                    1, 0, 1,
                    1, 1, 0), ncol = 3, nrow = 3, byrow = TRUE)
  res = wrapper.rgcca(X = data, design = design, tau = c(1, 1, 0), ncomp = 2)
  expect_equal(res$loadings$Y[1,1], -0.74960222)
})

# Test 4: Check if AVE (Average Variance Explained) is calculated
test_that("wrapper.rgcca returns AVE", {
  data(nutrimouse)
  Y = unmap(nutrimouse$diet)
  data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)
  design = matrix(c(0, 1, 1,
                    1, 0, 1,
                    1, 1, 0), ncol = 3, nrow = 3, byrow = TRUE)
  res = wrapper.rgcca(X = data, design = design, tau = c(1, 1, 0), ncomp = 2)
  expect_true("AVE" %in% names(res))
  expect_equal(length(res$AVE), 3)  # One AVE for each block
})

# Test 5: Check if all outputs are returned when 'all.outputs' is TRUE
test_that("wrapper.rgcca returns all outputs", {
  data(nutrimouse)
  Y = unmap(nutrimouse$diet)
  data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)
  design = matrix(c(0, 1, 1,
                    1, 0, 1,
                    1, 1, 0), ncol = 3, nrow = 3, byrow = TRUE)
  res = wrapper.rgcca(X = data, design = design, tau = c(1, 1, 0), ncomp = 2, all.outputs = TRUE)
  expect_true("variates" %in% names(res))
  expect_true("loadings" %in% names(res))
  expect_true("loadings.star" %in% names(res))
  expect_true("crit" %in% names(res))
  expect_true("names" %in% names(res))
})

# Test 6: Check if 'near.zero.var' impacts the results
test_that("wrapper.rgcca respects near.zero.var", {
  data(nutrimouse)
  Y = unmap(nutrimouse$diet)
  data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)
  design = matrix(c(0, 1, 1,
                    1, 0, 1,
                    1, 1, 0), ncol = 3, nrow = 3, byrow = TRUE)
  res = wrapper.rgcca(X = data, design = design, tau = c(1, 1, 0), ncomp = 2, near.zero.var = TRUE)
  expect_s3_class(res, "sparse.rgcca")
})
