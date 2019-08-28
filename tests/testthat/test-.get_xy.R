context(".get_xy")
## TODO more descriptive and segmented tests for all functions
test_that(".get_xy works", {
  expect_error(object = .get_xy(mc = list(data = mae_data, X = NULL, Y = Ya)), class = "inv_XY")
  expect_error(object = .get_xy(mc = list(data = mae_data, X = Xa, Y = "foo")), class = "inv_XY")
  
  expect_true(is(.get_xy(mc = list(data = mae_data, X = Xa, Y = Ya))$Y, "matrix"))
  expect_true(is(.get_xy(mc = list(data = mae_data, X = Xa, Y = Ya))$X, "matrix"))
  expect_error(object = .get_xy(mc = list(data = mae_data, X=Xa, Y=Ya), DA=TRUE), class = "inv_XY")
})