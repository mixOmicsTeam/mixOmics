context(".get_xy")

test_that(".get_xy works", {
  expect_error(object = .get_xy(mc = list(data = mae_data, X = NULL, Y = Ya)), class = "inv_XY")
  expect_error(object = .get_xy(mc = list(data = mae_data, X = Xa, Y = "foo")), class = "inv_XY")
  expect_true(is(.get_xy(mc = list(data = mae_data, X = Xa, Y = Ya))$Y, "matrix"))
  expect_true(is(.get_xy(mc = list(data = mae_data, X = Xa, Y = Ya))$X, "matrix"))
})


context(".plsMethodsHelper")

test_that(".plsMethodsHelper works", {
  expect_true(is(.plsMethodsHelper(mc = list(data=mae_data, formula=as.formula(paste0(Xa,'~', Ya))))$X, "matrix"))
  expect_error(object = .plsMethodsHelper(mc = list(data=mae_data, formula=as.formula(paste0(Xa,'~', "foo")))), class = "inv_XY")
})