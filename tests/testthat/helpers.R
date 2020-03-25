#' Expect identical ignoreing $call slot
#'
#' @param object an R object
#' @param expected an R object
#' @param ... other args passed to expect_identical
#' @keywords internal
#' @return Logical, TRUE if all elements expect $call are identical
.almost_identical <- function(object, expected, ...) {
    require(testthat)
    object$call <- expected$call <- NULL
    expect_identical(object, expected, ...)
}


#' Expect equal numerics after rounding
#'
#' @param numeric_value numeric, outcome of test
#' @param expected numeric, reference to check
#' @param digits integer, decimals to round
#' @return logical, if TRUE, values are almost equal
#' @keywords internal
#' @examples
#' .numerically_close(3.414, 3.14, digits =2)
.expect_numerically_close <- function(numeric_value, expected, digits = 2) {
    require(testthat)
    expect_equal(round(numeric_value, digits = digits), round(expected, digits = digits))
}
