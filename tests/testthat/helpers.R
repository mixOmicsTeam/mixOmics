#' Expect identical ignoreing $call slot
#'
#' @param object an R object
#' @param expected an R object
#' @param ... other args passed to expect_identical
#'
#' @return Logical, TRUE if all elements expect $call are identical
#'
.almost_identical <- function(object, expected, ...) {
    require(testthat)
    object$call <- expected$call <- NULL
    expect_identical(object, expected, ...)
}
