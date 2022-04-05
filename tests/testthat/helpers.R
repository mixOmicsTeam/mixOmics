#' Expect identical ignoring $call slot
#'
#' @param object an R object
#' @param expected an R object
#' @param ... other args passed to expect_identical
#' @keywords internal
#' @return Logical, TRUE if all elements expect $call are identical
.almost_identical <- function(object, expected, ...) {
    object$call <- expected$call <- NULL
    expect_identical(object, expected, ...)
}


#' Expect equal numeric after rounding
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


.load_url <- function (url, ..., sha1 = NULL) {
  # based very closely on code for devtools::source_url
  stopifnot(is.character(url), length(url) == 1)
  temp_file <- tempfile()
  on.exit(unlink(temp_file))
  request <- httr::GET(url)
  httr::stop_for_status(request)
  writeBin(httr::content(request, type = "raw"), temp_file)
  file_sha1 <- digest::digest(file = temp_file, algo = "sha1")
  if (is.null(sha1)) {
    message("SHA-1 hash of file is ", file_sha1)
  }
  else {
    if (nchar(sha1) < 6) {
      stop("Supplied SHA-1 hash is too short (must be at least 6 characters)")
    }
    file_sha1 <- substr(file_sha1, 1, nchar(sha1))
    if (!identical(file_sha1, sha1)) {
      stop("SHA-1 hash of downloaded file (", file_sha1, 
           ")\n  does not match expected value (", sha1, 
           ")", call. = FALSE)
    }
  }
  load(temp_file, envir = .GlobalEnv)
}
