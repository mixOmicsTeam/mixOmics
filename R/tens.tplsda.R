# ==============================================================================
# Tensor plsda generalization; developed @ Melbourne Integrative Genomics based
# on Kilmer's tensor m-product algebra
# ==============================================================================

#' Check if y input is a vector, matrix or tensor that can be inferred as fixed
#' class labels for each sample across all time points; if so return these
#' classes in consistent form as a vector.
#'
#' @author Brendan Lu
#' @keywords internal
.y_to_vec <- function(y) {
  # input is already a vector
  if (is.null(dim(y))) {
    return(y)
  }
  # input is a matrix with a single column
  if (length(dim(y)) == 2 && dim(y)[2] != 1) {
    return(c(y))
  }
  # input is a tensor with a single column
  if (length(dim(y)) == 3 && dim(y)[2] == 1 && dim(y)[3] != 1) {
    return(c(y))
  }
  stop("
    Please ensure y input is either a vector, or a matrix / tensor with one
    vertical column representing outcome categories for each sample

    If specifying categories across time points, please specify
    `multilevel = TRUE` in function call and pass in a matrix or tensor with
    the number of columns corresponding to the number of time points
  ")
}

#' Check if y input is a matrix or tensor that can be inferred as repeated class
#' measurements for each sample; if so, return these classes in a consistent
#' form as a n x 1 x t tensor.
#'
#' Note: We do not have to check that the dimensions are compatible with x
#' tensor input, this check will be done in the call to tpls.
#'
#' @author Brendan Lu
#' @keywords internal
.y_to_tens <- function(y) {
  # input is a n x t matrix
  if (length(dim(y)) == 2) {
    n <- dim(y)[1]
    t <- dim(y)[2]
    return(array(y, dim = c(n, 1, t)))
  }
  # input is a n x 1 x t tensor
  if (length(dim(y)) == 3 && dim(y)[2] == 1) {
    return(y)
  }
  stop("
    When `multilevel = TRUE` in tplsda, please ensure that y input is either a
    n x t matrix, or n x 1 x t column tensor, representing the classes for
    each sample across the t repeated measurements.
  ")
}

#' Run tensor PLSDA-like analysis
#'
#' Developed @ Melbourne Integrative Genomics
#'
#' @author Brendan Lu
#' @export
tplsda <- function(
  x,
  y,
  multilevel = FALSE,
  ncomp = NULL,
  m = NULL,
  minv = NULL,
  center = TRUE,
  bpparam = NULL
) {
  # most parameter checking is done in tpls call, only the checks that are
  # unique / necessary for tplsda are done here
  if (length(dim(x)) != 3) {
    stop("Please ensure x input tensor is an order-3 array")
  }

  # BLTODO: may need to vendor unmap.R into big copy-paste script for Saritha
  # or library(mixOmics) at the top would work neatly too
  if (multilevel) {
    y_tens <- .y_to_tens(y)
    n <- dim(y_tens)[1]
    t <- dim(y_tens)[3]

    # we need to determine all the unique levels that show up across all
    # repeated observations for all samples
    factors_list <- lapply(array(seq_len(t)), function(i) factor(y_tens[, , i]))
    unique_levels <- unique(unlist(lapply(factors_list, levels)))

    y <- array(0, dim = c(n, length(unique_levels), t))
    for (i in seq_len(t)) {
      # unmap conveniently allows us to manually control the total sets of
      # groups (class measurements) our input sample is drawn from
      y[, , i] <- unmap(factors_list[[i]], unique_levels)
    }
  } else {
    y_vec <- .y_to_vec(y)
    y_vec <- factor(y_vec)
    if (nlevels(y_vec) == 1) {
      stop("y should contain 2 or more distinct classes")
    }
    y_mat <- unmap(y_vec)
    # bltodo: WANT colnames(y_mat) = levels(y_vec), and be able to propogate
    # names into tpls - need to work on preserving names for tpls and tpca
    # when mixOmics2 is rolled out to end users (see existing plsda)

    # fill y matrix into a tensor with the same t as x tensor input
    y <- array(y_mat, dim = c(dim(y_mat), dim(x)[3]))
  }

  # note the rest of parameter validation will be done in tpls call
  return(invisible(tpls(
    x,
    y,
    ncomp = ncomp,
    m = m,
    minv = minv,
    mode = "regression",
    center = center,
    bpparam = bpparam
  )))
}
