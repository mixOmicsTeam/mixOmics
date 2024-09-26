# ==============================================================================
# Tensor plsda generalization; developed @ Melbourne Integrative Genomics based
# on Kilmer's tensor m-product algebra
# ==============================================================================

#' Check if y input is a vector, matrix or tensor that can be inferred as
#' classes across all time points, if so return these classes in consistent form
#' as a vector.
#'
#' @keywords internal
.y_to_vec <- function(y) {
  if (is.null(dim(y))) {
    return(y)
  }
  if (length(dim(y)) == 2 && dim(y)[2] != 1) {
    return(c(y))
  }
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

.y_to_tens <- function() {
  NULL
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

  if (multilevel) {

  } else {
    y_vec <- .y_to_vec(y)
    if (dim(x)[3] > 1) {
      y_vec <- factor(y_vec)
      if (nlevels(y_vec) == 1) {
        stop("y should contain 2 or more distinct classes")
      }

      y_mat <- unmap(y_vec)
      # BLTODO: WANT colnames(y_mat) = levels(y_vec), and be able to propogate
      # names into tpls - need to work on preserving names for tpls and tpca
      # when mixOmics2 is rolled out to end users

      # if y specifies the same classes across all time points, fill y into a
      # tensor to match the number of time points in t
      y <- array(y_mat, dim = c(dim(y_mat), dim(x)[3]))
    }
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
