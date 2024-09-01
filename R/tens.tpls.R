# ==============================================================================
# Tensor pls generalization; developed @ Melbourne Integrative Genomics based
# on Kilmer's tensor m-product algebra
# ==============================================================================

#' @description Convert a singular values tensor (in compressed matrix form) to
#' a set of indices corresponding to the (column,face) pairs of the top `ncomp`
#' singular values. NEEDS singular values to be in matrix form.
#' @author Brendan Lu
#' @keywords internal
.obtain_k_t_flatten_sort <- function(s_mat, ncomp) {
  # bltodo: if ultimately only use once just place it in function body directly
  return(
    .unravel_index(
      order(as.vector(s_mat))[1:ncomp],
      dim(s_mat)
    )
  )
}

#' @author Brendan Lu
#' @export
tpls <- function(
  x,
  y,
  ncomp = NULL,
  m = NULL,
  minv = NULL,
  mode = "canonical",
  center = TRUE,
  matrix_output = TRUE, # FALSE only takes effect for tsvdm method?
  bpparam = NULL
) {
  # allowed modes mirror sklearn's 2D PLS models
  # scikit-learn.org/stable/modules/cross_decomposition.html#cross-decomposition
  allowed_modes <- c("canonical", "regression", "tsvdm")
  if (!(mode %in% allowed_modes)) {
    stop(paste(
      "Please ensure mode is one of: ",
      paste(allowed_modes, collapse = ", ")
    ))
  }

  if (length(dim(x)) != 3) {
    stop("Please ensure x input tensor is an order-3 array")
  } else {
    n <- dim(x)[1]
    p <- dim(x)[2]
    t <- dim(x)[3]
  }

  if (length(dim(y)) != 3) {
    stop("Please ensure y input tensor is an order-3 array")
  } else {
    n2 <- dim(x)[1]
    q <- dim(x)[2]
    t2 <- dim(x)[3]
  }

  if (n != n2) {
    stop("Please ensure x and y tensor inputs have matching number of samples")
  }

  if (t != t2) {
    stop("Please ensure x and y tensor inputs have matching number of time 
    points")
  }

  k <- min(n, p, q)
  # maximum number of non-zero entries in the f-diagonal singular values tensor
  # of tsvdm(XtY)
  max_rank <- k * t

  .stop_invalid_transform_input(m, minv)

  # use dctii as default transform if user does not specify an explicit one
  if (is.null(m)) {
    transforms <- dctii_m_transforms(t, bpparam = bpparam)
    m <- transforms$m
    minv <- transforms$minv
  }

  if (center) {
    mean_slice_x <- apply(x, c(2, 3), mean)
    mean_slice_y <- apply(y, c(2, 3), mean)
    x <- sweep(x, c(2, 3), STATS = mean_slice_x, FUN = "-")
    y <- sweep(y, c(2, 3), STATS = mean_slice_y, FUN = "-")
  }

  # compute tsvdm of XtY based on tensor facewise transpose and facewise prod.
  tsvdm_decomposition_xty <- tsvdm(
    ft(m(x)) %fp% m(y),
    transform = FALSE,
    full_frontal_slices = FALSE,
    svals_matrix_form = TRUE,
    facewise_truncate = k
  )

  # process ncomp input, much simpler than tpca as we only accept integer input
  # bltodo: investigate non integer input options? explained variance semantics?
  if (is.null(ncomp)) {
    ncomp <- max_rank
  } else if (ncomp %% 1 == 0) {
    stopifnot(ncomp > 0 && ncomp <= max_rank)
  } else {
    stop("Please input an integer or NULL for ncomp parameter")
  }

  # bltodo: tpls explained variance?
  # abstracted away into function as we do not need to concern ourselves with
  # any notion of explained variance
  k_t_flatten_sort <- .obtain_k_t_flatten_sort(tsvdm_decomposition_xty$s, ncomp)

  if (mode == "tsvdm") {
    # simplest algorithm - just uses everything from the tsvdm call with no
    # deflation
    x_loadings <- tsvdm_decomposition_xty$u
    y_loadings <- tsvdm_decomposition_xty$v

    # project the scaled / transformed x, y data onto the loadings
    x_projected <- x %fp% x_loadings
    y_projected <- y %fp% y_loadings

    if (matrix_output) {
      x_loadings <- .extract_tensor_columns(x_loadings, k_t_flatten_sort)
      y_loadings <- .extract_tensor_columns(y_loadings, k_t_flatten_sort)
      x_projected <- .extract_tensor_columns(x_projected, k_t_flatten_sort)
      y_projected <- .extract_tensor_columns(y_projected, k_t_flatten_sort)
    }

    return(invisible(list(
      ncomp = ncomp,
      x = x,
      y = y,
      x_loadings = x_loadings,
      y_loadings = y_loadings,
      x_projected = x_projected,
      y_projected = y_projected
    )))
  }
}
