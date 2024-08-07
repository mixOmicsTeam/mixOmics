# ==============================================================================
# Tensor pca based on Mor's TCAM algorithm
# ==============================================================================

#' @description R implementation of np.unravel_index. NOTE: currently only works
#' for 1D to 2D column-major conversion.
#' @keywords internal
.unravel_index <- function(indices, dim) {
  nrows <- dim[1]
  return(
    list(
      row_indxs = sapply(indices, function(x) (x - 1) %% nrows + 1),
      col_indxs = sapply(indices, function(x) (x - 1) %/% nrows + 1)
    )
  )
}

#' @description Helper function to convert compressed matrix form of the
#' singular values into a sparse tensor; saves another costly call to tsvdm
#' when computing the transform.
#' @param mat Matrix s.t. each column contains the f-diagonal singular values
#' @param dim Dimension of output tensor
#' @keywords internal
.singular_vals_mat_to_tens <- function(mat, dim) {
  n <- dim[1]
  p <- dim[2]
  t <- dim[3]
  k <- min(n, p)

  tens <- array(0, dim = dim)
  for (i in seq_len(t)) {
    tens[1:k, 1:k, i] <- diag(mat[, i])
  }
  return(tens)
}

#' @description Tensor analogue of PCA introduced by Mor et al. (2022) based on
#' Kilmer's m-product algebra and tsvdm
#' @export
tpca <- function(
  x,
  ncomp = NULL,
  m = NULL,
  minv = NULL,
  center = TRUE,
  bpparam = NULL
) {
  if (length(dim(x)) != 3) {
    stop("Please ensure input tensor is an order-3 array.")
  } else {
    n <- dim(x)[1]
    p <- dim(x)[2]
    t <- dim(x)[3]
    k <- min(n, p)
  }

  if (center) {
    x <- apply(x, c(1, 3), mean)
  }

  # the idea here is to just pass whatever m and minv input has been provided 
  # in this function call to be handled in the tsvdm call
  # bltodo: adopt this approach in tred
  tsvdm_decomposition <- tsvdm(
    x, m, minv,
    keep_hats = TRUE,
    svals_matrix_form = TRUE,
    bpparam = bpparam
  )

  # flatten out Fortran column major style, then get sort order
  singular_values <- as.vector(tsvdm_decomposition$shat)
  k_t_flatten_sort <- order(singular_values, decreasing = TRUE)
  singular_values <- singular_values[k_t_flatten_sort]

  # get explained variance
  # bltodo: check relevant tensor theory for this!
  squared_singular_values <- singular_values ^ 2
  total_var <- sum(squared_singular_values)
  explained_variance_ratio <- squared_singular_values / total_var

  # process n_components input
  if (is.null(ncomp)) {
    ncomp <- k * t
  } else if (ncomp > 0 && ncomp < 1) {
    # pick ncomp to be minimum number of integer components to explain the
    # inputted variance
    ratio_cumsum <- cumsum(explained_variance_ratio)
    ncomp <- findInterval(ncomp, ratio_cumsum) + 1
  } else if (ncomp %% 1 == 0) {
    stopifnot(ncomp > 0 && ncomp <= (k * t))
  } else {
    stop("Please input an integer, 0 < float < 1, or NULL for ncomp parameter")
  }

  # convert the argsort indexes back into the two dimensional indexes
  # corresponding to the singular values matrix
  # NOTE: these are the collection if i_h's and j_h's in Mor et al. (2022)
  k_t_flatten_sort <- .unravel_index(
    k_t_flatten_sort,
    dim(tsvdm_decomposition$shat)
  )

  # compute the values of the projected data
  # bltodo: benchmark this against new_data *_M vhat
  x_projected <- facewise_product(
    tsvdm_decomposition$uhat,
    .singular_vals_mat_to_tens(tsvdm_decomposition$shat)
  )

  # bltodo: compute rho and return?
  return(invisible(list(
    ncomp = ncomp,
    x = x_projected,
    variates = x_projected,
    explained_variance = explained_variance_ratio[1:ncomp]
  )))
}
