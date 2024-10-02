# ==============================================================================
# Tensor pca based on Mor's TCAM algorithm
# ==============================================================================

#' R implementation of np.unravel_index. NOTE: currently only works for 1D to 2D
#' column-major conversion, and returns a list of 2D indices. Returns a matrix
#' output of length(indices) columns, with two rows. The first row corresponds
#' to the sorted k indices, and the second row contains the sorted t indices.
#'
#' @author Brendan Lu
#' @keywords internal
.unravel_index <- function(indices, dim) {
  nrows <- dim[1]
  return(sapply(
    indices,
    FUN = function(x) {
      c(
        (x - 1) %% nrows + 1, # transformed k tensor position
        (x - 1) %/% nrows + 1 # transformed t tensor position
      )
    }
  ))
}

#' Extract tensor columns specified by .unravel_index output.
#' Effectively achieves:
#'
#'     self.loadings_matrix_ = hatV[
#'         :, self._k_t_flatten_sort[0], self._k_t_flatten_sort[1]
#'     ].T
#'
#' type of indexing in Numpy.
#'
#' This is the 'tensor compression' based on the ordered indices specified
#' by each column of k_t_indices as described in Mor's TCAM.
#'
#' @author Brendan Lu
#' @keywords internal
.extract_tensor_columns <- function(tensor, k_t_indices) {
  return(apply(
    k_t_indices, 2,
    FUN = function(index_column) tensor[, index_column[1], index_column[2]]
  ))
}

#' Helper function to convert compressed matrix form of the singular values into
#' a sparse tensor.
#'
#' @param mat Matrix s.t. each column contains the f-diagonal singular values
#' @param dim Dimension of output tensor
#' @author Brendan Lu
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

#' Tensor PCA-like dimensionality reduction
#'
#' Tensor analogue of PCA introduced by Mor et al. (2022) based on Kilmer's
#' m-product algebra and tsvdm.
#'
#' @param x Tensor input.
#' @param ncomp The estimated number of components. ncomp can be explicitly set
#' using an integer value, 0 < float value < 1, or left as NULL which will
#' default to the maximum ncomp value possible.
#' @param m A function which applies an orthogonal tensor tubal transform.
#' @param minv The inverse of m.
#' @param center If set to false, the data tensor will not be centralized into
#' Mean Deviation Form (see Mor et al. 2022). By default, the mean horizontal
#' slice of the input tensor(s) are subtracted, so that all of the horizontal
#' slices sum to 0, analgous to centering matrix data.
#' @param matrix_output Collect the top singular vectors across the tensor,
#' organized by the magnitude of the corresponding singular vector, and place
#' them into the columns of a matrix. Corresponds to the 'matrix compression'
#' type of scheme described in Mor et al. 2022.
#' @param bpparam A \linkS4class{BiocParallelParam} object indicating the type
#' of parallelisation. Does not have any effect if transform functions
#' explicitly set using \code{m}, \code{minv}.
#' @author Brendan Lu
#' @export
tpca <- function(
  x,
  ncomp = NULL,
  m = NULL,
  minv = NULL,
  center = TRUE,
  matrix_output = TRUE,
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

  # use dctii as default transform if user does not specify an explicit one
  validated_transforms <- .stop_invalid_transform_input(m, minv, t, bpparam)
  m <- validated_transforms$m
  minv <- validated_transforms$minv

  # bltodo: add scaling as well?
  if (center) {
    mean_slice <- apply(x, c(2, 3), mean)
    x <- sweep(x, c(2, 3), STATS = mean_slice, FUN = "-")
  }

  # NOTE: passing in bpparam configuration is not needed, as any
  # bpparam specification will already take effect above in the transforms
  tsvdm_decomposition <- tsvdm(
    x, m, minv,
    keep_hats = TRUE,
    svals_matrix_form = TRUE
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

  # process ncomp input
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
    k_t_flatten_sort[1:ncomp],
    dim(tsvdm_decomposition$shat)
  )

  # compute the values of the projected data
  # bltodo: benchmark this against \eqn{new_data *_M vhat}
  # bltodo: zero out uhat first? svd truncates it to rank already, but can
  # further truncate to ncomp columns?
  x_projected <- tsvdm_decomposition$uhat %fp%
    .singular_vals_mat_to_tens(tsvdm_decomposition$shat, dim = c(n, p, t))

  loadings <- tsvdm_decomposition$vhat

  # extract the columns in compressed form
  if (matrix_output) {
    x_projected <- .extract_tensor_columns(x_projected, k_t_flatten_sort)
    loadings <- .extract_tensor_columns(
      tsvdm_decomposition$vhat,
      k_t_flatten_sort
    )
  }

  # BLTODO: (to add in) zero out and compute rho as well
  # see _rank_q_truncation_zero_out() in `tred`
  # need to use this for tensor output to be appropriate

  return(invisible(list(
    ncomp = ncomp,
    x = x,
    loadings = loadings,
    variates = x_projected, # bltodo: maybe just adopt a different name here?
    explained_variance = explained_variance_ratio[1:ncomp]
  )))
}
