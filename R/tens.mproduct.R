# ==============================================================================
# utilities and functions for Kilmer's m product
# ==============================================================================

#' @description Apply a function across the last dimension of an input vector,
#' matrix, or tensor. This function defines both a parallel algorithm using
#' BiocParrallel, and also a simple \code{apply} algorithm. Even on Windows
#' Machines, setting \code{bpparam = BiocParallel::SerialParam()} offers a
#' notable speedup for larger 3D array inputs.
#' @param x Numerical array input.
#' @param mat Function which defines the tubal transform.
#' @param bpparam A \linkS4class{BiocParallelParam} object indicating the type
#' of parallelisation.
#' @return A tensor of the same size under the specified tubal transform,
#' denoted \hat{x}.
#' @author Saritha Kodikara, Brendan Lu
#' @keywords internal
.apply_mat_transform <- function(x, mat, bpparam) {
  if (length(dim(x)) == 1) {
    return(mat %*% x)
  } else if (length(dim(x)) == 2) {
    return(mat %*% x)
  } else if (length(dim(x)) == 3) {
    n <- dim(x)[1]
    p <- dim(x)[2]
    t <- dim(x)[3]
    if (is.null(bpparam)) {
      # apply algorithm --------------------------------------------------------
      transformed_slices <- apply(
        aperm(x, c(3, 1, 2)), c(2, 3), # permuted tensor with mode 3 in front
        function(slice) mat %*% slice # left multiply by transform mat
      )
      # permute back to original orientation
      return(aperm(array(transformed_slices, dim = c(t, n, p)), c(2, 3, 1)))
      #-------------------------------------------------------------------------
    } else {
      # BiocParallel algorithm -------------------------------------------------
      transformed_slices <- BiocParallel::bplapply(
        lapply(seq_len(p), function(i) t(x[, i, ])), # a list of t x n matrices
        function(slice) mat %*% slice, # left multiply by transform mat
        BPPARAM = bpparam
      )
      # unlist and cast into array with p facewise t x n matrices
      # then appropriately rotate to get original n x p x t matrix
      return(
        aperm(array(unlist(transformed_slices), dim = c(t, n, p)), c(2, 3, 1))
      )
      #-------------------------------------------------------------------------
    }
  } else {
    # error: some array >3D has been inputted
    stop(
      "Only order 1 (vector), 2 (matrix), 3 (3D tensor) arrays are supported"
    )
  }
}

#' @description Returns functions \code{m} and \code{m_inv} which apply tubal
#' transforms defined by the matrix \code{m_mat}.
#' @param mat Function which defines the tubal transform.
#' @param m_inv_mat Function which defines inverse tubal transform
#' @param bpparam A \linkS4class{BiocParallelParam} object indicating the type
#' of parallelisation.
#' @return
#' \item{m}{A function which applies the matrix m_mat along the last dimension
#' of a given numerical input array. For 3D tensors it performs the mode-3
#' product.}
#' \item{minv}{The inverse of m.}
#' @author Brendan Lu
#' @export
matrix_to_m_transforms <- function(
  m_mat,
  m_inv_mat = NULL,
  bpparam = NULL
) {
  # error: non-matrix input
  stopifnot(length(dim(m_mat)) == 2)
  # error: non-square matrix input
  stopifnot(dim(m_mat)[1] == dim(m_mat)[2])

  # invert m_mat if m_inv_mat not specified
  if (is.null(m_inv_mat)) {
    # error: non-invertible input
    m_inv_mat <- solve(m_mat)
  } else {
    # error: specified matrices are not the same size
    stopifnot(identical(dim(m_mat), dim(m_inv_mat)))
  }

  return(list(
    m = function(x) .apply_mat_transform(x, m_mat, bpparam = bpparam),
    minv = function(x) .apply_mat_transform(x, m_inv_mat, bpparam = bpparam)
  ))
}

#' @description Returns functions \code{m} and \code{m_inv} which apply tubal
#' transforms defined by the Discrete Cosine Transform (DCT-II variant). This
#' is equivalent to Scipys DCTI-ii algorithm with \code{norm='ortho'}.
#' @param t The length of the transform.
#' @param bpparam A \linkS4class{BiocParallelParam} object indicating the type
#' of parallelisation.
#' @return
#' \item{m}{A function which applies the dct-ii along the last dimension of a
#' given numerical input array. For 3D tensors it performs the mode-3 product
#' with the DCT matrix.}
#' \item{minv}{The inverse of m}
#' @export
dctii_m_transforms <- function(t, bpparam = NULL) {
  return(matrix_to_m_transforms(m_mat = gsignal::dctmtx(t), bpparam = bpparam))
}

#' @description Compute Kilmer's facewise product. Note that the for-loop
#' implementation is relatively fast, and very readable. There's also a
#' BiocParralel implementation here, but it lacks significant benchmarking
#' results.
#' @keywords internal
.binary_facewise <- function(a, b, bpparam) {
  na <- dim(a)[1]
  pa <- dim(a)[2]
  ta <- dim(a)[3]

  nb <- dim(b)[1]
  pb <- dim(b)[2]
  tb <- dim(b)[3]

  # error: different t for each input
  stopifnot(ta == tb)
  t <- ta
  # error: non-conforming facewise dimensions
  stopifnot(pa == nb)

  if (is.null(bpparam)) {
    # for-loop algorithm -------------------------------------------------------
    fp_ab <- array(0, dim = c(na, pb, t))
    for (i in 1:t) {
      fp_ab[, , i] <- a[, , i] %*% b[, , i]
    }
    return(fp_ab)
    #---------------------------------------------------------------------------
  } else {
    # BiocParallel algorithm ---------------------------------------------------
    return(
      simplify2array(
        BiocParallel::bplapply(
          array(1:t),
          FUN = function(i) a[, , i] %*% b[, , i],
          BPPARAM = bpparam
        )
      )
    )
    #---------------------------------------------------------------------------
  }
}

#' @description Compute Kilmer's facewise product cumulatively across any
#' arbitrary number of tensor inputs.
#' @param ... Arbitrary number of numerical tensor inputs.
#' @param bpparam A \linkS4class{BiocParallelParam} object indicating the type
#' of parallelisation.
#' @return Cumulative facewise product.
#' @export
facewise_product <- function(..., bpparam = NULL) {
  return(
    Reduce(
      function(a, b) .binary_facewise(a, b, bpparam = bpparam),
      list(...)
    )
  )
}
