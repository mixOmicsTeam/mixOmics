# ==============================================================================
# Kilmer's t-SVD decomposition for order-3 tensors
# ==============================================================================

#' Tensor SVD-like decomposition algorithm
#'
#' Return the t-SVDM decomposition from Kilmer et al. (2021).
#'
#' @param x Tensor input.
#' @param m A function which applies an orthogonal tensor tubal transform.
#' @param minv The inverse of m.
#' @param transform Set to FALSE if the input is already in a m-transformed
#' space, and you just want a facewise svd returned (this is used in tpls).
#' @param keep_hats Set to TRUE if you do not want minv to be applied to the
#' output (this is used in tpca), i.e. return outputs in the hat-space.
#' @param svals_matrix_form Setting to TRUE will return a compressed version of
#' the singular values tensor S, where the singular values of each f-diagonal
#' frontal slice becomes the column of a matrix, with t columns in total.
#' @param bpparam A \linkS4class{BiocParallelParam} object indicating the type
#' of parallelisation. Does not have any effect if transform functions
#' explicitly set using \code{m}, \code{minv}.
#' @author Brendan Lu
#' @export
tsvdm <- function(
  x,
  m = NULL,
  minv = NULL,
  transform = TRUE, # differs to tred, control initial transform to m-space
  keep_hats = FALSE,
  svals_matrix_form = FALSE,
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

  if (transform) {
    # use dctii as default transform if user does not specify an explicit one
    validated_transforms <- .stop_invalid_transform_input(m, minv, t, bpparam)
    m <- validated_transforms$m
    minv <- validated_transforms$minv

    x <- m(x)
  }

  u <- array(0, dim = c(n, k, t))
  v <- array(0, dim = c(p, k, t))

  # bltodo: less ugly way to write this?
  if (svals_matrix_form) {
    s <- array(0, dim = c(k, t))
    for (i in seq_len(t)) {
      facewise_svd <- svd(x[, , i])
      u[, , i] <- facewise_svd$u
      s[, i] <- facewise_svd$d
      v[, , i] <- facewise_svd$v
    }
  } else {
    s <- array(0, dim = c(k, k, t))
    for (i in seq_len(t)) {
      facewise_svd <- svd(x[, , i])
      u[, , i] <- facewise_svd$u
      s[, , i] <- diag(facewise_svd$d)
      v[, , i] <- facewise_svd$v
    }
  }

  if (transform && keep_hats) {
    # make clear returning in hat space
    return(list(uhat = u, shat = s, vhat = v))
  } else {
    if (transform) {
      # minv will work on s regardless of what form it is in
      output <- lapply(list(u, s, v), minv)
      names(output) <- c("u", "s", "v")
      return(output)
    } else {
      output <- list(u = u, s = s, v = v)
      return(output)
    }
  }
}
