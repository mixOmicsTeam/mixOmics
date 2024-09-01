# ==============================================================================
# Kilmer's t-SVD decomposition for order-3 tensors
# ==============================================================================

#' @description Return the t-SVDM decomposition from Kilmer et al. (2021).
#' @author Brendan Lu
#' @export
tsvdm <- function(
  x,
  m = NULL,
  minv = NULL,
  transform = TRUE, # differs to tred, control initial transform to m-space
  keep_hats = FALSE,
  full_frontal_slices = TRUE,
  svals_matrix_form = FALSE,
  facewise_truncate = NULL,
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

  if (!is.null(facewise_truncate)) {
    if (is.integer(facewise_truncate)) {
      k <- facewise_truncate
    } else {
      stop("Please input an integer or NULL for facewise_truncate parameter")
    }
  }

  if (transform) {
    .stop_invalid_transform_input(m, minv)

    # use dctii as default transform if user does not specify an explicit one
    if (is.null(m)) {
      transforms <- dctii_m_transforms(t, bpparam = bpparam)
      m <- transforms$m
      minv <- transforms$minv
    }

    x <- m(x)
  }

  if (full_frontal_slices) {
    u <- array(0, dim = c(n, n, t))
    v <- array(0, dim = c(p, p, t))
  } else {
    u <- array(0, dim = c(n, k, t))
    v <- array(0, dim = c(p, k, t))
  }

  if (svals_matrix_form) {
    s <- array(0, dim = c(k, t))
    for (i in seq_len(t)) {
      facewise_svd <- svd(x[, , i])
      u[, 1:k, i] <- facewise_svd$u
      s[, i] <- facewise_svd$d
      v[, 1:k, i] <- facewise_svd$v
    }
  } else {
    s <- array(0, dim = c(n, p, t))
    for (i in seq_len(t)) {
      facewise_svd <- svd(x[, , i])
      u[, 1:k, i] <- facewise_svd$u
      s[1:k, 1:k, i] <- diag(facewise_svd$d)
      v[, 1:k, i] <- facewise_svd$v
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
