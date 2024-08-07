# ==============================================================================
# Kilmer's t-SVD decomposition for order-3 tensors
# ==============================================================================

#' @description Return the t-SVDM decomposition from Kilmer et al. (2021).
#' @export
tsvdm <- function(
  x,
  m = NULL,
  minv = NULL,
  keep_hats = FALSE,
  full_frontal_slices = TRUE,
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
  # bltodo: bit clumsy repeated code from m_product() function body
  if (
    xor(is.function(m), is.function(minv)) ||
      xor(is.null(m), is.null(minv))
  ) {
    stop(
      "If explicitly defined, both m and its inverse must be defined as 
      functions."
    )
  }
  # use dctii as default transform if user does not specify an explicit one
  if (is.null(m)) {
    transforms <- dctii_m_transforms(t, bpparam = bpparam)
    m <- transforms$m
    minv <- transforms$minv
  }

  xhat <- m(x)
  uhat <- array(0, dim = c(n, n, t))
  vhat <- array(0, dim = c(p, p, t))

  if (svals_matrix_form) {
    shat <- array(0, dim = c(k, t))
    for (i in seq_len(t)) {
      facewise_svd <- svd(xhat[, , i])
      uhat[, , i] <- facewise_svd$u
      shat[, i] <- facewise_svd$d
      vhat[, , i] <- facewise_svd$v
    }
  } else {
    shat <- array(0, dim = c(n, p, t))
    for (i in seq_len(t)) {
      facewise_svd <- svd(xhat[, , i])
      uhat[, , i] <- facewise_svd$u
      shat[1:k, 1:k, i] <- diag(facewise_svd$d)
      vhat[, , i] <- facewise_svd$v
    }
  }

  if (keep_hats) {
    return(list(uhat = uhat, shat = shat, vhat = vhat))
  } else {
    output <- lapply(list(uhat, shat, vhat), minv)
    names(output) <- c("u", "s", "v")
    return(output)
  }
}
