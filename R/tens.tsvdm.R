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
