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

#' @description Get the k, t index corresponding to the largest singular value.
#' Mirrors .unravel_index() but more efficient as we only care about the largest
#' singular value in tpls.
#' @author Brendan Lu
#' @keywords internal
.obtain_k_t_top <- function(s_mat) {
  flat_index <- which.max(as.vector(s_mat))
  nrows <- dim(s_mat)[1]
  return(c(
    (flat_index - 1) %% nrows + 1, # transformed k tensor position
    (flat_index - 1) %/% nrows + 1 # transformed t tensor position
  ))
}

#' @author Brendan Lu
#' @export
tpls <- function(
  x,
  y,
  ncomp = NULL,
  m = NULL,
  minv = NULL,
  mode = "regression",
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
    n2 <- dim(y)[1]
    q <- dim(y)[2]
    t2 <- dim(y)[3]
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

  # project x and y into 'hat-space', expensive computation so do it once and
  # save into variable
  xhat <- m(x)
  yhat <- m(y)

  # process ncomp input, much simpler than tpca as we only accept integer input
  # bltodo: investigate non integer input options? explained variance semantics?
  if (is.null(ncomp)) {
    ncomp <- max_rank
  } else if (ncomp %% 1 == 0) {
    stopifnot(ncomp > 0 && ncomp <= max_rank)
  } else {
    stop("Please input an integer or NULL for ncomp parameter")
  }

  if (mode == "tsvdm") {
    # simplest algorithm - just uses everything from the tsvdm call of XtY based
    # without any deflation steps
    tsvdm_decomposition_xty <- tsvdm(
      ft(xhat) %fp% yhat,
      transform = FALSE,
      svals_matrix_form = TRUE
    )

    # loadings are just the u,v tensors from the tsvdm call
    x_loadings <- tsvdm_decomposition_xty$u
    y_loadings <- tsvdm_decomposition_xty$v

    # project the scaled / transformed x, y data onto the loadings
    x_projected <- xhat %fp% x_loadings
    y_projected <- yhat %fp% y_loadings

    if (matrix_output) {
      # bltodo: tpls explained variance?
      k_t_flatten_sort <- .obtain_k_t_flatten_sort(
        tsvdm_decomposition_xty$s,
        ncomp
      )
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

  } else if (mode == "canonical" || mode == "regression") {
    # preallocate output arrays which we will fill during the iterative process
    x_loadings <- array(0, dim = c(p, ncomp))
    y_loadings <- array(0, dim = c(q, ncomp))
    x_projected <- array(0, dim = c(n, ncomp))
    y_projected <- array(0, dim = c(n, ncomp))

    for (i in seq_len(ncomp)) {
      # compute tsvdm
      # BLTODO: tensor crossprod to speed up?
      tsvdm_decomposition_xty <- tsvdm(
        ft(xhat) %fp% yhat,
        transform = FALSE,
        svals_matrix_form = TRUE
      )

      # get indices corresponding to largest singular value
      k_t_top <- .obtain_k_t_top(tsvdm_decomposition_xty$s)
      curr_x_loadings <- tsvdm_decomposition_xty$u[, k_t_top[1], k_t_top[2]]
      curr_y_loadings <- tsvdm_decomposition_xty$v[, k_t_top[1], k_t_top[2]]

      # note that only one face of xhat and yhat is relevant per iteration
      # using k_t_top[2] we can basically just work in matrix world for a bit
      curr_x_projected <- xhat[, , k_t_top[2]] %*% curr_x_loadings
      curr_y_projected <- yhat[, , k_t_top[2]] %*% curr_y_loadings

      # perform deflation
      # the calculation of the x regression coefficient and deflation step of
      # xhat is the same regardless of pls mode
      curr_x_reg_coef <- crossprod(xhat[, , k_t_top[2]], curr_x_projected) /
        as.numeric(crossprod(curr_x_projected, curr_x_projected))

      xhat[, , k_t_top[2]] <- xhat[, , k_t_top[2]] -
        tcrossprod(curr_x_projected, curr_x_reg_coef)

      # the calculation of the y regression coefficient and deflation step of
      # yhat differs depending on the pls mode
      if (mode == "canonical") {
        curr_y_reg_coef <- crossprod(yhat[, , k_t_top[2]], curr_y_projected) /
          as.numeric(crossprod(curr_y_projected, curr_y_projected))

        yhat[, , k_t_top[2]] <- yhat[, , k_t_top[2]] -
          tcrossprod(curr_y_projected, curr_y_reg_coef)
      }

      if (mode == "regression") {
        curr_y_reg_coef <- crossprod(yhat[, , k_t_top[2]], curr_x_projected) /
          as.numeric(crossprod(curr_x_projected, curr_x_projected))

        yhat[, , k_t_top[2]] <- yhat[, , k_t_top[2]] -
          tcrossprod(curr_x_projected, curr_y_reg_coef)
      }

      x_loadings[, i] <- curr_x_loadings
      y_loadings[, i] <- curr_y_loadings
      x_projected[, i] <- curr_x_projected
      y_projected[, i] <- curr_y_projected
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

  } else {
    stop("Unexpected error in tpls, check 'mode' parameter input")
  }
}
