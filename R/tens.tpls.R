# ==============================================================================
# Tensor pls generalization; developed @ Melbourne Integrative Genomics based
# on Kilmer's tensor m-product algebra
# ==============================================================================

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
  bpparam = NULL
) {
  allowed_modes <- c("canonical", "regression")
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

  if (length(dim(x)) != 3) {
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

  tsvdm_decomposition_xty <- tsvdm(
    facewise_product(ft(x), y, m = m, minv = minv),
    transform = FALSE, # by default tsvdm will transform input into m space
    full_frontal_slices = FALSE,
    svals_matrix_form = TRUE
  )

  # process ncomp input, much simple than tpca as we only accept integer input
  # bltodo: investigate non integer input options? explained variance semantics?
  if (is.null(ncomp)) {
    ncomp <- min(n, p, q) * t
  }

  for (i in seq_len(ncomp)) {
    
  }
}
