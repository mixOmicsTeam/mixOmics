# ==============================================================================
# Tensor block pls generalization; developed @ Melbourne Integrative Genomics
# based on Kilmer's tensor m-product algebra
# ==============================================================================

#' Tensor block PLS
#'
#' Developed @ Melbourne Integrative Genomics.
#'
#' @param x A list of tensor inputs (termed 'blocks') measured on the same
#' samples.
#' @param y Tensor input.
#' @param ncomp The estimated number of components. ncomp must be explicitly set
#' as an integer in tpls.
#' @param design Numeric matrix of size length(x) x length(x) with values
#' between 0 and 1; the i,j entry in the matrix indicates the strength of the
#' relationship to be modelled between the i-th and j-th blocks. A value of 0
#' indicates no relationship, with 1 being the maximum. Alternatively, one can
#' input "null" for a fully disconnected design, or "full" for a fully connected
#' design, or a single scalar value between 0 and which will designate all
#' off-diagonal elements of a fully connected design.
#' @param m A function which applies an orthogonal tensor tubal transform.
#' @param minv The inverse of m.
#' @param mode Currently supports tensor analogues of canonical, regression,
#' and svd PLS modes. Defaults to "regression" mode.
#' @param center If set to false, the data tensor will not be centralized into
#' Mean Deviation Form (see Mor et al. 2022). By default, the mean horizontal
#' slice of the input tensor(s) are subtracted, so that all of the horizontal
#' slices sum to 0, analgous to centering matrix data.
#' @param bpparam A \linkS4class{BiocParallelParam} object indicating the type
#' of parallelisation. Does not have any effect if transform functions
#' explicitly set using \code{m}, \code{minv}.
#' @author Brendan Lu
#' @export
block_tpls <- function(
  x,
  y,
  ncomp = NULL,
  design,
  m = NULL,
  minv = NULL,
  mode = "regression",
  center = TRUE,
  bpparam = NULL
) {
  # bltodo: add docs link to existing block.pls for user reference
  NULL
}
