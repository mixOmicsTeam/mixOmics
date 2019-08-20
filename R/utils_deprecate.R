#' deprecate unnamed first matrix argument
#' 
#' Deprecate warning for pca(matrix) {unnamed first argument which is
#' takes as data and not X} which now should be pca(data=NULL, X=matrix)
#' @noRd
.deprecate_ufma <- function(
  message='mixOmics arguments have changed and this code will not work
  in near future (you can use named arguments to avoid this),
  please see documentation', .subclass='deprecated'){
  .warning(message = message,
           .subclass = .subclass)
} 

