.check_plsda <- function(mc) {
  if (is.null(mc$multilevel))
  {
    if (is.null(mc$Y))
      stop("'Y' has to be something else than NULL.")
    
    if (is.null(dim(mc$Y)))
    {
      mc$Y = factor(mc$Y)
    } else {
      stop("'Y' should be a factor or a class vector.")
    }
    
    if (nlevels(mc$Y) == 1)
      stop("'Y' should be a factor with more than one level")
    
    Y.mat = unmap(mc$Y)
    colnames(Y.mat) = levels(mc$Y)
    mc$Y <- Y.mat
    
  } else {
    # we expect a vector or a 2-columns matrix in 'Y' and the repeated measurements in 'mc$multilevel'
    mc$multilevel = data.frame(mc$multilevel)
    
    if ((nrow(X) != nrow(mc$multilevel)))
      stop("unequal number of rows in 'X' and 'mc$multilevel'.")
    
    if (ncol(mc$multilevel) != 1)
      stop("'mc$multilevel' should have a single column for the repeated measurements, other factors should be included in 'Y'.")
    
    if (!is.null(ncol(mc$Y)) && !ncol(mc$Y) %in% c(0,1,2))# mc$multilevel 1 or 2 factors
      stop("'Y' should either be a factor, a single column data.frame containing a factor, or a 2-columns data.frame containing 2 factors.")
    mc$multilevel = data.frame(mc$multilevel, mc$Y)
    mc$multilevel[, 1] = as.numeric(factor(mc$multilevel[, 1])) # we want numbers for the repeated measurements
  }
  mc
}