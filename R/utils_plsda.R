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
    rownames(Y.mat) <- rownames(mc$X)
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

.check_plsda_block <- function(mc) {
  Y <- mc$Y
  X <- mc$X
  indY <- mc$indY
  
  if (!isNULL(Y))
  {
    if (is.null(dim(Y)))
    {
      Y = factor(Y)
    } else {
      stop("'Y' should be a factor or a class vector.")
    }
    
    if (nlevels(Y) == 1)
      stop("'Y' should be a factor with more than one level")
    
    Y.input = Y
    Y = unmap(Y)
    colnames(Y) = levels(Y.input)
    rownames(Y) = rownames(X[[1]])
  } else if (!isNULL(indY)) {
    temp = X[[indY]]
    #not called Y to not be an input of the wrapper.sparse.mint.block
    if (is.null(dim(temp)))
    {
      temp = factor(temp)
    } else {
      stop("'Y' should be a factor or a class vector.")
    }
    
    if (nlevels(temp) == 1)
      stop("'X[[indY]]' should be a factor with more than one level")
    
    Y.input = temp
    X[[indY]] = unmap(temp)
    colnames(X[[indY]]) = levels(Y.input)
    rownames(X[[indY]]) = rownames(X[-indY][[1]])
    
  } else if (isNULL(indY)) {
    stop("Either 'Y' or 'indY' is needed")
  }
  mc[c("X", "Y", "indY", "Y.input")] <- list(X, Y, indY, Y.input)
  mc
}