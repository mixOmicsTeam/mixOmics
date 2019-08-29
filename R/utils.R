### TODO drop deprecated utils

## ----------- .formula_checker -----------
#Check formula's class and ensure non-NULL X and Y are not provided with it
#
#'Gets a call which includes a $formula entry and expects it to be of form
#\code{Y~X} or \code{Y~X1+X2+...} (for block), it also checks that $X and $Y are NULL.
#
#@param mc A call list
#
#@return Exception handler
#
#' @noRd
.formula_checker <- function(mc, block=FALSE){
  fl <- vapply(as.list(mc$formula)[-1L], length, 1) ## formula list lengths
  
  if (isTRUE(block)) {
    ## check formula is Y~X1+X2+...
    if ( length(fl) != 2 || fl[1] != 1 || fl[2] < 2)
      .inv_bformula()
  } else {
    ## check formula is Y~X
    if ( length(fl) != 2 || fl[1] != 1 || fl[2] != 1)
      .inv_sformula()
  }
  ## X and Y must be NULL
  if (!all(vapply(mc[c("X", "Y")], is.null, TRUE)))
    .inv_signature()
}

## ----------- .pcaMethodsHelper  ----------- 
# Adjust methods for pca family and call internal
#
# @param mc The call object containing 'data', 'X', ... .
# @param fun The function originally called, pca, ipca, spca, ... .
# @param pframe How many levels higher should the arguments be evaluated.
# Assumes the internal name in parent.frame is .fun (.pca, .ipca, ... .)
#
# @return Evaluated call to internal (.fun)
#' @noRd
.pcaMethodsHelper <- function(mc, fun='pca'){
  mc <- .get_x(mc)
  mc$data <- NULL ## remove data as not needed in internal
  mc[[1L]] <- as.name(sprintf(".%s", fun)) ## add internal to the call's function
  eval.parent(mc) ## evaluate the call
}

## ----------- .pcaEntryChecker ----------- 
#' PCA Entry Checker
#'
#' checks X, ncomp, scale, center, max.iter, and tol and return possibly
#' adjusted call. for IPCA, 'center' argument does
#' not apply (is always TRUE) and NA's are not allowed.
#' Note: if the error classes are too many, you can batch replace and simplify.
#' @noRd
.pcaEntryChecker <-
  function(mc,fun='pca') {
    if(class(try(mc$ret.call)) != "logical") stop("'ret.call' must be logical")
    ## data.frame that specifies which arguments need checks in each function
    checks_df <- data.frame(
      row.names = c('check.keepX', 'check.center', 'check.NA'),
      pca = c(FALSE, TRUE, FALSE),
      spca = c(TRUE, TRUE, FALSE),
      ipca = c(FALSE, FALSE, TRUE),
      sipca = c(TRUE, FALSE, TRUE))
    
  #-------- X
  if ( is.data.frame(mc$X) )
    mc$X <- as.matrix(mc$X)
  
  if (!(is.matrix(mc$X) && is.numeric(mc$X)))
    .stop(.subclass = 'inv_X', message = "'X' must be a numeric matrix.")
  
  if ( any(is.infinite(mc$X)) )
     .stop(message = "infinite values in data not allowed",
           .subclass = "inv_X")
     

  ## checking NAs, only for ICA
  if (checks_df['check.NA', fun]) {
    if (any(is.na(mc$X)))
      .stop(.subclass = 'inv_X',message = "data cannot contain missing values
      for 'ipca' functions")
  }


  #-------- ncomp
  if (is.null(mc$ncomp))
    mc$ncomp <- as.integer(min(dim(mc$X)))

  if ( !is.numeric(mc$ncomp) || mc$ncomp < 1 || !is.finite(mc$ncomp))
    .stop(.subclass = 'inv_ncomp', 
          message = "'ncomp' must be an integer not greater than min(dim(X))")

  if (mc$ncomp > min(dim(mc$X)))
    .stop(.subclass = 'inv_ncomp', 
          message = "'ncomp' must be an integer not greater than min(dim(X))")

  #-------- keepX
  if ( checks_df['check.keepX', fun] ) {
    if (is.null(mc$keepX))
      mc$keepX <- rep(ncol(mc$X), mc$ncomp)
    if (length(mc$keepX) != mc$ncomp)
      .stop(message = paste0("length of 'keepX' must be equal to 'ncomp', 
                             that is ", mc$ncomp, "."))
    if (any(mc$keepX > ncol(mc$X)))
      .stop(message = "each component of 'keepX' must be lower or equal 
            to smallest dimension of 'X'")
    }

  #-------- center
  ## no centering option in ICA
  if ( checks_df['check.center', fun] ) {
    if (!is.logical(mc$center))
    {
      if (!is.numeric(mc$center) || (length(mc$center) != ncol(mc$X)))
        .stop(message = "'center' should be either a logical value or a numeric
              vector of length equal to the number of columns of 'X'.")
    }
  }


  #-------- scale
  if (!is.logical(mc$scale))
  {
    if (!is.numeric(mc$scale) || (length(mc$scale) != ncol(mc$X)))
      .stop(message = "'scale' should be either a logical value or a 
            numeric vector of length equal to the number of columns of 'X'.")
  }

  #-------- max.iter
  if (is.null(mc$max.iter) || !is.numeric(mc$max.iter) || 
      mc$max.iter < 1 || !is.finite(mc$max.iter))
    .stop("invalid value for 'max.iter'.")

  #-------- tol
  if (is.null(mc$tol) || !is.numeric(mc$tol) || mc$tol < 0 || !is.finite(mc$tol))
    .stop("invalid value for 'tol'.")
  
  return(mc)
  }

## ----------- .call_return ----------- 
## TODO keep only one returner function
## function to add the call to the result if asked
## The aim is to reduce the methods code
## Handling it in here to avoid unnecessary copy of the fully evaluated call
## unless user asks for it. No full copy is made unless ret.call=TRUE.
.call_return <-
  function(result = list(), ## result from internal
           mcr, ## always match.call()
           fun.name = 'pca', ## name of the function
           pframe = 3L) {  ## how many frames above to evaluate in?
    
  if ( isTRUE(mcr$ret.call) ) {
    {
      mcr[[1]] = as.name(fun.name)
      mcr[-1L] <- lapply(mcr[-1L], function(x) eval.parent(x, n = pframe))
    }
    ## makes expect_identical() easy to have call as first arg and just
    ## drop it using pca.res[-1]
    
    result <- structure(c(list(call = mcr), result), class = class(result))
  }
    return(result)
}

## ----------- .check_sig_ANY ----------- 
# Handle arguments possibly from old code
#
# Arguments pased to 'ANY' can only include numeric X and Y.
# The X, Y arguments should be named, otherwise they'll be taken as data
# 
# 
# @param mc matched call of form (data=MAE, X=X_name, Y=Y_name, formula=NULL)
#' @noRd
.check_sig_ANY <- function(mc, fun = "block.pls") {
  ## TODO ensure all referrals to help by this handler are doc backed.
  ## ---- if data matches what we expect X to be, it might be old code
  if ( class(try(mc$data)) %in% c("data.frame", "matrix", "list") ) {
    .stop(message = sprintf("mixOmics arguments have changed.
              Please carefully read the documentation and try to use named 
              arguments such as %s(X=list(), ...) as opposed to %s(mat, ...).",
                            fun, fun) ,
          .subclass = "defunct")
  }
  ## ---- if still data is not NULL and not above, it's messed up
  if ( !isNULL(mc$data)) { ## data must be NULL
    .stop(message = "data should be a MultiAssayExperiment class, or NULL",
          .subclass = "inv_signature")
  }
  ## ---- to ensure formula class never makes it to 'ANY' - dear Lord!
  if ( !isNULL(mc$formula) ) { ## formula must be NULL
    msg <- sprintf("With numerical X and Y, formula should not be provided. 
              See ?%s", fun)
    .stop(message = msg,
          .subclass = "inv_signature")
  }
}

## ----------- .matched_samples ----------- 
#' Keep onlymatched samples in MAE object
#'
#' Creating this so I would not have to import every time
#' 
#' 
# @param mae 
#' @importFrom MultiAssayExperiment MatchedAssayExperiment
#' @noRd
.matched_samples <- function(mae) {
  message("Keeping matching samples in MultiAssayExperiment object ...")
  mae <- MatchedAssayExperiment(mae)
}



