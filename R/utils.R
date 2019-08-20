### TODO drop deprecated utils
## ----------- .get_xy ----------- 
#' Create matrices from names of assays and a data object
#' 
#' get MAE data and a call list containing either character X and Y,
#' or a formula of form Y~X and check for validity of X (assay name) and 
#' Y (assay/coldata name) and return matrices of X and Y in mc$X and mc$Y, 
#' getting rid of  $data and/or $formula args
#' 
#' @param mc A list with at least (data=MAE, X=X_name, Y=Y_name)=
#' @noRd
#' @importFrom MultiAssayExperiment complete.cases
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment assay

.get_xy <- function(mc){ 
  mcc <- mc ## copy before change
  ## ensure X and Y are character
  tryCatch({
    if(any(class(c(mc$X, mc$Y))!="character")) 
      stop("mc must have character X and Y", call. = FALSE)
  }, error=function(e) 
    message("oops! something went wrong! Please inputs again or contact us if you had no luck!"))
  
  ## ensure X and Y are valid
  if(!mc$X %in% names(assays(mc$data))) 
    .stop(.subclass = "inv_XY",
          message = " 'X' is not a valid assay from 'data'")
  
  ## --- if 'Y' is a colData
  if(mc$Y %in% names(colData(mc$data))) {
    if(mc$Y %in% names(assays(mc$data))) 
      .stop(.subclass = "inv_XY", 
            message = paste0(mcc$Y, " matches to both colData and assay
                             in 'data', change its name in one and continue."))
    ## ----- if Y is a colData column subset it using X samples
    ## keep X assay's colData & change DataFrame to data.frame
    Xcoldata <- suppressMessages( as.data.frame(colData(mc$data[,,mc$X]))) 
    mc$Y <- Xcoldata[,mc$Y] ## keep the coldata desired for Y
    ## if Y not numeric coerce
    if (! typeof(mc$Y) %in% c("numeric", "integer")) {
      if (typeof(mc$Y) == "factor") {
        warning(paste0("The column data ", sQuote(mcc$Y)," is a factor,
                       coercing to a numeric vector with names..."))
        mc$Y <- structure(as.numeric(mc$Y),
                          names = as.character(mc$Y), 
                          class = "numeric")
        
      } else if(typeof(mc$Y)=="character"){
        ## if Y is a character colData and the number of unique terms are less than total,
        ## coerce it to factor and then numeric with a warning
        if(length(unique(mc$Y)) <  length(mc$Y) ){
          .warning("char_Y", message = paste0("The column data ",mcc$Y, " is character vector, coercing to factor and then named numeric for pls"))

          mc$Y <- structure(as.numeric(as.factor(mc$Y)),
                            names=mc$Y, class="numeric")
        } else {
          .stop(.subclass = "inv_XY", message = paste0(" 'Y' is not a numeric/integer column (or a factor coercible to numeric)"))
        }
      }

    }
    ## if all is well with Y
    mc$X <- assay(mc$data, mc$X)
    ## ----- If Y is assay name
  } else if(mc$Y %in% names(assays(mc$data))){
    mc$data <- mc$data[,complete.cases(mc$data[,,c(mc$X, mc$Y)])] ## keep complete data b/w two assays
    mc$X <- assay(mc$data, mc$X)
    mc$Y <- assay(mc$data, mc$Y)
  } else {.stop(.subclass = "inv_XY", message = paste0(sQuote(mcc$Y), " is not an assay or column data from the MAE object" ))}
  mc$X <- t(as.matrix(mc$X))
  mc$Y <- as.matrix(mc$Y)

  if(!1 %in% dim(mc$Y)){ ## if Y is matrix transpose it
    mc$Y <- t(mc$Y)
  }
  return(mc)
}

## ----------- .sformula_checker -----------  ## single-RHS formula checker
#' Check formula's class and ensure non-NULL X and Y are not provided with it
#'
#' Gets a call which includes a $formula entry and expects it to be of form
#' \code{Y~X}, it also checks that $X and $Y are NULL.
#' 
#' @param mc A call list
#'
#' @return Exception handler
#'
#' @noRd
.sformula_checker <- function(mc){
  if(class(try(mc$formula))!="formula")
    .inv_sformula()
  ## check formula is Y~X
  if(any(sapply(as.list(mc$formula), length)!=1))
    .inv_sformula()
  ## X and Y must be NULL
  if(!all(sapply(mc[c("X", "Y")], is.null)))
    .inv_signature()
}

## ----------- .plsMethodsHelper ----------- 
####  get call list including potentiall X, Y, formula, and data and retain only valid X and Y
.plsMethodsHelper <- function(mc){
  mc[c('data', 'formula')] <- lapply( mc[c('data', 'formula')], eval.parent)
  mcc <- mc ## copy so can change mc but keep the call for .get_xy
  # expectedArgs <- c('X', 'Y', 'formula', 'data')
  # mc[expectedArgs] <- lapply(mc[expectedArgs], eval.parent)
  
  ##============================= if data
  if (!is.null(try(mc$data)
  )) {
    ## ensure it's MAE class
    if (class(try(mc$data)
    )  !=  "MultiAssayExperiment") {
      .inv_mae(mc$data)
    }
    
    ##--------------- if data & formula≠NULL
    ##--- i) if (data,formula) given change it to X and Y matrices
    if (class(try(mc$formula)
    )  !=  "NULL") {
      .sformula_checker(mc = mc)
      mc[c("X", "Y")] <- as.character(as.list(mc$formula)[3:2])
      mc <- .get_xy(mc, mcc)
    }
    ##--------------- if data & formula=NULL
    else {
      ## check X and Y exist
      if(any(sapply(mc[c("X", "Y")], function(xy) {class(try(xy))=="NULL"})))
        .inv_assay()
      ## in case they're stored in variables
      mc[c("X", "Y")] <- lapply( mc[c("X", "Y")], eval.parent)
      ## ensure it is a single character
      if(any(sapply( mc[c("X", "Y")], length)!=1))
        .stop(.subclass = "inv_XY", message = "'X' and 'Y' must be assay names from 'data'")
      mc <- .get_xy(mc, mcc)
    }
    ##--- if data, X and Y , expect X and Y to be assays and change them to matrices
    # else if(class(try(mc$formula))!="NULL"){
    #   ## if formula not a fomrula class, expect it to be NULL and X and Y to be assay/colData
    #   if(class(try(mc$formula))!="NULL")
    #     .stop(.subclass = "inv_formula", message = "'formula' must be a formula object of form Y~X")
    #
    # }

  }
  ##============================= if data=NULL and formula≠NULL
  else if (class(try(mc$formula))!="NULL"){
    mc$formula <- as.formula(mc$formula)
    .sformula_checker(mc=mc)
    mc[c('Y','X')] <- as.list(mc$formula)[2:3]
  }
  mc[c("X", "Y")] <- lapply( mc[c("X", "Y")], eval.parent)
  mc$data <- mc$formula <- NULL
  return(mc)
}

## ----------- .check_data_assay ----------- 
#' Check validity of (data = data, X = assay) in call
#' @param mc A call object
#' @param vcs Valid classes for 'data'
#' @noRd
#' @importFrom SummarizedExperiment assay assays
.check_data_assay <- function(mc, ## the call
                              vcs = c("MultiAssayExperiment",
                                      "SummarizedExperiment",
                                      "SingleCellExperiment")){
  data <- eval.parent(mc$data, 2L)
  X <-    eval.parent(mc$X, 2L)
  ## check that data is provided
  if (missing(data) || is.null(data)) {
    .stop(.subclass = "inv_data",
          "'X' is character but 'data' containing 'X' assay not provided. See ?pca.")
    
    ## check it is of valid class
  } else if ( !class(try(data)) %in%  vcs) {
    .stop(.subclass = "inv_data",
          message = paste0("'data' must be of class: ",
                           paste0(vcs, collapse = ", or ")))
  }
  
  ## if X is not a valid assay throw appropriate error
  if (!X %in% tryCatch(names(assays(data)), error = function(e) e)) {
    .stop(.subclass = "inv_assay",
          message = "'X' is not a valid assay name from 'data', it should be one of: ",
          paste0(names(assays(data)), collapse = ", "))
  }
}

## ----------- .pcaMethodsHelper  ----------- 
#' Adjust methods for pca family and call internal
#'
#' @param mc The call object containing 'data', 'X', ... .
#' @param fun The function originally called, pca, ipca, spca, ... .
#' Assumes the internal name in parent.frame is .fun (.pca, .ipca, ... .)
#'
#' @return Evaluated call to internal (.fun)
#' @noRd
.pcaMethodsHelper <- function(mc, fun='pca'){
  mc[-1L] <- lapply(mc[-1L], function(x) eval.parent(x, n = 2L))
  .check_data_assay(mc)
  # mcr <- mc ## to return as output
  # mcr[[1L]] <- as.name(fun) ## returned call to have 'pca' as function
  mc$X <- t(assay(mc$data, mc$X)) ## change X into matrix of assay name
  mc$data <- NULL ## remove data as not needed in internal
  mc[[1L]] <- as.name(sprintf(".%s", fun)) ## add internal to the call's function
  result <- eval.parent(mc) ## evaluate the call
  # if(isTRUE(mc$ret.call))
  # result$call <- mcr ## replace the returned call from internal by the current one
  return(result)
}

## ----------- .pcaEntryChecker ----------- 
#' PCA Entry Checker
#'
#' checks X, ncomp, scale, center, max.iter, and tol. Makes adjustments to
#' them in the parent environment if necessary. for IPCA, 'center' argument does
#' not apply (is always TRUE) and NA's are not allowed.
#' Note: it possibly replaces X, ncomp, and max.iter in parent.frame using <<-.
#' Note: if the error classes are too many, you can batch replace and simplify.
#' @noRd
.pcaEntryChecker <-
  function(mc,fun='pca') {
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
  if (is.null(mc$max.iter) ||!is.numeric(mc$max.iter) || 
      mc$max.iter < 1 || !is.finite(mc$max.iter))
    .stop("invalid value for 'max.iter'.")

  #-------- tol
  if (is.null(mc$tol) || !is.numeric(mc$tol) || mc$tol < 0 || !is.finite(mc$tol))
    .stop("invalid value for 'tol'.")
  
  return(mc)
  }

## ----------- .call_return  ----------- 
## function to add the call to the result if asked
.call_return <-
  function(result = list(), ## result from internal
           ret.call = FALSE, ## whether call should be returned
           mcr, ## always match.call()
           fun.name = 'pca') { ## name of the function
    
  if (isTRUE(ret.call)) {
    {
      mcr[[1]] = as.name(fun.name)
      mcr[-1L] <- lapply(mcr[-1L], function(x) eval.parent(x, n = 2))
    }
    ## makes expect_identical() easy to have call as first arg and just
    ## drop it using pca.res[-1]
    result <- c(list(call = mcr), result)
  }
    return(result)
}

## ----------- .switch_arg_names ----------- 
## customised from scran/R/utils_other.R
.switch_arg_names <- function(old.val, new.val, msg=NULL) {
  if (!is.null(old.val)) {
    old.arg <- deparse(substitute(old.val))
    new.arg <- deparse(substitute(new.val))
    if ( is.null(msg) ) {
      msg <- sprintf("The use of %s as data matrix is deprecated and is only used for assay names
                     from'data' now, use %s instead.
                     See documentation for details.", old.arg, new.arg)
    }
    .Deprecated(new = new.arg, old = old.arg)
    old.val
  } else {
    new.val
  }
}

