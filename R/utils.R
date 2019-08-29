### TODO drop deprecated utils

## ----------- .formula_checker -----------
#' Check formula's class and ensure non-NULL X and Y are not provided with it
#'
#' Gets a call which includes a $formula entry and expects it to be of form
#' \code{Y~X} or \code{Y~X1+X2+...} (for block), it also checks that $X and $Y are NULL.
#' 
#' @param mc A call list
#'
#' @return Exception handler
#'
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

## ----------- .check_data_assay ----------- 
#' Check validity of (data = data, X = assay) in call
#' @param mc A call object
#' @param vcs Valid classes for 'data'
#' @noRd
#' @importFrom SummarizedExperiment assay assays
# .check_data_assay <- function(mc, ## the call
#                               vcs = c("MultiAssayExperiment",
#                                       "SummarizedExperiment",
#                                       "SingleCellExperiment")){
#   data <- eval.parent(mc$data, 2L)
#   X <-    eval.parent(mc$X, 2L)
#   ## check that data is provided
#   if (isNULL(data)) {
#     .stop(.subclass = "inv_data",
#           "'X' is character but 'data' containing 'X' assay not provided. See ?pca.")
#     
#     ## check it is of valid class
#   } else if ( !class(try(data)) %in%  vcs) {
#     .stop(.subclass = "inv_data",
#           message = paste0("'data' must be of class: ",
#                            paste0(vcs, collapse = ", or ")))
#   }
#   
#   ## if X is not a valid assay throw appropriate error
#   if (!X %in% tryCatch(names(assays(data)), error = function(e) e)) {
#     .stop(.subclass = "inv_assay",
#           message = "'X' is not a valid assay name from 'data', it should be one of: ",
#           paste0(names(assays(data)), collapse = ", "))
#   }
# }

## ----------- .pcaMethodsHelper  ----------- 
#' Adjust methods for pca family and call internal
#'
#' @param mc The call object containing 'data', 'X', ... .
#' @param fun The function originally called, pca, ipca, spca, ... .
#' @param pframe How many levels higher should the arguments be evaluated.
#' Assumes the internal name in parent.frame is .fun (.pca, .ipca, ... .)
#'
#' @return Evaluated call to internal (.fun)
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

## ----------- .ret_call ----------- 
## function to add the call to the result if asked
# .ret_call <-
#   function(result = list(), ## result from internal
#            mcr, ## always evaluated mc
#            fun.name = 'pca') { ## name of the function
#     
#     if (isTRUE(mcr$ret.call) ) {
#       {
#         mcr[[1]] = as.name(fun.name)
#       }
#       ## makes expect_identical() easy to have call as first arg and just
#       ## drop it using pca.res[-1]
#       
#       result <- structure(c(list(call = mcr), result), class = class(result))
#     }
#     return(result)
#   }
# ## ----------- .switch_arg_names ----------- #TODO drop it?
# ## customised from scran/R/utils_other.R
# .switch_arg_names <- function(old.val, new.val, msg=NULL) {
#   if (!is.null(old.val)) {
#     old.arg <- deparse(substitute(old.val))
#     new.arg <- deparse(substitute(new.val))
#     if ( is.null(msg) ) {
#       msg <- sprintf("The use of %s as data matrix is deprecated and is only used for assay names
#                      from'data' now, use %s instead.
#                      See documentation for details.", old.arg, new.arg)
#     }
#     .Deprecated(new = new.arg, old = old.arg)
#     old.val
#   } else {
#     new.val
#   }
# }


## ----------- .call_internal ----------- 
## function to add the call to the result if asked
# .call_internal <-
#   function(mc, ## always match.call()
#            fun.name = 'pca', ## name of the function
#            pframe = 3L) {  ## how many frames above to evaluate in?
#     mc[-1L] <- lapply(mc[-1L], eval.parent)
#     mc$data <- mc$formula <- NULL 
#     mc[[1L]] <- as.name(sprintf(".%s", fun.name))
#     result <- eval(mc)
#     
#     if ( isTRUE(mc$ret.call) ) {
#       {
#         mcr[[1]] = as.name(fun.name)
#         mcr[-1L] <- lapply(mcr[-1L], function(x) eval.parent(x, n = pframe))
#       }
#       ## makes expect_identical() easy to have call as first arg and just
#       ## drop it using pca.res[-1]
#       
#       result <- structure(c(list(call = mcr), result), class = class(result))
#     }
#     return(result)
#   }


## ----------- .block_get_xy ----------- 
#' Create matrices from names of assays and a data object
#' 
#' get MAE data and a call list containing either character X and Y,
#' or a formula of form Y~X1+X2+... and check for validity of X (assay names) 
#' and Y and return list of matrices of X and matrix Y in mc$X and mc$Y, 
#' getting rid of  $data and/or $formula args
#' 
#' @param mc A list with at least (data=MAE, formula = Y~X1+X2+...)
#' @noRd
## TODO does it?
#' @importFrom SummarizedExperiment colData
#' @importFrom MultiAssayExperiment experiments
#' @importFrom SummarizedExperiment assay

.block_get_xy <- function(mc){
  ## ---- if data and X , Y given, expect X,Y
  if (is.null(mc$formula)) { 
      if (isNULL(mc$X) || isNULL(mc$Y)) {
        .stop("X and/or Y is NULL, should be vector and character of assay names, respectively.")
      }
    ## if a X is list, unlist it
      if (is(mc$X, "list")) { 
        mc$X <- unlist(mc$X)
      }
    if (!is.null(mc$indY)) { ## if indY provided
      mc$Y <- mc$X[indY]
      mc$X <- mc$X[-indY]
      }
      assay.names <- c(mc$Y, mc$X)
      ## expect it to be character
      if (!is(assay.names, "character"))
        .stop("X and Y must be characters")
    
  } 
  ## ----- if data and formula given, expect NULL X and Y, and set them from formula
  else {  ## formula
    .formula_checker(mc, block = TRUE) ## check formula validity
    ## esnure X and Y are NULL
    if ( !(isNULL(mc$X) && isNULL(mc$Y)) )
      .stop(message = "Where 'data' and 'formula' are provided 'X' and 'Y' should be NULL.", 
            .subclass = "inv_signature")
    
    assay.names <- rownames(attr(terms(mc$formula), "factors"))
  }
  ## ----- check assay names and create the matrices
  missing.assays <-  assay.names[!assay.names %in% names(experiments(mc$data))]
  if (length(missing.assays))
    .stop(paste("Not valid assay name(s) from data: ", paste(missing.assays, collapse = ", ")), .subclass = "inv_XY")
  ## get X and Y in X first
  mc$X <- list()
  for (asy in assay.names) {
    mc$X[[asy]] <- as.matrix(t(assay(mc$data, asy)))
  }
  
  ## separate Y
  mc$Y <- mc$X[[1]]
  mc$X <- mc$X[-1]
  
  mc$data <- mc$formula <- NULL
  mc
}

## ----------- .blockDA_get_xy ----------- 
#' Create matrices from names of assays and a data object
#' 
#' get MAE data and a call list containing either character X and Y,
#' or a formula of form Y~X1+X2+... and check for validity of X (assay names) 
#' and Y and return list of matrices of X and matrix Y in mc$X and mc$Y, 
#' getting rid of  $data and/or $formula args
#' 
#' @param mc A list with at least (data=MAE, formula = Y~X1+X2+...)
#' @noRd
#' @importFrom MultiAssayExperiment MatchedAssayExperiment
## TODO does it?
#' @importFrom SummarizedExperiment colData
#' @importFrom MultiAssayExperiment experiments
#' @importFrom SummarizedExperiment assay

# .blockDA_get_xy <- function(mc){
#   ## ---- if data and X , Y given, expect X,Y
#   # TODO make this internal and re-use
#   if (is.null(mc$formula)) { 
#     if (isNULL(mc$X) || isNULL(mc$Y)) {
#       .stop("X should be vector of assay names, and Y a column data name from data.")
#     }
#     ## if a X is list, unlist it
#     if (is(mc$X, "list")) { 
#       mc$X <- unlist(mc$X)
#     }
#     ## expect it to be character
#     if (!is(c(mc$X, mc$Y), "character"))
#       .stop("X and Y must be characters")
#     Y.input <- mc$Y
#     X.input <- mc$X
#   } 
#   ## ----- if data and formula given, expect NULL X and Y, and set them from formula
#   else {  ## formula
#     .formula_checker(mc, block = TRUE) ## check formula validity
#     ## esnure X and Y are NULL
#     if ( !(isNULL(mc$X) && isNULL(mc$Y)) )
#       .stop(message = "Where 'data' and 'formula' are provided 'X' and 'Y' should be NULL.", 
#             .subclass = "inv_signature")
#     
#     trms <- rownames(attr(terms(mc$formula), "factors"))
#     Y.input <- trms[1]
#     X.input <- trms[-1]
#   }
#   ## ----- check assay names and create the matrices
#   missing.assays <-  X.input[!X.input %in% names(experiments(mc$data))]
#   missing.coldata <-  !(Y.input %in% names(colData(mc$data)))
#   if (length(missing.assays))
#     .stop(paste("Not valid assay name(s) from data: ", paste(missing.assays, collapse = ", ")), .subclass = "inv_XY")
#   
#   if (isTRUE(missing.coldata))
#     .stop(sprintf("LHS of formula should be a valid column data from 
#                   'data', (i.e. one of: %s) . %s is not.", 
#                   paste(names(colData(mc$data)), collapse = ", " ), sQuote(Y.input)), 
#           .subclass = "inv_XY")
#   
#   ## get matched samples
#   mc$data <- MatchedAssayExperiment(mc$data)
#   ## get X 
#   mc$X <- list()
#   for (asy in X.input) {
#     mc$X[[asy]] <- as.matrix(t(assay(mc$data, asy)))
#   }
#   
#   ## separate Y
#   mc$Y <- colData(mc$data)[,Y.in]
#   
#   mc$data <- mc$formula <- NULL
#   mc
# }



## ----------- .check_sig_ANY ----------- 
#' Handle arguments possibly from old code
#'
#' Arguments pased to 'ANY' can only include numeric X and Y.
#' The X, Y arguments should be named, otherwise they'll be taken as data
#' 
#' 
#' @param mc matched call of form (data=MAE, X=X_name, Y=Y_name, formula=NULL)
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
#' @param mae 
#' @importFrom MultiAssayExperiment MatchedAssayExperiment
#' @noRd
.matched_samples <- function(mae) {
  message("Keeping matching samples in MultiAssayExperiment object ...")
  mae <- MatchedAssayExperiment(mae)
}



