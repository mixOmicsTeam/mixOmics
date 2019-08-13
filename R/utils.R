## ------- put the character b/w two single quotation marks ----
squote <- function(char) paste0("'",char,"'")

##-----------------------------------------------------------------------------------------------------------------------
## ---- check that EXP can be evaluated and returned, if not throw error with error_message error  ----
.checkExpr <- function(EXP, error_message = paste0(as.expression(substitute(EXP)), " cannot be evaluated")){
  formals(stop)$call. <-FALSE
  out <- tryCatch(EXP, error = function(e) e)
  if ("simpleError" %in% class(out)) stop(error_message)
  return(out)
}

##-----------------------------------------------------------------------------------------------------------------------
## ---- get a MAE object, assay name/index, and the call list, and return the data matrix for MAE methods ----
## args.list contains X entry
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment assays
.getDM <- function(X, assay){ ## MAE to data.matrix
  args.list <- match.call()[-1]
  ## ---------- get the assay name from either name  or index provided and check
  if (is.numeric(assay)){
    if(assay-floor(assay)!=0) .stop(.subclass = "inv_xy", message = paste0(assay, " is not a valid assay index. Use an integer."))
    if(assay<1 | assay >length(assays(X))) .stop(.subclass = "inv_xy", message = paste0("assay index must be positive integer smaller than or equal to the number of assays in ",
                                                                                                     args.list["X"]," (i.e. 1:",length(assays(X)),")"))
    assay <- names(assays(X))[assay]
  } else if(is.null(assay)){
    .stop(.subclass = "inv_xy", message = paste0("Please provide an assay from object ", args.list["X"]))
    } else if(is.na(assay)){
      .stop(.subclass = "inv_xy", message = "assay cannot be NA ")
    } else if(!is.character(assay)){
    ## if 'assay' it none of acceptable forms
    .stop(.subclass = "inv_xy", message = paste0("'assay' must be either an assay name or index from ", args.list["X"]))
  }

  ## ---------- use assay name to get the data matrix
  if(! assay %in% names(assays(X)))
    .stop(.subclass = "inv_xy", message = paste0(assay, " is not a valid assay from ","'",args.list["X"],"'"))
  ## transpose and create a data matrix
  X <- .checkExpr(EXP = data.matrix(t(assay(X,assay))),
                                        error_message = paste0("could not create a data matrix from assay '", assay, "' in ", args.list["X"]) )
  if(!is.numeric(as.matrix(X))) .stop(.subclass = "inv_xy", message = paste0("The ", assay, " assay contains non-numeric values"))
  return(X)
}


## ------------------------------------------
## ----- get MAE data and a call list containing character X and Y, check for validity of X (assay name) and Y (assay/coldata name)
## and return matrices of X and Y in mc$X and mc$Y, getting rid of data and/or formula args
#' @importFrom MultiAssayExperiment complete.cases
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment assay

.names2mat <- function(mc, mcc){ ## mc is list(data=MAE, X=X_name, Y=Y_name)
  tryCatch({
    if(any(class(c(mc$X, mc$Y))!="character")) stop("mc must have character X and Y", call. = FALSE)
  }, error=function(e) message("oops! something went wrong! please check X and Y again or contact us if you had no luck!"))
  if(!mc$X %in% names(assays(mc$data))) .stop(.subclass = "inv_xy", message = " 'X' is not a valid assay from 'data'")
  ## --- if 'Y' is a colData
  if(mc$Y %in% names(colData(mc$data))){
    if(mc$Y %in% names(assays(mc$data))) .stop(.subclass = "inv_xy", message = paste0(mcc$Y, " matches to both colData and assay in 'data', change its name in one and continue."))
    ## ----- if Y is a colData column subset it using X samples
    Xcoldata <- suppressMessages( as.data.frame(colData(mc$data[,,mc$X]))) ## keep X assay coldata - DataFrame to data.frame
    mc$Y <- Xcoldata[,mc$Y] ## keep the coldata desired for Y
    ## if Y not numeric
    if(! typeof(mc$Y) %in% c("numeric","integer")){
      if(typeof(mc$Y)=="factor") {
        warning(paste0("The column data ", squote(mcc$Y)," is a factor, coercing to a numeric vector with names..."))
        mc$Y <- structure(as.numeric(mc$Y),
                          names=as.character(mc$Y), class="numeric")
      } else if(typeof(mc$Y)=="character"){
        ## if Y is a character colData and the number of unique terms are less than total,
        ## coerce it to factor and then numeric with a warning
        if(length(unique(mc$Y)) <  length(mc$Y) ){
          .warning("char_Y", message = paste0("The column data ",mcc$Y, " is character vector, coercing to factor and then named numeric for pls"))

          mc$Y <- structure(as.numeric(as.factor(mc$Y)),
                            names=mc$Y, class="numeric")
        } else {
          .stop(.subclass = "inv_xy", message = paste0(" 'Y' is not a numeric/integer column (or a factor coercible to numeric)"))
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
  } else {.stop(.subclass = "inv_xy", message = paste0(squote(mcc$Y), " is not an assay or column data from the MAE object" ))}
  mc$X <- t(as.matrix(mc$X))
  mc$Y <- as.matrix(mc$Y)

  if(!1 %in% dim(mc$Y)){ ## if Y is matrix transpose it
    mc$Y <- t(mc$Y)
  }
  return(mc)
}

##------------------------------------------
## -----  get call list including potentiall X, Y, formula, and data and retain only valid X and Y
.getXY <- function(mc){
  mc[c('data', 'formula')]<- lapply( mc[c('data', 'formula')], eval.parent)
  mc$formula <- eval.parent(mc$formula)
  mcc <- mc ## copy so can change mc but keep the call for .names2mat
  # expectedArgs <- c('X', 'Y', 'formula', 'data')
  # mc[expectedArgs] <- lapply(mc[expectedArgs], eval.parent)
  ## function to check formula's class and ensure non-NULL X and Y are not provided with it
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
      mc <- .names2mat(mc, mcc)
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
        .stop(.subclass = "inv_xy", message = "'X' and 'Y' must be assay names from 'data'")
      mc <- .names2mat(mc, mcc)
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

##------------------------------------------
## ----- function to pls methods arguments plus the methods call mode and create internal level arguments

####################################################################################
## ---------- function to handle MAE methods for pca family #TODO
# .pcaMethods <- function(ml, fun='pca'){
#   ## if assay is not valid throw appropriate error
#   arg.ind <- match(names(formals(.pca)), names(mli), 0L)
#   ml$assay <- eval.parent(ml$assay)
#   ml$X <- eval.parent(ml$X)
#   if(!ml$assay %in% tryCatch(names(assays(ml$X)), error=function(e)e)) .inv_assay()
#   ml[[1L]] <- as.name(sprintf('.%s',fun))
#   ml$X <- t(assay(ml$X, ml$assay))
#   ## evaluate the call in the parent.frame and return
#   result <- eval(ml, parent.frame())
#   return(result)
# }

.pcaMethodsHelper <- function(mc, fun='pca'){
  mc[-1L] <- lapply(mc[-1L], function(x) eval.parent(x, n = 2L))
  .check_data_assay(mc)
  mcr <- mc ## to return as output
  mcr[[1L]] <- as.name(fun) ## returned call to have 'pca' as function
  mc$X <- t(assay(mc$data, mc$X)) ## change X into matrix of assay name
  mc$data <- NULL ## remove data as not needed in internal
  mc[[1L]] <- as.name(sprintf(".%s", fun)) ## add internal to the call's function
  result <- eval.parent(mc) ## evaluate the call
  result$call <- mcr ## replace the returned call from internal by the current one
  return(result)
}

#' PCA Entry Checker
#'
#' checks X, ncomp, scale, center, max.iter, and tol. Makes adjustments to
#' them in the parent environment if necessary. for IPCA, 'center' argument does
#' not apply (is always TRUE) and NA's are not allowed.
#' Note: it possibly replaces X, ncomp, and max.iter in parent.frame using <<-.
#' Note: if the error classes are too many, you can batch replace and simplify.
#' @noRd
.pcaEntryChecker <-
  function(mc,
           check.keepX = FALSE,
           check.center = TRUE,
           check.NA = FALSE) {


  #-------- X
  if (is.data.frame(mc$X))
    mc$X <- as.matrix(mc$X)

  if (!(is.matrix(mc$X) && is.numeric(mc$X)))
    .stop(.subclass = 'inv_matrix', message = "'X' must be a numeric matrix.")

  if (any(apply(mc$X, 1, is.infinite)))
    .stop(.subclass = 'inv_matrix', message = "infinite values in 'X'.")

  ## checking NAs, only for ICA
  if (check.NA) {
    if (any(is.na(mc$X)))
      .stop(.subclass = 'inv_X',message = "missing values in 'X'.")
  }


  #-------- ncomp
  if (is.null(mc$ncomp))
    mc$ncomp <- as.integer(min(dim(mc$X)))

  if ( !is.numeric(mc$ncomp) || mc$ncomp < 1 || !is.finite(mc$ncomp))
    .stop(.subclass = 'inv_ncomp', message = "'ncomp' must be an integer not greater than min(dim(X))")

  if (mc$ncomp > min(dim(mc$X)))
    .stop(.subclass = 'inv_ncomp', message = "'ncomp' must be an integer not greater than min(dim(X))")

  #-------- keepX
  if (check.keepX) {
    if (is.null(mc$keepX))
      mc$keepX <- rep(ncol(mc$X), mc$ncomp)
    if (length(mc$keepX) != mc$ncomp)
      .stop(message = paste0("length of 'keepX' must be equal to 'ncomp', that is ", mc$ncomp, "."))
    if (any(mc$keepX > ncol(mc$X)))
      .stop(message = "each component of 'keepX' must be lower or equal to smallest dimension of 'X'")
    keepX <<- mc$keepX ## write to parent.frame
    }

  #-------- center
  ## no centering option in ICA
  if (check.center) {
    if (!is.logical(mc$center))
    {
      if (!is.numeric(mc$center) || (length(mc$center) != ncol(mc$X)))
        .stop(message = "'center' should be either a logical value or a numeric vector of length equal to the number of columns of 'X'.")
    }
  }


  #-------- scale
  if (!is.logical(mc$scale))
  {
    if (!is.numeric(mc$scale) || (length(mc$scale) != ncol(mc$X)))
      .stop(message = "'scale' should be either a logical value or a numeric vector of length equal to the number of columns of 'X'.")
  }

  #-------- max.iter
  if (is.null(mc$max.iter) || !is.numeric(mc$max.iter) || mc$max.iter < 1 || !is.finite(mc$max.iter))
    .stop("invalid value for 'max.iter'.")

  #-------- tol
  if (is.null(mc$tol) || !is.numeric(mc$tol) || mc$tol < 0 || !is.finite(mc$tol))
    .stop("invalid value for 'tol'.")

  ## write in parent.frame
  X <<- mc$X
  ncomp <<- mc$ncomp
  max.iter <<- as.integer(mc$max.iter)
}
