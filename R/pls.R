#### pls: perform a PLS
#### this function is a particular setting of .mintBlock,
#### the formatting of the input is checked in .mintWrapper
## ----------- Description ----------- 
#' Partial Least Squares (PLS) Regression
#'
#' Function to perform Partial Least Squares (PLS) regression.
#'
#' \code{pls} function fit PLS models with \eqn{1, \ldots ,}\code{ncomp}
#' components. Multi-response models are fully supported. The \code{X} and
#' \code{Y} datasets can contain missing values.
#'
#' The type of algorithm to use is specified with the \code{mode} argument.
#' Four PLS algorithms are available: PLS regression \code{("regression")}, PLS
#' canonical analysis \code{("canonical")}, redundancy analysis
#' \code{("invariant")} and the classical PLS algorithm \code{("classic")} (see
#' References). Different modes relate on how the Y matrix is deflated across
#' the iterations of the algorithms - i.e. the different components.
#'
#' - Regression mode: the Y matrix is deflated with respect to the information
#' extracted/modelled from the local regression on X. Here the goal is to
#' predict Y from X (Y and X play an asymmetric role). Consequently the latent
#' variables computed to predict Y from X are different from those computed to
#' predict X from Y.
#'
#' - Canonical mode: the Y matrix is deflated to the information
#' extracted/modelled from the local regression on Y. Here X and Y play a
#' symmetric role and the goal is similar to a Canonical Correlation type of
#' analysis.
#'
#' - Invariant mode: the Y matrix is not deflated
#'
#' - Classic mode: is similar to a regression mode. It gives identical results
#' for the variates and loadings associated to the X data set, but differences
#' for the loadings vectors associated to the Y data set (different
#' normalisations are used). Classic mode is the PLS2 model as defined by
#' Tenenhaus (1998), Chap 9.
#'
#' Note that in all cases the results are the same on the first component as
#' deflation only starts after component 1.
#'
#' The estimation of the missing values can be performed by the reconstitution
#' of the data matrix using the \code{nipals} function. Otherwise, missing
#' values are handled by casewise deletion in the \code{pls} function without
#' having to delete the rows with missing data.
#'
#' logratio transform and multilevel analysis are performed sequentially as
#' internal pre-processing step, through \code{\link{logratio.transfo}} and
#' \code{\link{withinVariation}} respectively.
#' 
## ----------- Parameters ----------- 
#' @inheritParams pca
#' @param Y Numeric vector or matrix of responses (for multi-response models),
#' or name of such an \code{assay} or
#' \code{colData} from \code{data}. \code{NA}s are allowed.
#' @param mode Character string. What type of algorithm to use, (partially)
#' matching one of \code{"regression"}, \code{"canonical"}, \code{"invariant"}
#' or \code{"classic"}. See Details.
#' @param near.zero.var Boolean, see the internal \code{\link{nearZeroVar}}
#' function (should be set to TRUE in particular for data with many zero
#' values). Setting this argument to FALSE (when appropriate) will speed up the
#' computations. Default value is FALSE
#' @param logratio One of ('none','CLR'). Default to 'none'
#' @param multilevel Design matrix for repeated measurement analysis, where
#' multlevel decomposition is required. For a one factor decomposition, the
#' repeated measures on each individual, i.e. the individuals ID is input as
#' the first column. For a 2 level factor decomposition then 2nd AND 3rd
#' columns indicate those factors. See examples in \code{?spls}).
#' @param all.outputs Boolean. Computation can be faster when some specific
#' (and non-essential) outputs are not calculated. Default = \code{TRUE}.

#'  see examples.

## ----------- Value -----------
#' @return \code{pls} returns an object of class \code{"mixo_pls"}, a list that
#' contains the following components:
#'
#' \item{X}{the centered and standardized original predictor matrix.}
#' \item{Y}{the centered and standardized original response vector or matrix.}
#' \item{ncomp}{the number of components included in the model.}
#' \item{mode}{the algorithm used to fit the model.} \item{variates}{list
#' containing the variates.} \item{loadings}{list containing the estimated
#' loadings for the \eqn{X} and \eqn{Y} variates.} \item{names}{list containing
#' the names to be used for individuals and variables.} \item{tol}{the
#' tolerance used in the iterative algorithm, used for subsequent S3 methods}
#' \item{iter}{Number of iterations of the algorthm for each component}
#' \item{max.iter}{the maximum number of iterations, used for subsequent S3
#' methods} \item{nzv}{list containing the zero- or near-zero predictors
#' information.} \item{scale}{whether scaling was applied per predictor.}
#' \item{logratio}{whether log ratio transformation for relative proportion
#' data was applied, and if so, which type of transformation.}
#' \item{explained_variance}{amount of variance explained per component (note
#' that contrary to PCA, this amount may not decrease as the aim of the method
#' is not to maximise the variance, but the covariance between data sets).}
#' \item{input.X}{numeric matrix of predictors in X that was input, before any
#' saling / logratio / multilevel transformation.} \item{mat.c}{matrix of
#' coefficients from the regression of X / residual matrices X on the
#' X-variates, to be used internally by \code{predict}.}
#' \item{defl.matrix}{residual matrices X for each dimension.}
#' 
## ----------- Ref ----------- 
#' @author Sébastien Déjean, Ignacio González, Kim-Anh Lê Cao, Al J Abadi.
#' @seealso \code{\link{spls}}, \code{\link{summary}}, \code{\link{plotIndiv}},
#' \code{\link{plotVar}}, \code{\link{predict}}, \code{\link{perf}} and
#' http://www.mixOmics.org for more details.
#' @references Tenenhaus, M. (1998). \emph{La regression PLS: theorie et
#' pratique}. Paris: Editions Technic.
#'
#' Wold H. (1966). Estimation of principal components and related models by
#' iterative least squares. In: Krishnaiah, P. R. (editors), \emph{Multivariate
#' Analysis}. Academic Press, N.Y., 391-420.
#'
#' Abdi H (2010). Partial least squares regression and projection on latent
#' structure regression (PLS Regression). \emph{Wiley Interdisciplinary
#' Reviews: Computational Statistics}, 2(1), 97-106.
#' @keywords regression multivariate
#' 
## ----------- Examples ----------- 
#' @example examples/pls-example.R
## setting the document name here so internal would not force the wrong name
#' @name pls
NULL
## ----------- Internal ----------- 
.pls = function(X=NULL,
                Y=NULL,
                ncomp = 2,
                scale = TRUE,
                mode = c("regression", "canonical", "invariant", "classic"),
                tol = 1e-06,
                max.iter = 100,
                near.zero.var = FALSE,
                logratio = "none",
                multilevel = NULL,
                all.outputs = TRUE,
                formula=NULL,
                ret.call=FALSE){
    mc <- as.list(match.call()[-1])
    
    ## make sure mode matches given arguments, and if it is not provided put as the first one in the definition
    mc$mode <- .matchArg(mode)
    ## if formula or data is given, process arguments to match default pls
    if(any(c("formula", "data") %in% names(mc))) {
        mc <- .plsMethodsHelper(mc = mc)
    }
    mc$DA <- FALSE
    # # call to '.mintWrapper'
    result <- do.call(.mintWrapper, mc)
    # choose the desired output from 'result'
    result = list(
        X = result$A[-result$indY][[1]],
        Y = result$A[result$indY][[1]],
        ncomp = result$ncomp,
        mode = result$mode,
        variates = result$variates,
        loadings = result$loadings,
        loadings.star = result$loadings.star,
        names = result$names,
        tol = result$tol,
        iter = result$iter,
        max.iter = result$max.iter,
        nzv = result$nzv,
        scale = scale,
        logratio = logratio,
        explained_variance = result$explained_variance,
        input.X = result$input.X,
        mat.c = result$mat.c#,
        #defl.matrix = result$defl.matrix
    )
    
    class(result) = c("mixo_pls")
    # output if multilevel analysis
    if (!is.null(multilevel))
    {
        result$multilevel = multilevel
        class(result) = c("mixo_mlpls",class(result))
    }
    ## keeping this for ease of re-using code
    if ( isTRUE(ret.call) ) { ## return call
        mcr <- match.call()
        mcr[-1] <- lapply(mcr[-1], eval)
        mcr[[1L]] <- quote(pls)
        result <- c(call = mcr, result)
    }
    return(invisible(result))
    
}
## ----------- Generic ----------- 
#' @param formula A \code{formula} object of form \code{LHS ~ RHS} where 
#' \code{LHS} and \code{RHS} are either matrices or assay names (no quotation)
#'  from \code{data}. 
#'  \code{X} and \code{Y} must be \code{NULL} when \code{formula} is used. 
#'  See Details.
#' @param ... Aguments passed to the generic.
#' @export
#' @rdname pls
setGeneric('pls', function(data=NULL, X=NULL, Y=NULL, formula=NULL, ...) standardGeneric('pls'))

## ----------- Methods ----------- 
#### ANY ####
#' @param ... Aguments passed to the generic.
#' @export
#' @rdname pls
setMethod('pls', 'ANY', function(data=NULL, X=NULL, Y=NULL, formula=NULL, ..., ret.call=FALSE) {
    mget(names(formals()), sys.frame(sys.nframe())) ## just to evaluate
    
    ## legacy code
    if ( class(try(data)) %in% c("data.frame", "matrix") )
        .stop(message = "mixOmics arguments have changed.
              Please carefully read the documentation and try to use named 
              arguments such as pls(X=mat, ...) as opposed to pls(mat, ...).",
              .subclass = "defunct")
    
    
    if ( !(missing(data) || is.null(data))) { ## data must be NULL
        .stop(message = "data should be a MultiAssayExperiment class, or NULL",
              .subclass = "inv_signature")
    }
    if ( !(missing(formula) || is.null(formula)) ) { ## formula must be NULL
        .stop(message = "With numerical X and Y, formula should not be provided. See ?pls",
              .subclass = "inv_signature")
    }

    mc <- match.call()
    mc[-1L] <- lapply(mc[-1L], eval.parent)
    mc$ret.call <- mc$data <- mc$formula <- NULL 
    mc[[1L]] <- quote(.pls)
    result <- eval(mc)
    .call_return(result, ret.call, mcr = match.call(), fun.name = 'pls')
})

#### signature(data = 'MultiAssayExperiment', formula = "formula") ####
## expect X and Y to be NULL
#' @export
#' @rdname pls
setMethod('pls', signature(data = 'MultiAssayExperiment', formula = 'formula'), 
          function(data=NULL, X=NULL, Y=NULL, formula=NULL, ..., ret.call=FALSE) {
              mget(names(formals()), sys.frame(sys.nframe())) ## just to evaluate
              ## X and Y NULL or missing
              if ( !((missing(X) || is.null(X)) && (missing(Y) || is.null(Y)) ) )
                  .stop(message = "Where 'data' and 'formula' are provided 'X' and 'Y' should be NULL.", 
                        .subclass = "inv_signature")
              mc <- match.call()
              mc[-1L] <- lapply(mc[-1L], eval.parent)
              .sformula_checker(mc) ## check formula validity
              mc[c('Y', 'X')] <- as.character(formula[2:3])
              mc <- .get_xy(mc = mc)
              mc$ret.call <- mc$data <- mc$formula <- NULL 
              mc[[1L]] <- quote(.pls)
              result <- eval(mc)
              .call_return(result, ret.call, mcr = match.call(), fun.name = 'pls')
          })


#### signature(data != 'MultiAssayExperiment', formula = "formula") ####
## expect X and Y to be NULL
#' @export
#' @rdname pls
setMethod('pls', signature(formula = 'formula'), 
          function(data=NULL, X=NULL, Y=NULL, formula=NULL, ..., ret.call=FALSE) {
              mget(names(formals()), sys.frame(sys.nframe())) ## just to evaluate
              if ( !(missing(data) || is.null(data))) { ## data must be NULL
                  .stop(message = "data should be a MultiAssayExperiment class, or NULL",
                        .subclass = "inv_signature")
              }
              ## X and Y NULL or missing
              if ( !((missing(X) || is.null(X)) && (missing(Y) || is.null(Y)) ) )
                  .stop(message = "'formula' must not be provided with
                        'X' and 'Y'", .subclass = "inv_signature")
              mc <- match.call()
              mc[-1L] <- lapply(mc[-1L], eval.parent)
              .sformula_checker(mc) ## check formula validity
              mc$X <- eval.parent(as.list(formula)[[3]], n = 2)
              mc$Y <- eval.parent(as.list(formula)[[2]], n = 2)
              # mc <- .get_xy(mc = mc)
              mc$ret.call <- mc$data <- mc$formula <- NULL 
              mc[[1L]] <- quote(.pls)
              result <- eval(mc)
              .call_return(result, ret.call, mcr = match.call(), fun.name = 'pls')
          })


#### signature(data = 'MultiAssayExperiment', formula != "formula") ####
## expect X and Y to be valid characters
#' @export
#' @rdname pls
setMethod('pls', signature(data = 'MultiAssayExperiment'), 
          function(data=NULL, X=NULL, Y=NULL, formula=NULL, ..., ret.call=FALSE) {
              mget(names(formals()), sys.frame(sys.nframe())) ## just to evaluate
              ## X and Y NULL or missing
              if ( !(missing(formula) || is.null(formula)) ) { ## formula must be NULL
                  .stop(message = "With numerical X and Y, formula should not be provided. See ?pls", 
                        .subclass = "inv_signature")
              }
              
              mc <- match.call()
              mc[-1L] <- lapply(mc[-1L], eval.parent)
              mc <- .get_xy(mc = mc)
              mc$ret.call <- mc$data <- mc$formula <- NULL 
              mc[[1L]] <- quote(.pls)
              result <- eval(mc)
              .call_return(result, ret.call, mcr = match.call(), fun.name = 'pls')
          })
