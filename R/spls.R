# ========================================================================================================
# spls: perform a sparse PLS
# this function is a particular setting of .mintBlock, the formatting of the input is checked in .mintWrapper
# ========================================================================================================
## ----------- Description ----------- 
#' Sparse Partial Least Squares (sPLS)
#'
#' Function to perform sparse Partial Least Squares (sPLS). The sPLS approach
#' combines both integration and variable selection simultaneously on two data
#' sets in a one-step strategy.
#'
#' \code{spls} function fit sPLS models with \eqn{1, \ldots ,}\code{ncomp}
#' components. Multi-response models are fully supported. The \code{X} and
#' \code{Y} datasets can contain missing values.
#'
#' The type of algorithm to use is specified with the \code{mode} argument.
#' Four PLS algorithms are available: PLS regression \code{("regression")}, PLS
#' canonical analysis \code{("canonical")}, redundancy analysis
#' \code{("invariant")} and the classical PLS algorithm \code{("classic")} (see
#' References and \code{?pls} for more details).
#'
#' The estimation of the missing values can be performed by the reconstitution
#' of the data matrix using the \code{nipals} function. Otherwise, missing
#' values are handled by casewise deletion in the \code{spls} function without
#' having to delete the rows with missing data.
#'
#' logratio transform and multilevel analysis are performed sequentially as
#' internal pre-processing step, through \code{\link{logratio.transfo}} and
#' \code{\link{withinVariation}} respectively.
#'
#' Multilevel sPLS enables the integration of data measured on two different
#' data sets on the same individuals. This approach differs from multilevel
#' sPLS-DA as the aim is to select subsets of variables from both data sets
#' that are highly positively or negatively correlated across samples. The
#' approach is unsupervised, i.e. no prior knowledge about the sample groups is
#' included.
## ----------- Parameters ----------- 
#' @inheritParams pls
#' @param keepX Numeric vector of length \code{ncomp}, the number of variables
#' to keep in \eqn{X}-loadings. By default all variables are kept in the model.
#' @param keepY Numeric vector of length \code{ncomp}, the number of variables
#' to keep in \eqn{Y}-loadings. By default all variables are kept in the model.
#'  
## ----------- Value -----------
#' @return \code{spls} returns an object of class \code{"spls"}, a list that
#' contains the following components:
#'
#' \item{X}{the centered and standardized original predictor matrix.}
#' \item{Y}{the centered and standardized original response vector or matrix.}
#' \item{ncomp}{the number of components included in the model.}
#' \item{mode}{the algorithm used to fit the model.} \item{keepX}{number of
#' \eqn{X} variables kept in the model on each component.} \item{keepY}{number
#' of \eqn{Y} variables kept in the model on each component.}
#' \item{variates}{list containing the variates.} \item{loadings}{list
#' containing the estimated loadings for the \eqn{X} and \eqn{Y} variates.}
#' \item{names}{list containing the names to be used for individuals and
#' variables.} \item{tol}{the tolerance used in the iterative algorithm, used
#' for subsequent S3 methods} \item{iter}{Number of iterations of the algorthm
#' for each component} \item{max.iter}{the maximum number of iterations, used
#' for subsequent S3 methods} \item{nzv}{list containing the zero- or near-zero
#' predictors information.} \item{scale}{whether scaling was applied per
#' predictor.} \item{logratio}{whether log ratio transformation for relative
#' proportion data was applied, and if so, which type of transformation.}
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
#' @seealso \code{\link{pls}}, \code{\link{summary}}, \code{\link{plotIndiv}},
#' \code{\link{plotVar}}, \code{\link{cim}}, \code{\link{network}},
#' \code{\link{predict}}, \code{\link{perf}} and http://www.mixOmics.org for
#' more details.
#' @references Sparse PLS: canonical and regression modes:
#'
#' Lê Cao, K.-A., Martin, P.G.P., Robert-Granie, C. and Besse, P. (2009).
#' Sparse canonical methods for biological data integration: application to a
#' cross-platform study. \emph{BMC Bioinformatics} \bold{10}:34.
#'
#' Lê Cao, K.-A., Rossouw, D., Robert-Granie, C. and Besse, P. (2008). A sparse
#' PLS for variable selection when integrating Omics data. \emph{Statistical
#' Applications in Genetics and Molecular Biology} \bold{7}, article 35.
#'
#' Sparse SVD: Shen, H. and Huang, J. Z. (2008). Sparse principal component
#' analysis via regularized low rank matrix approximation. \emph{Journal of
#' Multivariate Analysis} \bold{99}, 1015-1034.
#'
#' PLS methods: Tenenhaus, M. (1998). \emph{La regression PLS: theorie et
#' pratique}. Paris: Editions Technic. Chapters 9 and 11.
#'
#' Abdi H (2010). Partial least squares regression and projection on latent
#' structure regression (PLS Regression). \emph{Wiley Interdisciplinary
#' Reviews: Computational Statistics}, 2(1), 97-106.
#'
#' Wold H. (1966). Estimation of principal components and related models by
#' iterative least squares. In: Krishnaiah, P. R. (editors), \emph{Multivariate
#' Analysis}. Academic Press, N.Y., 391-420.
#'
#' On multilevel analysis:
#'
#' Liquet, B., Lê Cao, K.-A., Hocini, H. and Thiebaut, R. (2012) A novel
#' approach for biomarker selection and the integration of repeated measures
#' experiments from two platforms. \emph{BMC Bioinformatics} \bold{13}:325.
#'
#' Westerhuis, J. A., van Velzen, E. J., Hoefsloot, H. C., and Smilde, A. K.
#' (2010). Multivariate paired data analysis: multilevel PLSDA versus OPLSDA.
#' \emph{Metabolomics}, \bold{6}(1), 119-128.
#' @keywords regression multivariate
#' 
## ----------- Examples ----------- 
#' @example examples/spls-example.R
## setting the document name here so internal would not force the wrong name
#' @name spls
NULL
## ----------- Internal ----------- 
.spls = function(X=NULL,
                 Y=NULL,
                 ncomp = 2,
                 mode = c("regression", "canonical", "invariant", "classic"),
                 keepX=NULL,
                 keepY=NULL,
                 scale = TRUE,
                 tol = 1e-06,
                 max.iter = 100,
                 near.zero.var = FALSE,
                 logratio = "none",   # one of "none", "CLR"
                 multilevel = NULL,
                 all.outputs = TRUE,
                 data=NULL,
                 formula=NULL)
{
    mc <- as.list(match.call()[-1])
    
    ## make sure mode matches given arguments, and if it is not provided put as the first one in the definition
    mc$mode <- .matchArg(mode)
    ## if formula or data is given, process arguments to match default pls
    if(any(c("formula", "data") %in% names(mc))){
        mc <- .plsMethodsHelper(mc=mc)
    }
    mc$DA <- FALSE
    # # call to '.mintWrapper'
    result <- do.call(.mintWrapper, mc)
    # choose the desired output from 'result'
    out = list(
        call = match.call(),
        X = result$A[-result$indY][[1]],
        Y = result$A[result$indY][[1]],
        ncomp = result$ncomp,
        mode = result$mode,
        keepX = result$keepX,
        keepY = result$keepY,
        variates = result$variates,
        loadings = result$loadings,
        loadings.star = result$loadings.star,
        names = result$names,
        tol = result$tol,iter = result$iter,
        max.iter = result$max.iter,
        nzv = result$nzv,
        scale = scale,
        logratio = logratio,
        explained_variance = result$explained_variance,
        input.X = result$input.X,
        mat.c = result$mat.c
    )
    
    
    class(out) = c("mixo_spls")
    # output if multilevel analysis
    if (!is.null(multilevel))
    {
        out$multilevel = multilevel
        class(out) = c("mixo_mlspls",class(out))
    }
    
    return(invisible(out))
}
## ----------- Generic ----------- 
#' @export
#' @rdname spls
setGeneric('spls', function(data=NULL, X=NULL, Y=NULL, formula=NULL, ...) standardGeneric('spls'))

## ----------- Methods ----------- 
#### ANY ####
#' @export
#' @rdname spls
setMethod('spls', 'ANY', function(data=NULL, X=NULL, Y=NULL, formula=NULL, ...) {
    mget(names(formals()), sys.frame(sys.nframe())) ## just to evaluate
    
    ## legacy code
    if ( class(try(data)) %in% c("data.frame", "matrix") )
        .stop(message = "mixOmics arguments have changed.
              Please carefully read the documentation and try to use named 
              arguments such as spls(X=mat, ...) as opposed to spls(mat, ...).",
              .subclass = "defunct")
    
    if ( !(missing(data) || is.null(data))) { ## data must be NULL
        .stop(message = "data should be a MultiAssayExperiment class, or NULL",
              .subclass = "inv_signature")
    }
    if ( !(missing(formula) || is.null(formula)) ) { ## formula must be NULL
        .stop(message = "With numerical X and Y, formula should not be provided. 
              See ?spls",
              .subclass = "inv_signature")
    }
    
    mc <- match.call()
    mc[-1L] <- lapply(mc[-1L], eval.parent)
    mc$data <- mc$formula <- NULL 
    mc[[1L]] <- quote(.spls)
    result <- eval(mc)
    .call_return(result, mc$ret.call, mcr = match.call(), fun.name = 'spls')
})

#### signature(data = 'MultiAssayExperiment', formula = "formula") ####
## expect X and Y to be NULL
#' @export
#' @rdname spls
setMethod('spls', signature(data = 'MultiAssayExperiment', formula = 'formula'), 
          function(data=NULL, X=NULL, Y=NULL, formula=NULL, ...) {
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
              mc$data <- mc$formula <- NULL 
              mc[[1L]] <- quote(.spls)
              result <- eval(mc)
              .call_return(result, mc$ret.call, mcr = match.call(), fun.name = 'spls')
          })


#### signature(data != 'MultiAssayExperiment', formula = "formula") ####
## expect X and Y to be NULL
#' @export
#' @rdname spls
setMethod('spls', signature(formula = 'formula'), 
          function(data=NULL, X=NULL, Y=NULL, formula=NULL, ...) {
              mget(names(formals()), sys.frame(sys.nframe())) ## just to evaluate
              mc <- match.call()
              mc[-1L] <- lapply(mc[-1L], eval.parent)
              .sformula_checker(mc) ## check formula validity
              mc$X <- eval.parent(as.list(formula)[[3]], n = 2)
              mc$Y <- eval.parent(as.list(formula)[[2]], n = 2)
              # mc <- .get_xy(mc = mc)
              mc$data <- mc$formula <- NULL 
              mc[[1L]] <- quote(.spls)
              result <- eval(mc)
              .call_return(result, mc$ret.call, mcr = match.call(), fun.name = 'spls')
          })


#### signature(data = 'MultiAssayExperiment', formula != "formula") ####
## expect X and Y to be valid characters
#' @export
#' @rdname spls
setMethod('spls', signature(data = 'MultiAssayExperiment'), 
          function(data=NULL, X=NULL, Y=NULL, formula=NULL, ...) {
              mget(names(formals()), sys.frame(sys.nframe())) ## just to evaluate
              ## X and Y NULL or missing
              if ( !(missing(formula) || is.null(formula)) ) { ## formula must be NULL
                  .stop(message = "With numerical X and Y, formula should not be provided. See ?spls", 
                        .subclass = "inv_signature")
              }
              
              mc <- match.call()
              mc[-1L] <- lapply(mc[-1L], eval.parent)
              mc <- .get_xy(mc = mc)
              mc$data <- mc$formula <- NULL 
              mc[[1L]] <- quote(.spls)
              result <- eval(mc)
              .call_return(result, mc$ret.call, mcr = match.call(), fun.name = 'spls')
          })
