# ========================================================================================================
# plsda: perform a PLS-DA
# this function is a particular setting of .mintBlock, the formatting of the input is checked in .mintWrapper
# ========================================================================================================
## ----------- Description ----------- 
#' Partial Least Squares Discriminant Analysis (PLS-DA).
#'
#' Function to perform standard Partial Least Squares regression to classify
#' samples.
#'
#' \code{plsda} function fit PLS models with \eqn{1,...,}\code{ncomp}
#' components to the factor or class vector \code{Y}. The appropriate indicator
#' matrix is created.
#'
#' logratio transform and multilevel analysis are performed sequentially as
#' internal pre-processing step, through \code{\link{logratio.transfo}} and
#' \code{\link{withinVariation}} respectively.
#'
#' Logratio can only be applied if the data do not contain any 0 value (for
#' count data, we thus advise the normalise raw data with a 1 offset).
#'
#' More details about the PLS modes in \code{?pls}.
#' 
## ----------- Parameters ----------- 
#' @inheritParams pls
#' @param Y a factor or a class vector for the discrete outcome.
#' @param multilevel sample information for multilevel decomposition for
#' repeated measurements. A numeric matrix or data frame indicating the
#' repeated measures on each individual, i.e. the individuals ID. See examples
#' in \code{?splsda}.
#'  
## ----------- Value -----------
#' @return \code{plsda} returns an object of class \code{"plsda"}, a list that
#' contains the following components:
#'
#' \item{X}{the centered and standardized original predictor matrix.}
#' \item{Y}{the centered and standardized indicator response vector or matrix.}
#' \item{ind.mat}{the indicator matrix.} \item{ncomp}{the number of components
#' included in the model.} \item{variates}{list containing the \code{X} and
#' \code{Y} variates.} \item{loadings}{list containing the estimated loadings
#' for the variates.} \item{names}{list containing the names to be used for
#' individuals and variables.} \item{nzv}{list containing the zero- or
#' near-zero predictors information.} \item{tol}{the tolerance used in the
#' iterative algorithm, used for subsequent S3 methods} \item{max.iter}{the
#' maximum number of iterations, used for subsequent S3 methods}
#' \item{iter}{Number of iterations of the algorthm for each component}
#' \item{explained_variance}{amount of variance explained per component (note
#' that contrary to PCA, this amount may not decrease as the aim of the method
#' is not to maximise the variance, but the covariance between X and the dummy
#' matrix Y).} \item{mat.c}{matrix of coefficients from the regression of X /
#' residual matrices X on the X-variates, to be used internally by
#' \code{predict}.} \item{defl.matrix}{residual matrices X for each dimension.}
#' 
## ----------- Ref ----------- 
#' @author Ignacio González, Kim-Anh Lê Cao, Al J Abadi.
#' @seealso \code{\link{splsda}}, \code{\link{summary}},
#' \code{\link{plotIndiv}}, \code{\link{plotVar}}, \code{\link{predict}},
#' \code{\link{perf}}, \code{\link{mint.block.plsda}},
#' \code{\link{block.plsda}} and http://mixOmics.org for more details.
#' @references On PLSDA: Barker M and Rayens W (2003). Partial least squares
#' for discrimination. \emph{Journal of Chemometrics} \bold{17}(3), 166-173.
#' Perez-Enciso, M. and Tenenhaus, M. (2003). Prediction of clinical outcome
#' with microarray data: a partial least squares discriminant analysis (PLS-DA)
#' approach. \emph{Human Genetics} \bold{112}, 581-592. Nguyen, D. V. and
#' Rocke, D. M. (2002). Tumor classification by partial least squares using
#' microarray gene expression data. \emph{Bioinformatics} \bold{18}, 39-50. On
#' log ratio transformation: Filzmoser, P., Hron, K., Reimann, C.: Principal
#' component analysis for compositional data with outliers. Environmetrics
#' 20(6), 621-632 (2009) Lê Cao K.-A., Costello ME, Lakis VA, Bartolo, F,Chua
#' XY, Brazeilles R, Rondeau P. MixMC: Multivariate insights into Microbial
#' Communities. PLoS ONE, 11(8): e0160169 (2016). On multilevel decomposition:
#' Westerhuis, J.A., van Velzen, E.J., Hoefsloot, H.C., Smilde, A.K.:
#' Multivariate paired data analysis: multilevel plsda versus oplsda.
#' Metabolomics 6(1), 119-128 (2010) Liquet, B., Lê Cao K.-A., Hocini, H.,
#' Thiebaut, R.: A novel approach for biomarker selection and the integration
#' of repeated measures experiments from two assays. BMC bioinformatics 13(1),
#' 325 (2012)
#' @keywords regression multivariate
#' 
## ----------- Examples ----------- 
#' @example examples/plsda-example.R
## setting the document name here so internal would not force the wrong name
#' @name plsda
NULL
## ----------- Internal ----------- 
.plsda = function(X=NULL,
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
                ret.call=FALSE){
  mc <- match.call() ## get the call args to adjust and pass to internal
  #### from previous version {
  #-- validation des arguments --#
  # most of the checks are done in the wrapper.mint.spls.hybrid function
  if (is.null(multilevel))
  {
    if (is.null(Y))
      stop("'Y' has to be something else than NULL.")
    
    if (is.null(dim(Y)))
    {
      Y = factor(Y)
    } else {
      stop("'Y' should be a factor or a class vector.")
    }
    
    if (nlevels(Y) == 1)
      stop("'Y' should be a factor with more than one level")
    
    Y.mat = unmap(Y)
    colnames(Y.mat) = levels(Y)
    
  } else {
    # we expect a vector or a 2-columns matrix in 'Y' and the repeated measurements in 'multilevel'
    multilevel = data.frame(multilevel)
    
    if ((nrow(X) != nrow(multilevel)))
      stop("unequal number of rows in 'X' and 'multilevel'.")
    
    if (ncol(multilevel) != 1)
      stop("'multilevel' should have a single column for the repeated measurements, other factors should be included in 'Y'.")
    
    if (!is.null(ncol(Y)) && !ncol(Y) %in% c(0,1,2))# multilevel 1 or 2 factors
      stop("'Y' should either be a factor, a single column data.frame containing a factor, or a 2-columns data.frame containing 2 factors.")
    multilevel = data.frame(multilevel, Y)
    multilevel[, 1] = as.numeric(factor(multilevel[, 1])) # we want numbers for the repeated measurements
    
    Y.mat = NULL
  }
  #### }
  mc$Y <- Y.mat ## matrix form of Y
  ## make sure mode matches given arguments, and if it is not provided put as the first one in the definition
  mc$mode <- .matchArg(mode)
  mc$DA <- TRUE
  mc$ret.call <- NULL ## not need by wrapper
  # # call to '.mintWrapper'
  mc[[1L]] <- quote(.mintWrapper)
  result <- eval(mc)
  # choose the desired output from 'result'
  result = structure(list(
    X = result$A[-result$indY][[1]],
    Y = if (is.null(multilevel))
    {
      Y
    } else {
      result$Y.factor
    },
    ind.mat = result$A[result$indY][[1]],
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
    explained_variance = result$explained_variance,#[-result$indY],
    input.X = result$input.X,
    mat.c = result$mat.c#,
  ), class = c("mixo_plsda","mixo_pls","DA"))
  
  # output if multilevel analysis
  if (!is.null(multilevel))
  {
    result$multilevel = multilevel
    class(result) = c("mixo_mlplsda",class(result))
  }
  ## keeping this cause we might revert S4 methods
  ## the return us handled in methods now
  ret.call <- FALSE
  if ( isTRUE(ret.call) ) { ## return call
    mcr <- match.call()
    mcr[-1] <- lapply(mcr[-1], eval)
    mcr[[1L]] <- quote(plsda)
    result <- c(call = mcr, result)
  }
  return(invisible(result))
  
}
## ----------- Generic ----------- 
#' @export
#' @rdname plsda
setGeneric('plsda', function(data=NULL, X=NULL, Y=NULL, formula=NULL, ...) standardGeneric('plsda'))

## ----------- Methods ----------- 
#### ANY ####
#' @export
#' @rdname plsda
setMethod('plsda', 'ANY', function(data=NULL, X=NULL, Y=NULL, formula=NULL, ...) {
  mget(names(formals()), sys.frame(sys.nframe())) ## just to evaluate
  
  ## legacy code
  if ( class(try(data)) %in% c("data.frame", "matrix") )
    .stop(message = "mixOmics arguments have changed.
              Please carefully read the documentation and try to use named 
              arguments such as plsda(X=mat, ...) as opposed to plsda(mat, ...).",
          .subclass = "defunct")
  
  if ( !(missing(data) || is.null(data))) { ## data must be NULL
    .stop(message = "data should be a MultiAssayExperiment class, or NULL",
          .subclass = "inv_signature")
  }
  if ( !(missing(formula) || is.null(formula)) ) { ## formula must be NULL
    .stop(message = "With numerical X and Y, formula should not be provided. 
              See ?plsda",
          .subclass = "inv_signature")
  }
  
  mc <- match.call()
  mc[-1L] <- lapply(mc[-1L], eval.parent)
  mc$data <- mc$formula <- NULL 
  mc[[1L]] <- quote(.plsda)
  result <- eval(mc)
  .call_return(result, mc$ret.call, mcr = match.call(), fun.name = 'plsda')
})

#### signature(data = 'MultiAssayExperiment', formula = "formula") ####
## expect X and Y to be NULL
#' @export
#' @rdname plsda
setMethod('plsda', signature(data = 'MultiAssayExperiment', formula = 'formula'), 
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
            mc[[1L]] <- quote(.plsda)
            result <- eval(mc)
            .call_return(result, mc$ret.call, mcr = match.call(), fun.name = 'plsda')
          })


#### signature(data != 'MultiAssayExperiment', formula = "formula") ####
## expect X and Y to be NULL
#' @export
#' @rdname plsda
setMethod('plsda', signature(formula = 'formula'), 
          function(data=NULL, X=NULL, Y=NULL, formula=NULL, ...) {
            mget(names(formals()), sys.frame(sys.nframe())) ## just to evaluate
            mc <- match.call()
            mc[-1L] <- lapply(mc[-1L], eval.parent)
            .sformula_checker(mc) ## check formula validity
            mc$X <- eval.parent(as.list(formula)[[3]], n = 2)
            mc$Y <- eval.parent(as.list(formula)[[2]], n = 2)
            # mc <- .get_xy(mc = mc)
            mc$data <- mc$formula <- NULL 
            mc[[1L]] <- quote(.plsda)
            result <- eval(mc)
            .call_return(result, mc$ret.call, mcr = match.call(), fun.name = 'plsda')
          })


#### signature(data = 'MultiAssayExperiment', formula != "formula") ####
## expect X and Y to be valid characters
#' @export
#' @rdname plsda
setMethod('plsda', signature(data = 'MultiAssayExperiment'), 
          function(data=NULL, X=NULL, Y=NULL, formula=NULL, ...) {
            mget(names(formals()), sys.frame(sys.nframe())) ## just to evaluate
            ## X and Y NULL or missing
            if ( !(missing(formula) || is.null(formula)) ) { ## formula must be NULL
              .stop(message = "With numerical X and Y, formula should not be provided. See ?plsda", 
                    .subclass = "inv_signature")
            }
            mc <- match.call()
            mc[-1L] <- lapply(mc[-1L], eval.parent)
            mc <- .get_xy(mc = mc)
            mc$data <- mc$formula <- NULL 
            mc[[1L]] <- quote(.plsda)
            result <- eval(mc)
            .call_return(result, mc$ret.call, mcr = match.call(), fun.name = 'plsda')
          })
