# ========================================================================== #
# block.pls: perform a horizontal PLS on a combination of datasets,
# input as a list in X
# this function is a particular setting of .mintBlock,
# the formatting of the input is checked in .mintWrapperBlock
# ========================================================================== #
## ----------- Description ----------- 
#' N-integration with Projection to Latent Structures models (PLS)
#'
#' Integration of multiple data sets measured on the same samples or
#' observations, ie. N-integration. The method is partly based on Generalised
#' Canonical Correlation Analysis.
#'
#' \code{block.pls} function fits a horizontal integration PLS model with a
#' specified number of components per block). An outcome needs to be provided,
#' either by \code{Y} or by its position \code{indY} in the list of blocks
#' \code{X}. Multi (continuous)response are supported. \code{X} and \code{Y}
#' can contain missing values. Missing values are handled by being disregarded
#' during the cross product computations in the algorithm \code{block.pls}
#' without having to delete rows with missing data. Alternatively, missing data
#' can be imputed prior using the \code{nipals} function.
#'
#' The type of algorithm to use is specified with the \code{mode} argument.
#' Four PLS algorithms are available: PLS regression \code{("regression")}, PLS
#' canonical analysis \code{("canonical")}, redundancy analysis
#' \code{("invariant")} and the classical PLS algorithm \code{("classic")} (see
#' References and \code{?pls} for more details).
#'
#' Note that our method is partly based on Generalised Canonical Correlation
#' Analysis and differs from the MB-PLS approaches proposed by Kowalski et al.,
#' 1989, J Chemom 3(1) and Westerhuis et al., 1998, J Chemom, 12(5).
#'
## ----------- Parameters ----------- 
#' @inheritParams pls
#' @param X A list of data sets (called 'blocks') measured on the same samples.
#' Data in the list should be arranged in matrices, samples x variables, with
#' samples order matching in all data sets.
#' @param Y Matrix response for a multivariate regression framework. Data
#' should be continuous variables (see block.plsda for supervised
#' classification and factor reponse)
#' @param indY To supply if Y is missing, indicates the position of the matrix
#' response in the list \code{X}
#' @param ncomp the number of components to include in the model. Default to 2.
#' Applies to all blocks.
#' @param design numeric matrix of size (number of blocks in X) x (number of
#' blocks in X) with values between 0 and 1. Each value indicates the strenght
#' of the relationship to be modelled between two blocks; a value of 0
#' indicates no relationship, 1 is the maximum value. If \code{Y} is provided
#' instead of \code{indY}, the \code{design} matrix is changed to include
#' relationships to \code{Y}.
#' @param scheme Either "horst", "factorial" or "centroid". Default =
#' \code{horst}, see reference.
#' @param init Mode of initialization use in the algorithm, either by Singular
#' Value Decompostion of the product of each block of X with Y ("svd") or each
#' block independently ("svd.single"). Default = \code{svd.single}.
#' 
## ----------- Value -----------
#' @return \code{spls} returns an object of class \code{"spls"}, a list that
#' contains the following components:
#'
#' @return \code{block.pls} returns an object of class \code{"block.pls"}, a
#' list that contains the following components:
#'
#' \item{X}{the centered and standardized original predictor matrix.}
#' \item{indY}{the position of the outcome Y in the output list X.}
#' \item{ncomp}{the number of components included in the model for each block.}
#' \item{mode}{the algorithm used to fit the model.} \item{variates}{list
#' containing the variates of each block of X.} \item{loadings}{list containing
#' the estimated loadings for the variates.} \item{names}{list containing the
#' names to be used for individuals and variables.} \item{nzv}{list containing
#' the zero- or near-zero predictors information.} \item{iter}{Number of
#' iterations of the algorthm for each component}
#' \item{explained_variance}{Percentage of explained variance for each
#' component and each block}
#' 
## ----------- Ref ----------- 
#' @author Florian Rohart, Benoit Gautier, Kim-Anh Lê Cao, Al J Abadi.
#' @seealso \code{\link{plotIndiv}}, \code{\link{plotArrow}},
#' \code{\link{plotLoadings}}, \code{\link{plotVar}}, \code{\link{predict}},
#' \code{\link{perf}}, \code{\link{selectVar}}, \code{\link{block.spls}},
#' \code{\link{block.plsda}} and http://www.mixOmics.org for more details.
#' @references Tenenhaus, M. (1998). \emph{La regression PLS: theorie et
#' pratique}. Paris: Editions Technic.
#'
#' Wold H. (1966). Estimation of principal components and related models by
#' iterative least squares. In: Krishnaiah, P. R. (editors), \emph{Multivariate
#' Analysis}. Academic Press, N.Y., 391-420.
#'
#' Tenenhaus A. and Tenenhaus M., (2011), Regularized Generalized Canonical
#' Correlation Analysis, Psychometrika, Vol. 76, Nr 2, pp 257-284.
#' @keywords regression multivariate
#' 
## ----------- Examples ----------- 
#' @example examples/block.pls-example.R
## setting the document name here so internal would not force the wrong name
#' @name block.pls
NULL
## ----------- Internal ----------- 
.block.pls = function(X=NULL,
                      Y=NULL,
                      indY=NULL,
                      ncomp = 2,
                      design=NULL,
                      scheme=c("horst", "factorial", "centroid"),
                      mode = c("regression", "canonical", "invariant", "classic"),
                      scale = TRUE,
                      init=c("svd.single", "svd") ,
                      tol = 1e-06,
                      max.iter = 100,
                      near.zero.var = FALSE,
                      all.outputs = TRUE,
                      ret.call=FALSE)
{
    mc <- match.call.defaults()
    mc$scheme <- .matchArg(scheme)
    mc$mode <- .matchArg(mode)
    mc$init <- .matchArg(init)
    mc$ret.call <- NULL ## drop as it is not used by internal wrapper
    # call to '.mintWrapperBlock'
    mc[[1L]] <- quote(.mintWrapperBlock)
    result <- eval(mc)
    
    # calculate weights for each dataset
    weights = .getWeights(result$variates, indY = result$indY)
    
    # choose the desired output from 'result'
    out <- list(X = result$A,
             indY = result$indY,
             ncomp = result$ncomp,
             mode = result$mode,
             variates = result$variates,
             loadings = result$loadings,
             crit = result$crit,
             AVE = result$AVE,
             names = result$names,
             init = result$init,
             tol = result$tol,
             iter = result$iter,
             max.iter = result$max.iter,
             nzv = result$nzv,
             scale = result$scale,
             design = result$design,
             scheme = result$scheme,
             weights = weights,
             explained_variance = result$explained_variance)
    
    # give a class
    class(out) = c("block.pls","sgcca")
    
    return(invisible(out))
    
}
## ----------- Generic ----------- 
#' @export
#' @rdname block.pls
setGeneric('block.pls', function(data=NULL, X=NULL, Y=NULL, formula=NULL, ...) 
    standardGeneric('block.pls'))

## ----------- Methods ----------- 
#### ANY ####
#' @export
#' @rdname block.pls
setMethod('block.pls', 'ANY', function(data=NULL, X=NULL, Y=NULL, formula=NULL, ...) {
    
    tryCatch(mget(names(formals()), sys.frame(sys.nframe())),
                   error = function(e) stop(e$message, call. = FALSE))
    mc <- match.call()
    mc[-1L] <- lapply(mc[-1L], eval)
    
    ## check signature format
    .check_sig_ANY(mc, fun = "block.pls")
    
    mc[[1L]] <- quote(.block.pls)
    result <- eval(mc)
    
    .call_return(result, match.call(), fun.name = 'block.pls')
})
#### signature(data = 'MultiAssayExperiment', formula != "formula") ####
#' @export
#' @rdname block.pls
setMethod('block.pls', signature(data = 'MultiAssayExperiment'), 
          function(data=NULL, X=NULL, Y=NULL, formula=NULL, ...) {
              tryCatch(mget(names(formals()), sys.frame(sys.nframe())),
                       error = function(e) stop(e$message, call. = FALSE))
              mc <- match.call()
              mc[-1] <- lapply(mc[-1], eval)
              mc <- .block_get_xy(mc = mc)
              mc$data <- mc$formula <- NULL 
              mc[[1L]] <- quote(.block.pls)
              result <- eval(mc)
              .call_return(result, mc$ret.call, mcr = match.call(), fun.name = 'block.pls')
          })


#### signature(data != 'MultiAssayExperiment', formula = "formula") ####
#' @export
#' @rdname block.pls
setMethod('block.pls', signature(formula = 'formula'), 
          function(data=NULL, X=NULL, Y=NULL, formula=NULL, ...) {
              tryCatch(mget(names(formals()), sys.frame(sys.nframe())),
                       error = function(e) stop(e$message, call. = FALSE))
              mc <- match.call()
              mc[-1L] <- lapply(mc[-1L], eval)
              .formula_checker(mc, block = TRUE) ## check formula validity
              mf <- stats::model.frame(mc$formula) ## THANK YOU stats::model.frame *cries*
              mc$Y <- as.matrix(mf[[1]])
              mc$X <- as.list(mf[-1])
              mc$data <- mc$formula <- NULL 
              mc[[1L]] <- quote(.block.pls)
              result <- eval(mc)
              .call_return(result, mc$ret.call, mcr = match.call(), fun.name = 'block.pls')
          })