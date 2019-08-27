# ========================================================================== #
# block.spls: perform a horizontal PLS-DA on a combination of datasets,
#   input as a list in X
#   this function is a particular setting of .mintBlock,
#   the formatting of the input is checked in .mintWrapperBlock
# ========================================================================== #
## ----------- Description ----------- 
#' N-integration and feature selection with sparse Projection to Latent
#' Structures models (sPLS)
#'
#' Integration of multiple data sets measured on the same samples or
#' observations, with variable selection in each data set, ie. N-integration.
#' The method is partly based on Generalised Canonical Correlation Analysis.
#'
#' \code{block.spls} function fits a horizontal sPLS model with a specified
#' number of components per block). An outcome needs to be provided, either by
#' \code{Y} or by its position \code{indY} in the list of blocks \code{X}.
#' Multi (continuous)response are supported. \code{X} and \code{Y} can contain
#' missing values. Missing values are handled by being disregarded during the
#' cross product computations in the algorithm \code{block.pls} without having
#' to delete rows with missing data. Alternatively, missing data can be imputed
#' prior using the \code{nipals} function.
#'
#' The type of algorithm to use is specified with the \code{mode} argument.
#' Four PLS algorithms are available: PLS regression \code{("regression")}, PLS
#' canonical analysis \code{("canonical")}, redundancy analysis
#' \code{("invariant")} and the classical PLS algorithm \code{("classic")} (see
#' References and \code{?pls} for more details).
#'
#' Note that our method is partly based on sparse Generalised Canonical
#' Correlation Analysis and differs from the MB-PLS approaches proposed by
#' Kowalski et al., 1989, J Chemom 3(1), Westerhuis et al., 1998, J Chemom,
#' 12(5) and sparse variants Li et al., 2012, Bioinformatics 28(19); Karaman et
#' al (2014), Metabolomics, 11(2); Kawaguchi et al., 2017, Biostatistics.
#'
#' Variable selection is performed on each component for each block of
#' \code{X}, and for \code{Y} if specified, via input parameter \code{keepX}
#' and \code{keepY}.
#'
#' Note that if \code{Y} is missing and \code{indY} is provided, then variable
#' selection on \code{Y} is performed by specifying the input parameter
#' directly in \code{keepX} (no \code{keepY} is needed).
#'
## ----------- Parameters ----------- 
#' @inheritParams block.spls
#' @param keepX A list of same length as X.  Each entry is the number of
#' variables to select in each of the blocks of X for each component. By
#' default all variables are kept in the model.
#' @param keepY Only if Y is provided. Each entry is the number of variables to
#' select in each of the blocks of Y for each component.
#' 
## ----------- Value -----------
#' @return \code{block.spls} returns an object of class \code{"block.spls"}, a
#' list that contains the following components:
#'
#' \item{X}{the centered and standardized original predictor matrix.}
#' \item{indY}{the position of the outcome Y in the output list X.}
#' \item{ncomp}{the number of components included in the model for each block.}
#' \item{mode}{the algorithm used to fit the model.} \item{keepX}{Number of
#' variables used to build each component of each block} \item{keepY}{Number of
#' variables used to build each component of Y} \item{variates}{list containing
#' the variates of each block of X.} \item{loadings}{list containing the
#' estimated loadings for the variates.} \item{names}{list containing the names
#' to be used for individuals and variables.} \item{nzv}{list containing the
#' zero- or near-zero predictors information.} \item{iter}{Number of iterations
#' of the algorthm for each component} \item{explained_variance}{Percentage of
#' explained variance for each component and each block}
#' 
## ----------- Ref ----------- 
#' @author Florian Rohart, Benoit Gautier, Kim-Anh Lê Cao, Al J Abadi.
#' @seealso \code{\link{plotIndiv}}, \code{\link{plotArrow}},
#' \code{\link{plotLoadings}}, \code{\link{plotVar}}, \code{\link{predict}},
#' \code{\link{perf}}, \code{\link{selectVar}}, \code{\link{block.pls}},
#' \code{\link{block.splsda}} and http://www.mixOmics.org for more details.
#' @references Tenenhaus, M. (1998). \emph{La regression PLS: theorie et
#' pratique}. Paris: Editions Technic.
#'
#' Wold H. (1966). Estimation of principal components and related models by
#' iterative least squares. In: Krishnaiah, P. R. (editors), \emph{Multivariate
#' Analysis}. Academic Press, N.Y., 391-420.
#'
#' Tenenhaus A. and Tenenhaus M., (2011), Regularized Generalized Canonical
#' Correlation Analysis, Psychometrika, Vol. 76, Nr 2, pp 257-284.
#'
#' Tenenhaus A., Philippe C., Guillemot V, Lê Cao K.A., Grill J, Frouin V.
#' Variable selection for generalized canonical correlation analysis.
#' \emph{Biostatistics}. kxu001
#' @keywords regression multivariate
#' 
## ----------- Examples ----------- 
#' @example examples/block.spls-example.R
## setting the document name here so internal would not force the wrong name
#' @name block.spls
NULL
## ----------- Internal ----------- 
.block.spls = function(X=NULL,
                       Y=NULL,
                       indY=NULL,
                       keepX=NULL,
                       keepY=NULL,
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
    class(out) = c("block.spls","sgcca")
    
    return(invisible(out))
    
}
## ----------- Generic ----------- 
#' @export
#' @rdname block.spls
setGeneric('block.spls', function(data=NULL, X=NULL, Y=NULL, formula=NULL, ...) 
    standardGeneric('block.spls'))

## ----------- Methods ----------- 
#### ANY ####
#' @export
#' @rdname block.spls
setMethod('block.spls', 'ANY', function(data=NULL, X=NULL, Y=NULL, formula=NULL, ...) {
    
    tryCatch(mget(names(formals()), sys.frame(sys.nframe())),
             error = function(e) stop(e$message, call. = FALSE))
    mc <- match.call()
    mc[-1L] <- lapply(mc[-1L], eval)
    
    ## check signature format
    .check_sig_ANY(mc, fun = "block.spls")
    
    mc[[1L]] <- quote(.block.spls)
    result <- eval(mc)
    
    .call_return(result, match.call(), fun.name = 'block.spls')
})

#### signature(data = 'MultiAssayExperiment', formula != "formula") ####
#' @export
#' @rdname block.spls
setMethod('block.spls', signature(data = 'MultiAssayExperiment'), 
          function(data=NULL, X=NULL, Y=NULL, formula=NULL, ...) {
              tryCatch(mget(names(formals()), sys.frame(sys.nframe())),
                       error = function(e) stop(e$message, call. = FALSE))
              mc <- match.call()
              mc[-1] <- lapply(mc[-1], eval)
              mc <- .block_get_xy(mc = mc)
              mc$data <- mc$formula <- NULL 
              mc[[1L]] <- quote(.block.spls)
              result <- eval(mc)
              .call_return(result, mc$ret.call, mcr = match.call(), fun.name = 'block.spls')
          })


#### signature(data != 'MultiAssayExperiment', formula = "formula") ####
#' @export
#' @rdname block.spls
setMethod('block.spls', signature(formula = 'formula'), 
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
              mc[[1L]] <- quote(.block.spls)
              result <- eval(mc)
              .call_return(result, mc$ret.call, mcr = match.call(), fun.name = 'block.spls')
          })