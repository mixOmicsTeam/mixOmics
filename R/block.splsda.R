# ========================================================================== #
# block.splsda: perform a horizontal sPLS-DA on a combination of datasets,
#   input as a list in X
#   this function is a particular setting of .mintBlock,
#   the formatting of the input is checked in .mintWrapperBlock
# ========================================================================== #
## ----------- Description ----------- 
#' N-integration and feature selection with Projection to Latent Structures
#' models (PLS) with sparse Discriminant Analysis
#'
#' Integration of multiple data sets measured on the same samples or
#' observations to classify a discrete outcome to classify a discrete outcome
#' and select features from each data set, ie. N-integration with sparse
#' Discriminant Analysis. The method is partly based on Generalised Canonical
#' Correlation Analysis.
#'
#'
#' \code{block.splsda} function fits a horizontal integration PLS-DA model with
#' a specified number of components per block). A factor indicating the
#' discrete outcome needs to be provided, either by \code{Y} or by its position
#' \code{indY} in the list of blocks \code{X}.
#'
#' \code{X} can contain missing values. Missing values are handled by being
#' disregarded during the cross product computations in the algorithm
#' \code{block.pls} without having to delete rows with missing data.
#' Alternatively, missing data can be imputed prior using the \code{nipals}
#' function.
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
#' Variable selection is performed on each component for each block of \code{X}
#' if specified, via input parameter \code{keepX}.
#'
## ----------- Parameters ----------- 
#' @inheritParams block.plsda
#' @param keepX A list of same length as X.  Each entry is the number of
#' variables to select in each of the blocks of X for each component. By
#' default all variables are kept in the model.
#' 
## ----------- Value -----------
#' @return \code{block.splsda} returns an object of class \code{"block.splsda",
#' "block.spls"}, a list that contains the following components:
#'
#' \item{X}{the centered and standardized original predictor matrix.}
#' \item{indY}{the position of the outcome Y in the output list X.}
#' \item{ncomp}{the number of components included in the model for each block.}
#' \item{mode}{the algorithm used to fit the model.} \item{keepX}{Number of
#' variables used to build each component of each block} \item{variates}{list
#' containing the variates of each block of X.} \item{loadings}{list containing
#' the estimated loadings for the variates.} \item{names}{list containing the
#' names to be used for individuals and variables.} \item{nzv}{list containing
#' the zero- or near-zero predictors information.} \item{iter}{Number of
#' iterations of the algorthm for each component} \item{weights}{Correlation
#' between the variate of each block and the variate of the outcome. Used to
#' weight predictions.} \item{explained_variance}{Percentage of explained
#' variance for each component and each block}
#' 
## ----------- Ref ----------- 
#' @author Florian Rohart, Benoit Gautier, Kim-Anh Lê Cao, Al J Abadi.
#' @seealso \code{\link{plotIndiv}}, \code{\link{plotArrow}},
#' \code{\link{plotLoadings}}, \code{\link{plotVar}}, \code{\link{predict}},
#' \code{\link{perf}}, \code{\link{selectVar}}, \code{\link{block.plsda}},
#' \code{\link{block.spls}} and http://www.mixOmics.org/mixDIABLO for more
#' details and examples.
#' @references On multiple integration with sPLS-DA and 4 data blocks:
#'
#' Singh A., Gautier B., Shannon C., Vacher M., Rohart F., Tebbutt S. and Lê
#' Cao K.A. (2016). DIABLO: multi omics integration for biomarker discovery.
#' BioRxiv available here:
#' \url{http://biorxiv.org/content/early/2016/08/03/067611}
#'
#' On data integration:
#'
#' Tenenhaus A., Philippe C., Guillemot V, Lê Cao K.A., Grill J, Frouin V.
#' Variable selection for generalized canonical correlation analysis.
#' \emph{Biostatistics}. kxu001
#'
#' Gunther O., Shin H., Ng R. T. , McMaster W. R., McManus B. M. , Keown P. A.
#' , Tebbutt S.J. , Lê Cao K-A. , (2014) Novel multivariate methods for
#' integration of genomics and proteomics data: Applications in a kidney
#' transplant rejection study, OMICS: A journal of integrative biology, 18(11),
#' 682-95.
#'
#' mixOmics article:
#'
#' Rohart F, Gautier B, Singh A, Lê Cao K-A. mixOmics: an R package for 'omics
#' feature selection and multiple data integration. PLoS Comput Biol 13(11):
#' e1005752
#' @keywords regression multivariate
#' 
## ----------- Examples ----------- 
#' @example examples/block.splsda-example.R
## setting the document name here so internal would not force the wrong name
#' @name block.splsda
NULL
## ----------- Internal ----------- 
.block.splsda = function(X=NULL,
                         Y=NULL,
                         indY=NULL,
                         ncomp = 2,
                         keepX=NULL,
                         design=NULL,
                         scheme=c("horst", "factorial", "centroid"),
                         mode = c("regression", "canonical", "invariant", "classic"),
                         scale = TRUE,
                         init=c("svd", "svd.single") ,
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
    
    ## check entries
    mc <- .check_plsda_block(mc)
    Y.input <- mc$Y.input
    mc$Y.input <- NULL ## not needed by wrapper
    # call to '.mintWrapperBlock'
    mc[[1L]] <- quote(.mintWrapperBlock)
    result <- eval(mc)
    
    # calculate weights for each dataset
    weights = .getWeights(result$variates, indY = result$indY)
    
    # choose the desired output from 'result'
    out <- list(X = result$A[-result$indY],
                Y = Y.input,
                ind.mat = result$A[result$indY][[1]],
                ncomp = result$ncomp,
                mode = result$mode,
                keepX = result$keepX[-result$indY],
                keepY = result$keepX[result$indY][[1]],
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
                indY = result$indY,
                weights = weights,
                explained_variance = result$explained_variance)
    
    # give a class
    class(out) = c("block.splsda","block.pls","sgccda","sgcca","DA")
    
    return(invisible(out))
    
}
## ----------- Generic ----------- 
#' @export
#' @rdname block.splsda
setGeneric('block.splsda', function(data=NULL, X=NULL, Y=NULL, formula=NULL, ...) 
    standardGeneric('block.splsda'))

## ----------- Methods ----------- 
#### ANY ####
#' @export
#' @rdname block.splsda
setMethod('block.splsda', 'ANY', function(data=NULL, X=NULL, Y=NULL, formula=NULL, ...) {
    
    tryCatch(mget(names(formals()), sys.frame(sys.nframe())),
             error = function(e) stop(e$message, call. = FALSE))
    mc <- match.call()
    mc[-1L] <- lapply(mc[-1L], eval)
    
    ## check signature format
    .check_sig_ANY(mc, fun = "block.splsda")
    
    mc[[1L]] <- quote(.block.splsda)
    result <- eval(mc)
    
    .call_return(result, match.call(), fun.name = 'block.splsda')
})

#### signature(data = 'MultiAssayExperiment') ####
#' @export
#' @rdname block.splsda
setMethod('block.splsda', signature(data = 'MultiAssayExperiment'), 
          function(data=NULL, X=NULL, Y=NULL, formula=NULL, ...) {
              tryCatch(mget(names(formals()), sys.frame(sys.nframe())),
                       error = function(e) stop(e$message, call. = FALSE))
              mc <- match.call()
              mc[-1] <- lapply(mc[-1], eval)
              mc$data <- .matched_samples(mc$data)
              mc <- .get_xy(mc = mc, DA = TRUE, block = TRUE)
              mc$data <- mc$formula <- NULL 
              mc[[1L]] <- quote(.block.splsda)
              result <- eval(mc)
              .call_return(result, mc$ret.call, mcr = match.call(), fun.name = 'block.splsda')
          })

#### signature(data = 'MatchedAssayExperiment') ####
## same as MultiAssayExperiment with different and no sample matching
#' @export
#' @rdname block.splsda
setMethod('block.splsda', signature(data = 'MatchedAssayExperiment'), 
          function(data=NULL, X=NULL, Y=NULL, formula=NULL, ...) {
              tryCatch(mget(names(formals()), sys.frame(sys.nframe())),
                       error = function(e) stop(e$message, call. = FALSE))
              mc <- match.call()
              mc[-1] <- lapply(mc[-1], eval)
              mc <- .get_xy(mc = mc, DA = TRUE, block = TRUE)
              mc$data <- mc$formula <- NULL 
              mc[[1L]] <- quote(.block.splsda)
              result <- eval(mc)
              .call_return(result, mc$ret.call, mcr = match.call(), fun.name = 'block.splsda')
          })

#### signature(data != 'MultiAssayExperiment', formula = "formula") ####
#' @export
#' @rdname block.splsda
setMethod('block.splsda', signature(formula = 'formula'), 
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
              mc[[1L]] <- quote(.block.splsda)
              result <- eval(mc)
              .call_return(result, mc$ret.call, mcr = match.call(), fun.name = 'block.splsda')
          })