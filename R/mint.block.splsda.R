# ========================================================================================================
# mint.block.splsda: perform a horizontal and vertical sPLS-DA on a combination of datasets, input as a list in X
# this function is a particular setting of internal_mint.block,
# the formatting of the input is checked in internal_wrapper.mint.block, which then call 'internal_mint.block'
# ========================================================================================================

#' NP-integration with Discriminant Analysis and variable selection
#' 
#' Function to integrate data sets measured on the same samples (N-integration)
#' and to combine multiple independent studies measured on the same variables
#' or predictors (P-integration) using variants of sparse multi-group and
#' generalised PLS-DA for supervised classification and variable selection.
#' 
#' The function fits sparse multi-group generalised PLS Discriminant Analysis
#' models with a specified number of \code{ncomp} components. A factor
#' indicating the discrete outcome needs to be provided, either by \code{Y} or
#' by its position \code{indY} in the list of blocks \code{X}.
#' 
#' \code{X} can contain missing values. Missing values are handled by being
#' disregarded during the cross product computations in the algorithm
#' \code{block.pls} without having to delete rows with missing data.
#' Alternatively, missing data can be imputed prior using the
#' \code{\link{impute.nipals}} function.
#' 
#' The type of algorithm to use is specified with the \code{mode} argument.
#' Four PLS algorithms are available: PLS regression \code{("regression")}, PLS
#' canonical analysis \code{("canonical")}, redundancy analysis
#' \code{("invariant")} and the classical PLS algorithm \code{("classic")} (see
#' References and more details in \code{?pls}).
#' 
#' @inheritParams mint.block.plsda
#' @inheritParams mint.block.spls
#' @return \code{mint.block.splsda} returns an object of class
#' \code{"mint.splsda", "block.splsda"}, a list that contains the following
#' components:
#' 
#' \item{X}{the centered and standardized original predictor matrix.}
#' \item{Y}{the centered and standardized original response vector or matrix.}
#' \item{ncomp}{the number of components included in the model for each block.}
#' \item{mode}{the algorithm used to fit the model.} \item{mat.c}{matrix of
#' coefficients from the regression of X / residual matrices X on the
#' X-variates, to be used internally by \code{predict}.} \item{variates}{list
#' containing the \eqn{X} and \eqn{Y} variates.} \item{loadings}{list
#' containing the estimated loadings for the variates.} \item{names}{list
#' containing the names to be used for individuals and variables.}
#' \item{nzv}{list containing the zero- or near-zero predictors information.}
#' \item{tol}{the tolerance used in the iterative algorithm, used for
#' subsequent S3 methods} \item{max.iter}{the maximum number of iterations,
#' used for subsequent S3 methods} \item{iter}{Number of iterations of the
#' algorithm for each component}
#' Note that the argument 'scheme' has now been hardcoded to 'horst' and 'init' to 'svd.single'. 
#' @author Florian Rohart, Benoit Gautier, Kim-Anh Lê Cao, Al J Abadi
#' @seealso \code{\link{spls}}, \code{\link{summary}}, \code{\link{plotIndiv}},
#' \code{\link{plotVar}}, \code{\link{predict}}, \code{\link{perf}},
#' \code{\link{mint.block.spls}}, \code{\link{mint.block.plsda}},
#' \code{\link{mint.block.pls}} and http://www.mixOmics.org/mixMINT for more
#' details.
#' @references On multi-group PLS: Rohart F, Eslami A, Matigian, N, Bougeard S,
#' Lê Cao K-A (2017). MINT: A multivariate integrative approach to identify a
#' reproducible biomarker signature across multiple experiments and platforms.
#' BMC Bioinformatics 18:128.
#' 
#' Eslami, A., Qannari, E. M., Kohler, A., and Bougeard, S. (2014). Algorithms
#' for multi-group PLS. J. Chemometrics, 28(3), 192-201.
#' 
#' On multiple integration with sparse PLSDA: Singh A., Gautier B., Shannon C.,
#' Vacher M., Rohart F., Tebbutt S. and Lê Cao K.A. (2016). DIABLO: multi omics
#' integration for biomarker discovery. BioRxiv available here:
#' \url{http://biorxiv.org/content/early/2016/08/03/067611}
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
#' @export
#' @examples
#' 
#' data(breast.TCGA)
#' 
#' # for the purpose of this example, we consider the training set as study1 and
#' # the test set as another independent study2.
#' study = c(rep("study1",150), rep("study2",70))
#' 
#' mrna = rbind(breast.TCGA$data.train$mrna, breast.TCGA$data.test$mrna)
#' mirna = rbind(breast.TCGA$data.train$mirna, breast.TCGA$data.test$mirna)
#' data = list(mrna = mrna, mirna = mirna)
#' 
#' Y = c(breast.TCGA$data.train$subtype, breast.TCGA$data.test$subtype)
#' 
#' res = mint.block.splsda(data,Y,study=study,
#' keepX = list(mrna=c(10,10), mirna=c(20,20)),ncomp=2)
#' 
#' res
#' 
#' 
mint.block.splsda <- function(X,
                              Y,
                              indY,
                              study,
                              ncomp = 2,
                              keepX,
                              design,
                              scale = TRUE,
                              tol = 1e-06,
                              max.iter = 100,
                              near.zero.var = FALSE,
                              all.outputs = TRUE
)
{
    if (!missing(Y))
    {
        if (is.null(dim(Y))) {
            Y = as.factor(Y)
        }  else {
            stop("'Y' should be a factor or a class vector.")
        }
        if (nlevels(Y) == 1)
            stop("'Y' should be a factor with more than one level")
        
        
        Y.input = Y
        Y = unmap(Y)
        colnames(Y) = paste0("Y", 1:ncol(Y))
        rownames(Y) = rownames(X[[1]])
        
    } else if (!missing(indY)) {
        temp = X[[indY]] #not called Y to not be an input of the wrapper.sparse.mint.block
        if (is.null(dim(temp))) {
            temp = as.factor(temp)
        }  else {
            stop("'Y' should be a factor or a class vector.")
        }
        if (nlevels(temp) == 1)
            stop("'X[[indY]]' should be a factor with more than one level")
        
        Y.input = temp
        X[[indY]] = unmap(temp)
        rownames(X[[indY]]) = rownames(X[[ifelse(indY == 1, 2, 1)]])
    } else if (missing(indY)) {
        stop("Either 'Y' or 'indY' is needed")
    }
    
    # call to 'internal_wrapper.mint.block'
    result <- internal_wrapper.mint.block(
        X = X,
        Y = Y,
        indY = indY,
        study = study,
        ncomp = ncomp,
        keepX = keepX,
        design = design,
        scheme = 'horst',
        mode = 'regression',
        scale = scale,
        init = 'svd.single',
        tol = tol,
        max.iter = max.iter,
        near.zero.var = near.zero.var,
        all.outputs = all.outputs,
        DA = TRUE
    )
    
    # calculate weights for each dataset
    weights = get.weights(result$variates, indY = result$indY)
    
    # choose the desired output from 'result'
    out = list(
        call = match.call(),
        X = result$A[-result$indY],
        Y = Y.input,
        indY = result$indY,
        weights = weights,
        ind.mat = result$A[result$indY][[1]],
        ncomp = result$ncomp,
        mode = result$mode,
        study = result$study,
        keepX = result$keepA[-result$indY],
        keepY = result$keepA[result$indY][[1]],
        variates = result$variates,
        loadings = result$loadings,
        variates.partial = result$variates.partial,
        loadings.partial = result$loadings.partial,
        names = result$names,
        init = result$init,
        tol = result$tol,
        iter = result$iter,
        max.iter = result$max.iter,
        nzv = result$nzv,
        scale = result$scale
    )
    
    class(out) = c("mint.block.splsda","mint.block.spls","block.spls","sgccda","sgcca","DA")
    return(invisible(out))
    
}
