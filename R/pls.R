# ========================================================================================================
# pls: perform a PLS
# this function is a particular setting of internal_mint.block, the formatting of the input is checked in internal_wrapper.mint
# ========================================================================================================

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
#' @param X Numeric matrix of predictors. \code{NA}s are allowed.
#' @param Y Numeric vector or matrix of responses (for multi-response models).
#' \code{NA}s are allowed.
#' @param ncomp Integer, the number of components to include in the model. Default to 2.
#' @param scale Logical. If scale = TRUE, each block is standardized to zero
#' means and unit variances (default: TRUE)
#' @param mode Character string. What type of algorithm to use, (partially)
#' matching one of \code{"regression"}, \code{"canonical"}, \code{"invariant"}
#' or \code{"classic"}. See Details.
#' @param tol Numeric, convergence stopping value.
#' @param max.iter Integer, the maximum number of iterations.
#' @param near.zero.var Logical, see the internal \code{\link{nearZeroVar}}
#' function (should be set to TRUE in particular for data with many zero
#' values). Setting this argument to FALSE (when appropriate) will speed up the
#' computations. Default value is FALSE.
#' @param logratio Character, one of ('none','CLR') specifies the log ratio transformation
#' to deal with compositional values that may arise from specific normalisation
#' in sequencing data. Default to 'none'.
#' @param multilevel Numeric, design matrix for repeated measurement analysis, where
#' multilevel decomposition is required. For a one factor decomposition, the
#' repeated measures on each individual, i.e. the individuals ID is input as
#' the first column. For a 2 level factor decomposition then 2nd AND 3rd
#' columns indicate those factors. See examples in \code{?spls}).
#' @param all.outputs Logical. Computation can be faster when some specific
#' (and non-essential) outputs are not calculated. Default = \code{TRUE}.
#' @return \code{pls} returns an object of class \code{"pls"}, a list that
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
#' \item{iter}{Number of iterations of the algorithm for each component}
#' \item{max.iter}{the maximum number of iterations, used for subsequent S3
#' methods} \item{nzv}{list containing the zero- or near-zero predictors
#' information.} \item{scale}{whether scaling was applied per predictor.}
#' \item{logratio}{whether log ratio transformation for relative proportion
#' data was applied, and if so, which type of transformation.}
#' \item{prop_expl_var}{The proportion of the variance explained by each
#' variate / component divided by the total variance in the \code{data} (after
#' removing the possible missing values) using the definition of 'redundancy'.
#' Note that contrary to \code{PCA}, this amount may not decrease in the
#' following components as the aim of the method is not to maximise the
#' variance, but the covariance between data sets (including the dummy matrix
#' representation of the outcome variable in case of the supervised
#' approaches).}
#' \item{input.X}{numeric matrix of predictors in X that was input, before any
#' scaling / logratio / multilevel transformation.} \item{mat.c}{matrix of
#' coefficients from the regression of X / residual matrices X on the
#' X-variates, to be used internally by \code{predict}.}
#' \item{defl.matrix}{residual matrices X for each dimension.}
#' @author Sébastien Déjean, Ignacio González, Florian Rohart, Kim-Anh Lê Cao,
#' Al J Abadi
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
#' @export
#' @examples
#' data(linnerud)
#' X <- linnerud$exercise
#' Y <- linnerud$physiological
#' linn.pls <- pls(X, Y, mode = "classic")
#' 
#' \dontrun{
#' data(liver.toxicity)
#' X <- liver.toxicity$gene
#' Y <- liver.toxicity$clinic
#' toxicity.pls <- pls(X, Y, ncomp = 3)
#' }
#' 
pls <- function(X,
                Y,
                ncomp = 2,
                scale = TRUE,
                mode = c("regression", "canonical", "invariant", "classic"),
                tol = 1e-06,
                max.iter = 100,
                near.zero.var = FALSE,
                logratio = "none",
                # one of "none", "CLR"
                multilevel = NULL,
                all.outputs = TRUE)
{
    # call to 'internal_wrapper.mint'
    result = internal_wrapper.mint(
        X = X,
        Y = Y,
        ncomp = ncomp,
        scale = scale,
        near.zero.var = near.zero.var,
        mode = mode,
        max.iter = max.iter,
        tol = tol,
        logratio = logratio,
        multilevel = multilevel,
        DA = FALSE,
        all.outputs = all.outputs
    )
    
    # choose the desired output from 'result'
    out = list(
        call = match.call(),
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
        prop_expl_var = result$prop_expl_var,
        input.X = result$input.X,
        mat.c = result$mat.c#,
        #defl.matrix = result$defl.matrix
    )
    
    class(out) = c("mixo_pls")
    # output if multilevel analysis
    if (!is.null(multilevel))
    {
        out$multilevel = multilevel
        class(out) = c("mixo_mlpls",class(out))
    }
    
    return(invisible(out))
    
}
