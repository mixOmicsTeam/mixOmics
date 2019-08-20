#############################################################################################################
# Authors:
#   Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Kim-Anh Le Cao, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Al J Abadi, Melbourne Integartive Genomics, The University of Melbourne, Australia
#
# created: 2009
# last modified: 2019
#
# Copyright (C) 2009
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#############################################################################################################


# ========================================================================================================
# spls: perform a sPLS
# this function is a particular setting of .mintBlock, the formatting of the input is checked in .mintWrapper
# ========================================================================================================

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
## --------------------------------------------------------------------------------------- parameters
#' @inheritParams pls
#' @param keepX Numeric vector of length \code{ncomp}, the number of variables
#' to keep in \eqn{X}-loadings. By default all variables are kept in the model.
#' @param keepY Numeric vector of length \code{ncomp}, the number of variables
#' to keep in \eqn{Y}-loadings. By default all variables are kept in the model.
## --------------------------------------------------------------------------------------- value
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
## ---------------------------------------------------------------------------------------
#' @author Sébastien Déjean, Ignacio González and Kim-Anh Lê Cao.
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
## --------------------------------------------------------------------------------------- examples
#' @example examples/spls-example.R

#' @export spls
spls = function(X,
Y,
ncomp = 2,
mode = c("regression", "canonical", "invariant", "classic"),
keepX,
keepY,
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


