#############################################################################################################
# Author :
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 22-04-2015
# last modified: 01-03-2016
#
# Copyright (C) 2015
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
# mint.pca: perform a vertical PCA on a combination of experiments, input as a matrix in X
# this function is a particular setting of .mintBlock,
# the formatting of the input is checked in .mintWrapper, which then call '.mintBlock'
# ========================================================================================================

# X: numeric matrix of predictors
# Y: numeric vector or matrix of responses
# ncomp: the number of components to include in the model. Default to 2.
# study: grouping factor indicating which samples are from the same study
# scale: boleean. If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
# tol: Convergence stopping value.
# max.iter: integer, the maximum number of iterations.
# near.zero.var: boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations








#' P-integration with Principal Component Analysis
#'
#' Function to integrate and combine multiple independent studies measured on
#' the same variables or predictors (P-integration) using a multigroup
#' Principal Component Analysis.
#'
#' \code{mint.pca} fits a vertical PCA model with \code{ncomp} components in
#' which several independent studies measured on the same variables are
#' integrated.  The \code{study} factor indicates the membership of each sample
#' in each study. We advise to only combine studies with more than 3 samples as
#' the function performs internal scaling per study.
#'
#' Missing values are handled by being disregarded during the cross product
#' computations in the algorithm without having to delete rows with missing
#' data. Alternatively, missing data can be imputed prior using the
#' \code{nipals} function.
#'
#' Useful graphical outputs are available, e.g. \code{\link{plotIndiv}},
#' \code{\link{plotLoadings}}, \code{\link{plotVar}}.
#'
#' @param X numeric matrix of predictors combining multiple independent studies
#' on the same set of predictors. \code{NA}s are allowed.
#' @param ncomp Number of components to include in the model (see Details).
#' Default to 2
#' @param study factor indicating the membership of each sample to each of the
#' studies being combined
#' @param scale boleean. If scale = TRUE, each block is standardized to zero
#' means and unit variances. Default = \code{TRUE}.
#' @param tol Convergence stopping value.
#' @param max.iter integer, the maximum number of iterations.
#' @return \code{mint.pca} returns an object of class \code{"mint.pca", "pca"},
#' a list that contains the following components:
#'
#' \item{X}{the centered and standardized original predictor matrix.}
#' \item{ncomp}{the number of components included in the model.}
#' \item{study}{The study grouping factor} \item{sdev}{the eigenvalues of the
#' covariance/correlation matrix, though the calculation is actually done with
#' the singular values of the data matrix or by using NIPALS.} \item{center,
#' scale}{the centering and scaling used, or \code{FALSE}.} \item{rotation}{the
#' matrix of variable loadings (i.e., a matrix whose columns contain the
#' eigenvectors).} \item{loadings}{same as 'rotation' to keep the mixOmics
#' spirit} \item{x}{the value of the rotated data (the centred (and scaled if
#' requested) data multiplied by the rotation/loadings matrix), also called the
#' principal components.} \item{variates}{same as 'x' to keep the mixOmics
#' spirit} \item{explained_variance}{explained variance from the multivariate
#' model, used for plotIndiv} \item{names}{list containing the names to be used
#' for individuals and variables.}
#' @author Florian Rohart, Kim-Anh Lê Cao
#' @seealso \code{\link{spls}}, \code{\link{summary}}, \code{\link{plotIndiv}},
#' \code{\link{plotVar}}, \code{\link{predict}}, \code{\link{perf}},
#' \code{\link{mint.spls}}, \code{\link{mint.plsda}}, \code{\link{mint.splsda}}
#' and http://www.mixOmics.org/mixMINT for more details.
#' @references Rohart F, Eslami A, Matigian, N, Bougeard S, Lê Cao K-A (2017).
#' MINT: A multivariate integrative approach to identify a reproducible
#' biomarker signature across multiple experiments and platforms. BMC
#' Bioinformatics 18:128.
#'
#' Eslami, A., Qannari, E. M., Kohler, A., and Bougeard, S. (2014). Algorithms
#' for multi-group PLS. J. Chemometrics, 28(3), 192-201.
#' @keywords regression multivariate
#' @examples
#'
#' res = mint.pca(X = stemcells$gene, ncomp = 3,
#' study = stemcells$study)
#'
#' plotIndiv(res, group = stemcells$celltype, legend=TRUE)
#'
#'
#' @export mint.pca
mint.pca = function(X,
ncomp = 2,
study,
scale = TRUE,
tol = 1e-06,
max.iter = 100)
{

    #-- checking general input parameters --------------------------------------#
    #---------------------------------------------------------------------------#

    #-- check that the user did not enter extra arguments
    arg.call = match.call()
    user.arg = names(arg.call)[-1]

    err = tryCatch(mget(names(formals()), sys.frame(sys.nframe())),
    error = function(e) e)

    if ("simpleError" %in% class(err))
    stop(err[[1]], ".", call. = FALSE)

    #-- X matrix
    if (is.data.frame(X))
    X = as.matrix(X)

    if (!is.matrix(X) || is.character(X))
    stop("'X' must be a numeric matrix.", call. = FALSE)

    if (any(apply(X, 1, is.infinite)))
    stop("infinite values in 'X'.", call. = FALSE)

    #-- put a names on the rows and columns of X --#
    X.names = colnames(X)
    if (is.null(X.names))
    X.names = paste("V", 1:ncol(X), sep = "")

    ind.names = rownames(X)
    if (is.null(ind.names))
    ind.names = 1:nrow(X)

    #-- ncomp
    if (is.null(ncomp))
    ncomp = min(nrow(X),ncol(X))

    ncomp = round(ncomp)

    if ( !is.numeric(ncomp) || ncomp < 1 || !is.finite(ncomp))
    stop("invalid value for 'ncomp'.", call. = FALSE)

    if (ncomp > min(ncol(X), nrow(X)))
    stop("use smaller 'ncomp'", call. = FALSE)


    #-- cheking scale
    if (!is.logical(scale))
    {
        if (!is.numeric(scale) || (length(scale) != ncol(X)))
        stop("'scale' should be either a logical value or a numeric vector of length equal to the number of columns of 'X'.",
        call. = FALSE)
    }

    #-- max.iter
    if (is.null(max.iter) || !is.numeric(max.iter) || max.iter < 1 || !is.finite(max.iter))
    stop("invalid value for 'max.iter'.", call. = FALSE)

    max.iter = round(max.iter)

    #-- tol
    if (is.null(tol) || !is.numeric(tol) || tol < 0 || !is.finite(tol))
    stop("invalid value for 'tol'.", call. = FALSE)

    #set the default study factor
    if (missing(study))
    {
        study = factor(rep(1,nrow(X)))
    } else {
        study = factor(study)
    }
    if (length(study) != nrow(X))
    stop(paste0("'study' must be a factor of length ",nrow(X),"."))

    if (any(table(study) <= 1))
    stop("At least one study has only one sample, please consider removing before calling the function again")
    if (any(table(study) < 5))
    warning("At least one study has less than 5 samples, mean centering might not do as expected")



    #-- end checking --#
    #------------------#

    # call to '.mintWrapper'
    mean_centered = mean_centering_per_study(data = X, study = study, scale = scale)
    X_mean_centered = as.matrix(mean_centered$concat.data)

    out = pca(X_mean_centered, ncomp = ncomp, max.iter = max.iter, tol = tol, scale = FALSE)

    # choose the desired output from 'result'
    out$study = study

    class(out) = c("mint.pca","pca")
    return(invisible(out))


}
