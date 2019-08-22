#############################################################################################################
# Author :
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 22-04-2015
# last modified: 05-10-2017
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
# mint.spls: perform a vertical sPLS on a combination of experiments, input as a matrix in X
# this function is a particular setting of .mintBlock,
# the formatting of the input is checked in .mintWrapper, which then call '.mintBlock'
# ========================================================================================================

# X: numeric matrix of predictors
# Y: numeric vector or matrix of responses
# ncomp: the number of components to include in the model. Default to 2.
# mode: input mode, one of "canonical", "classic", "invariant" or "regression". Default to "regression"
# study: grouping factor indicating which samples are from the same study
# keepX: number of \eqn{X} variables kept in the model on the last components.
# keepY: number of \eqn{Y} variables kept in the model on the last components.
# scale: boleean. If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
# tol: Convergence stopping value.
# max.iter: integer, the maximum number of iterations.
# near.zero.var: boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations
# all.outputs: calculation of non-essential outputs (e.g. explained variance, loadings.Astar, etc)








#' P-integration with variable selection
#'
#' Function to integrate and combine multiple independent studies measured on
#' the same variables or predictors (P-integration) using variants of
#' multi-group sparse PLS for variable selection (unsupervised analysis).
#'
#' \code{mint.spls} fits a vertical sparse PLS-DA models with \code{ncomp}
#' components in which several independent studies measured on the same
#' variables are integrated. The aim is to explain the continuous outcome
#' \code{Y} and selecting correlated features between both data sets \code{X}
#' and \code{Y}. The \code{study} factor indicates the membership of each
#' sample in each study. We advise to only combine studies with more than 3
#' samples as the function performs internal scaling per study.
#'
#' Multi (continuous)response are supported. \code{X} and \code{Y} can contain
#' missing values. Missing values are handled by being disregarded during the
#' cross product computations in the algorithm \code{mint.spls} without having
#' to delete rows with missing data. Alternatively, missing data can be imputed
#' prior using the \code{nipals} function.
#'
#' The type of algorithm to use is specified with the \code{mode} argument.
#' Four PLS algorithms are available: PLS regression \code{("regression")}, PLS
#' canonical analysis \code{("canonical")}, redundancy analysis
#' \code{("invariant")} and the classical PLS algorithm \code{("classic")} (see
#' References and more details in \code{?pls}).
#'
#' Variable selection is performed on each component for each block of
#' \code{X}, and for \code{Y} if specified, via input parameter \code{keepX}
#' and \code{keepY}.
#'
#' Useful graphical outputs are available, e.g. \code{\link{plotIndiv}},
#' \code{\link{plotLoadings}}, \code{\link{plotVar}}.
#'
#' @param X numeric matrix of predictors combining multiple independent studies
#' on the same set of predictors. \code{NA}s are allowed.
#' @param Y Matrix or vector response for a multivariate regression framework.
#' Data should be continuous variables (see \code{mint.splsda} for supervised
#' classification and factor reponse)
#' @param ncomp Number of components to include in the model. Default to 2
#' @param mode character string. What type of algorithm to use, (partially)
#' matching one of \code{"regression"}, \code{"canonical"}, \code{"invariant"}
#' or \code{"classic"}. See Details.
#' @param study grouping factor indicating which samples are from the same
#' study
#' @param keepX numeric vector indicating the number of variables to select in
#' \code{X} on each component. By default all variables are kept in the model.
#' @param keepY numeric vector indicating the number of variables to select in
#' \code{Y} on each component. By default all variables are kept in the model.
#' @param scale boleean. If scale = TRUE, each block is standardized to zero
#' means and unit variances. Default = \code{TRUE}.
#' @param tol Convergence stopping value.
#' @param max.iter integer, the maximum number of iterations.
#' @param near.zero.var boolean, see the internal \code{\link{nearZeroVar}}
#' function (should be set to TRUE in particular for data with many zero
#' values). Default = \code{FALSE}.
#' @param all.outputs boolean. Computation can be faster when some specific
#' (and non-essential) outputs are not calculated. Default = \code{TRUE}.
#' @return \code{mint.spls} returns an object of class
#' \code{"mint.spls","spls"}, a list that contains the following components:
#'
#' \item{X}{numeric matrix of predictors combining multiple independent studies
#' on the same set of predictors. \code{NA}s are allowed.} \item{Y}{the
#' centered and standardized original response vector or matrix.}
#' \item{ncomp}{the number of components included in the model.}
#' \item{study}{The study grouping factor} \item{mode}{the algorithm used to
#' fit the model.} \item{keepX}{Number of variables used to build each
#' component of X} \item{keepY}{Number of variables used to build each
#' component of Y} \item{variates}{list containing the variates of X - global
#' variates.} \item{loadings}{list containing the estimated loadings for the
#' variates - global loadings.} \item{variates.partial}{list containing the
#' variates of X relative to each study - partial variates.}
#' \item{loadings.partial}{list containing the estimated loadings for the
#' partial variates - partial loadings.} \item{names}{list containing the names
#' to be used for individuals and variables.} \item{nzv}{list containing the
#' zero- or near-zero predictors information.} \item{iter}{Number of iterations
#' of the algorthm for each component} \item{explained_variance}{Percentage of
#' explained variance for each component and each study (note that contrary to
#' PCA, this amount may not decrease as the aim of the method is not to
#' maximise the variance, but the covariance between data sets).}
#' @author Florian Rohart, Kim-Anh Lê Cao
#' @seealso \code{\link{spls}}, \code{\link{summary}}, \code{\link{plotIndiv}},
#' \code{\link{plotVar}}, \code{\link{predict}}, \code{\link{perf}},
#' \code{\link{mint.pls}}, \code{\link{mint.plsda}}, \code{\link{mint.splsda}}
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
#' \dontrun{
#' library(mixOmicsData)
#'
#' # for the purpose of this example, we artificially
#' # create a continuous response Y by taking gene 1.
#'
#' res = mint.spls(X = stemcells$gene[,-1], Y = stemcells$gene[,1], ncomp = 3,
#' keepX = c(10, 5, 15), study = stemcells$study)
#'
#' plotIndiv(res)
#'
#' #plot study-specific outputs for all studies
#' plotIndiv(res, study = "all.partial")
#'

#' #plot study-specific outputs for study "2"
#' plotIndiv(res, study = "2", col = 1:3, legend = TRUE)
#' }
#'
#' @export mint.spls
mint.spls = function(X,
Y,
ncomp = 2,
mode = c("regression", "canonical", "invariant", "classic"),
study,
keepX=rep(ncol(X), ncomp),
keepY=rep(ncol(Y), ncomp),
scale = TRUE,
tol = 1e-06,
max.iter = 100,
near.zero.var = FALSE,
all.outputs = TRUE)
{

    # call to '.mintWrapper'
    result = .mintWrapper(X = X, Y = Y, ncomp = ncomp, scale = scale, near.zero.var = near.zero.var, study = study, mode = mode,
    keepX = keepX, keepY = keepY,
    max.iter = max.iter, tol = tol, all.outputs = all.outputs)

    # choose the desired output from 'result'
    out = list(
        call = match.call(),
        X = result$A[-result$indY][[1]],
        Y = result$A[result$indY][[1]],
        ncomp = result$ncomp,
        study = result$study,
        mode = result$mode,
        keepX = result$keepX,
        keepY = result$keepY,
        variates = result$variates,
        loadings = result$loadings,
        variates.partial = result$variates.partial,
        loadings.partial = result$loadings.partial,
        names  =  result$names,
        tol = result$tol,
        iter = result$iter,
        max.iter = result$max.iter,
        nzv = result$nzv,
        scale = scale,
        explained_variance = result$explained_variance)

    class(out) = c("mint.spls","mixo_spls")
    return(invisible(out))

}
