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
# mint.plsda: perform a vertical PLS-DA on a combination of experiments, input as a matrix in X
# this function is a particular setting of .mintBlock,
# the formatting of the input is checked in .mintWrapper, which then call '.mintBlock'
# ========================================================================================================

# X: numeric matrix of predictors
# Y: a factor or a class vector for the discrete outcome
# ncomp: the number of components to include in the model. Default to 2.
# mode: input mode, one of "canonical", "classic", "invariant" or "regression". Default to "regression"
# study: grouping factor indicating which samples are from the same study
# scale: boleean. If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
# tol: Convergence stopping value.
# max.iter: integer, the maximum number of iterations.
# near.zero.var: boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations
# all.outputs: calculation of non-essential outputs (e.g. explained variance, loadings.Astar, etc)








#' P-integration with Projection to Latent Structures models (PLS) with
#' Discriminant Analysis
#'
#' Function to combine multiple independent studies measured on the same
#' variables or predictors (P-integration) using variants of multi-group PLS-DA
#' for supervised classification.
#'
#' \code{mint.plsda} function fits a vertical PLS-DA models with \code{ncomp}
#' components in which several independent studies measured on the same
#' variables are integrated. The aim is to classify the discrete outcome
#' \code{Y}. The \code{study} factor indicates the membership of each sample in
#' each study. We advise to only combine studies with more than 3 samples as
#' the function performs internal scaling per study, and where all outcome
#' categories are represented.
#'
#' \code{X} can contain missing values. Missing values are handled by being
#' disregarded during the cross product computations in the algorithm
#' \code{mint.plsda} without having to delete rows with missing data.
#' Alternatively, missing data can be imputed prior using the \code{nipals}
#' function.
#'
#' The type of algorithm to use is specified with the \code{mode} argument.
#' Four PLS algorithms are available: PLS regression \code{("regression")}, PLS
#' canonical analysis \code{("canonical")}, redundancy analysis
#' \code{("invariant")} and the classical PLS algorithm \code{("classic")} (see
#' References and more details in \code{?pls}).
#'
#' Useful graphical outputs are available, e.g. \code{\link{plotIndiv}},
#' \code{\link{plotLoadings}}, \code{\link{plotVar}}.
#'
#' @param X numeric matrix of predictors combining multiple independent studies
#' on the same set of predictors. \code{NA}s are allowed.
#' @param Y A factor or a class vector indicating the discrete outcome of each
#' sample.
#' @param ncomp Number of components to include in the model (see Details).
#' Default to 2
#' @param mode character string. What type of algorithm to use, (partially)
#' matching one of \code{"regression"} or \code{"canonical"}. See Details.
#' @param study factor indicating the membership of each sample to each of the
#' studies being combined
#' @param scale boleean. If scale = TRUE, each block is standardized to zero
#' means and unit variances. Default = \code{TRUE}.
#' @param tol Convergence stopping value.
#' @param max.iter integer, the maximum number of iterations.
#' @param near.zero.var boolean, see the internal \code{\link{nearZeroVar}}
#' function (should be set to TRUE in particular for data with many zero
#' values). Default = \code{FALSE}.
#' @param all.outputs boolean. Computation can be faster when some specific
#' (and non-essential) outputs are not calculated. Default = \code{TRUE}.
#' @return \code{mint.plsda} returns an object of class \code{"mint.plsda",
#' "plsda"}, a list that contains the following components:
#'
#' \item{X}{the centered and standardized original predictor matrix.}
#' \item{Y}{original factor} \item{ind.mat}{the centered and standardized
#' original response vector or matrix.} \item{ncomp}{the number of components
#' included in the model.} \item{study}{The study grouping factor}
#' \item{mode}{the algorithm used to fit the model.} \item{variates}{list
#' containing the variates of X - global variates.} \item{loadings}{list
#' containing the estimated loadings for the variates - global loadings.}
#' \item{variates.partial}{list containing the variates of X relative to each
#' study - partial variates.} \item{loadings.partial}{list containing the
#' estimated loadings for the partial variates - partial loadings.}
#' \item{names}{list containing the names to be used for individuals and
#' variables.} \item{nzv}{list containing the zero- or near-zero predictors
#' information.} \item{iter}{Number of iterations of the algorthm for each
#' component} \item{explained_variance}{Percentage of explained variance for
#' each component and each study (note that contrary to PCA, this amount may
#' not decrease as the aim of the method is not to maximise the variance, but
#' the covariance between X and the dummy matrix Y).}
#' @author Florian Rohart, Kim-Anh Lê Cao
#' @seealso \code{\link{spls}}, \code{\link{summary}}, \code{\link{plotIndiv}},
#' \code{\link{plotVar}}, \code{\link{predict}}, \code{\link{perf}},
#' \code{\link{mint.pls}}, \code{\link{mint.spls}}, \code{\link{mint.splsda}}
#' and http://www.mixOmics.org/mixMINT for more details.
#' @references Rohart F, Eslami A, Matigian, N, Bougeard S, Lê Cao K-A (2017).
#' MINT: A multivariate integrative approach to identify a reproducible
#' biomarker signature across multiple experiments and platforms. BMC
#' Bioinformatics 18:128.
#'
#' Eslami, A., Qannari, E. M., Kohler, A., and Bougeard, S. (2014). Algorithms
#' for multi-group PLS. J. Chemometrics, 28(3), 192-201.
#'
#' mixOmics article:
#'
#' Rohart F, Gautier B, Singh A, Lê Cao K-A. mixOmics: an R package for 'omics
#' feature selection and multiple data integration. PLoS Comput Biol 13(11):
#' e1005752
#' @keywords regression multivariate
#' @examples
#'
#' res = mint.plsda(X = stemcells$gene, Y = stemcells$celltype, ncomp = 3,
#' study = stemcells$study)
#'
#' plotIndiv(res)
#'
#' #plot study-specific outputs for all studies
#' plotIndiv(res, study = "all.partial")
#'
#' \dontrun{
#' #plot study-specific outputs for study "2"
#' plotIndiv(res, study = "2", col = 1:3, legend = TRUE)
#' }
#'
#' @export mint.plsda
mint.plsda = function(X,
Y,
ncomp = 2,
mode = c("regression", "canonical", "invariant", "classic"),
study,
scale = TRUE,
tol = 1e-06,
max.iter = 100,
near.zero.var = FALSE,
all.outputs = TRUE)
{


    #-- validation des arguments --#
    # most of the checks are done in '.mintWrapper'

    if (is.null(Y))
    stop("'Y' has to be something else than NULL.")

    if (is.null(dim(Y)))
    {
        Y = factor(Y)
    }  else {
        stop("'Y' should be a factor or a class vector.")
    }
    Y.mat = unmap(Y)
    colnames(Y.mat) = levels(Y)

    X = as.matrix(X)

    if (length(study) != nrow(X))
    stop(paste0("'study' must be a factor of length ",nrow(X),"."))

    if(sum(apply(table(Y,study)!=0,2,sum)==1) >0)
    stop("At least one study only contains a single level of the multi-levels outcome Y. The MINT algorithm cannot be computed.")

    if(sum(apply(table(Y,study)==0,2,sum)>0) >0)
    warning("At least one study does not contain all the levels of the outcome Y. The MINT algorithm might not perform as expected.")

    # call to '.mintWrapper'
    result = .mintWrapper(X = X, Y = Y.mat, study = study, ncomp = ncomp, scale = scale, near.zero.var = near.zero.var, mode = mode,
    max.iter = max.iter, tol = tol, all.outputs = all.outputs)


    # choose the desired output from 'result'
    out = list(
        call = match.call(),
        X = result$A[-result$indY][[1]],
        Y = Y,
        ind.mat = result$A[result$indY][[1]],
        ncomp = result$ncomp,
        study = result$study,
        mode = result$mode,
        variates = result$variates,
        loadings = result$loadings,
        variates.partial = result$variates.partial,
        loadings.partial = result$loadings.partial,
        names = result$names,
        tol = result$tol,
        iter = result$iter,
        max.iter = result$max.iter,
        nzv = result$nzv,
        scale = result$scale,
        explained_variance = result$explained_variance
        )

    class(out) = c("mint.plsda","mint.pls","mixo_pls","DA")
    return(invisible(out))


}
