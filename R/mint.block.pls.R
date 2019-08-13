#############################################################################################################
# Author :
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Benoit Gautier, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
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
# mint.block.pls: perform a horizontal and vertical PLS on a combination of datasets, input as a list in X
# this function is a particular setting of .mintBlock,
# the formatting of the input is checked in .mintWrapperBlock, which then call '.mintBlock'
# ========================================================================================================

# X: a list of data sets (called 'blocks') matching on the same samples. Data in the list should be arranged in samples x variables, with samples order matching in all data sets. \code{NA}s are not allowed.
# Y: outcome
# indY: to supply if Y is missing, indicate the position of the outcome in the list X.
# study: grouping factor indicating which samples are from the same study
# ncomp: numeric vector of length the number of blocks in \code{X}. The number of components to include in the model for each block (does not necessarily need to take the same value for each block). By default set to 2 per block.
# design: the input design.
# scheme: the input scheme, one of "horst", "factorial" or ""centroid". Default to "centroid"
# mode: input mode, one of "canonical", "classic", "invariant" or "regression". Default to "regression"
# scale: boleean. If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
# init: intialisation of the algorithm, one of "svd" or "svd.single". Default to "svd"
# tol: Convergence stopping value.
# max.iter: integer, the maximum number of iterations.
# near.zero.var: boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations
# all.outputs: calculation of non-essential outputs (e.g. explained variance, loadings.Astar, etc)








#' NP-integration
#'
#' Function to integrate data sets measured on the same samples (N-integration)
#' and to combine multiple independent studies measured on the same variables
#' or predictors (P-integration) using variants of multi-group and generalised
#' PLS (unsupervised analysis).
#'
#' The function fits multi-group generalised PLS models with a specified number
#' of \code{ncomp} components. An outcome needs to be provided, either by
#' \code{Y} or by its position \code{indY} in the list of blocks \code{X}.
#'
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
#' References and more details in \code{?pls}).
#'
#' @param X A list of data sets (called 'blocks') measured on the same samples.
#' Data in the list should be arranged in samples x variables, with samples
#' order matching in all data sets.
#' @param Y Matrix or vector response for a multivariate regression framework.
#' Data should be continuous variables (see \code{mint.block.plsda} for
#' supervised classification and factor reponse)
#' @param indY To be supplied if Y is missing, indicates the position of the
#' matrix / vector response in the list \code{X}
#' @param study factor indicating the membership of each sample to each of the
#' studies being combined
#' @param ncomp the number of components to include in the model. Default to 2.
#' @param design numeric matrix of size (number of blocks) x (number of blocks)
#' with only 0 or 1 values. A value of 1 (0) indicates a relationship (no
#' relationship) between the blocks to be modelled. If \code{Y} is provided
#' instead of \code{indY}, the \code{design} matrix is changed to include
#' relationships to \code{Y}.
#' @param scheme Either "horst", "factorial" or "centroid". Default =
#' \code{horst}, see reference.
#' @param mode character string. What type of algorithm to use, (partially)
#' matching one of \code{"regression"}, \code{"canonical"}, \code{"invariant"}
#' or \code{"classic"}. See Details.
#' @param scale boleean. If scale = TRUE, each block is standardized to zero
#' means and unit variances (default: TRUE)
#' @param init Mode of initialization use in the algorithm, either by Singular
#' Value Decompostion of the product of each block of X with Y ("svd") or each
#' block independently ("svd.single"). Default = \code{svd.single}.
#' @param tol Convergence stopping value.
#' @param max.iter integer, the maximum number of iterations.
#' @param near.zero.var boolean, see the internal \code{\link{nearZeroVar}}
#' function (should be set to TRUE in particular for data with many zero
#' values). Setting this argument to FALSE (when appropriate) will speed up the
#' computations. Default value is FALSE
#' @param all.outputs boolean. Computation can be faster when some specific
#' (and non-essential) outputs are not calculated. Default = \code{TRUE}.
#' @return \code{mint.block.pls} returns an object of class \code{"mint.pls",
#' "block.pls"}, a list that contains the following components:
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
#' algorthm for each component}
#' @author Florian Rohart, Benoit Gautier, Kim-Anh Lê Cao
#' @seealso \code{\link{spls}}, \code{\link{summary}}, \code{\link{plotIndiv}},
#' \code{\link{plotVar}}, \code{\link{predict}}, \code{\link{perf}},
#' \code{\link{mint.block.spls}}, \code{\link{mint.block.plsda}},
#' \code{\link{mint.block.splsda}} and http://www.mixOmics.org/mixMINT for more
#' details.
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
#'
#' # for the purpose of this example, we create data that fit in the context of
#' # this function.
#' # We consider the training set as study1 and the test set as another
#' # independent study2.
#'
#' study = c(rep("study1",150), rep("study2",70))
#'
#' # to put the data in the MINT format, we rbind the two studies
#' mrna = rbind(breast.TCGA$data.train$mrna, breast.TCGA$data.test$mrna)
#' mirna = rbind(breast.TCGA$data.train$mirna, breast.TCGA$data.test$mirna)
#'
#' # For the purpose of this example, we create a continuous response by
#' # taking the first mrna variable, and removing it from the data
#' Y = mrna[,1]
#' mrna = mrna[,-1]
#'
#' data = list(mrna = mrna, mirna = mirna)
#'
#' # we can now apply the function
#' res = mint.block.plsda(data, Y, study=study, ncomp=2)
#'
#' res
#'
#'
#' @export mint.block.pls
mint.block.pls = function(X,
Y,
indY,
study,
ncomp = 2,
design,
scheme,
mode,
scale = TRUE,
init ,
tol = 1e-06,
max.iter = 100,
near.zero.var = FALSE,
all.outputs = TRUE)
{
    # call to '.mintWrapperBlock'
    result = .mintWrapperBlock(X=X, Y=Y, indY=indY, study=study, ncomp=ncomp, design=design, scheme=scheme, mode=mode, scale=scale,
    init=init, tol=tol, max.iter=max.iter, near.zero.var=near.zero.var, all.outputs = all.outputs)

    # choose the desired output from 'result'
    out=list(
        call = match.call(),
        X = result$A,
        Y = result$A[[1]],
        ncomp = result$ncomp,
        mode = result$mode,
        study = result$study,
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
        scale = result$scale)

    class(out) = c("mint.block.pls","block.pls","sgcca")
    return(invisible(out))

}



