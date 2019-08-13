################################################################################
# Authors:
#   Florian Rohart,
#   Benoit Gautier,
#   Kim-Anh Le Cao,
#
# created: 22-04-2015
# last modified: 04-10-2017
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
################################################################################


# =============================================================================
# block.spls: perform a horizontal PLS-DA on a combination of datasets,
#   input as a list in X
#   this function is a particular setting of .mintBlock,
#   the formatting of the input is checked in .mintWrapperBlock
# =============================================================================

# X: a list of data sets (called 'blocks') matching on the same samples.
#   Data in the list should be arranged in samples x variables,
#   with samples order matching in all data sets. \code{NA}s are not allowed.
# Y: outcome.
# indY: to supply if Y is missing, indicate the position of the outcome in X.
# ncomp: numeric vector of length the number of blocks in \code{X}.
#   The number of components to include in the model for each block
#   (does not necessarily need to take the same value for each block).
#   By default set to 2 per block.
# keepX: A vector of same length as X.
#   Each entry keepX[i] is the number of X[[i]]-variables kept in the model.
# keepY: Only used if Y is provided.
#   Each entry keepY[i] is the number of Y-variables kept in the model.
# design: the input design.
# scheme: the input scheme, one of "horst", "factorial" or ""centroid".
#   Default to "centroid"
# mode: input mode, one of "canonical", "classic", "invariant" or "regression".
#   Default to "regression"
# scale: boleean. If scale = TRUE, each block is standardized to zero means
#   and unit variances (default: TRUE).
# init: intialisation of the algorithm, one of "svd" or "svd.single".
#   Default to "svd"
# tol: Convergence stopping value.
# max.iter: integer, the maximum number of iterations.
# near.zero.var: boolean, see the internal \code{\link{nearZeroVar}} function
#   (should be set to TRUE in particular for data with many zero values).
# all.outputs: calculation of non-essential outputs
#   (e.g. explained variance, loadings.Astar, etc)








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
#' @param X A list of data sets (called 'blocks') measured on the same samples.
#' Data in the list should be arranged in matrices, samples x variables, with
#' samples order matching in all data sets.
#' @param Y Matrix response for a multivariate regression framework. Data
#' should be continuous variables (see block.splsda for supervised
#' classification and factor reponse)
#' @param indY To supply if Y is missing, indicates the position of the matrix
#' response in the list \code{X}
#' @param ncomp the number of components to include in the model. Default to 2.
#' Applies to all blocks.
#' @param keepX A list of same length as X.  Each entry is the number of
#' variables to select in each of the blocks of X for each component. By
#' default all variables are kept in the model.
#' @param keepY Only if Y is provided. Each entry is the number of variables to
#' select in each of the blocks of Y for each component.
#' @param design numeric matrix of size (number of blocks in X) x (number of
#' blocks in X) with values between 0 and 1. Each value indicates the strenght
#' of the relationship to be modelled between two blocks; a value of 0
#' indicates no relationship, 1 is the maximum value. If \code{Y} is provided
#' instead of \code{indY}, the \code{design} matrix is changed to include
#' relationships to \code{Y}.
#' @param scheme Either "horst", "factorial" or "centroid". Default =
#' \code{horst}, see reference.
#' @param mode character string. What type of algorithm to use, (partially)
#' matching one of \code{"regression"}, \code{"canonical"}, \code{"invariant"}
#' or \code{"classic"}. See Details. Default = \code{regression}.
#' @param scale boleean. If scale = TRUE, each block is standardized to zero
#' means and unit variances. Default = \code{TRUE}.
#' @param init Mode of initialization use in the algorithm, either by Singular
#' Value Decompostion of the product of each block of X with Y ("svd") or each
#' block independently ("svd.single"). Default = \code{svd.single}.
#' @param tol Convergence stopping value.
#' @param max.iter integer, the maximum number of iterations.
#' @param near.zero.var boolean, see the internal \code{\link{nearZeroVar}}
#' function (should be set to TRUE in particular for data with many zero
#' values). Default = \code{FALSE}.
#' @param all.outputs boolean. Computation can be faster when some specific
#' (and non-essential) outputs are not calculated. Default = \code{TRUE}.
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
#' @author Florian Rohart, Benoit Gautier, Kim-Anh Lê Cao
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
#' @examples
#'
#'
#' # Example with multi omics TCGA study
#' # -----------------------------
#' # this is the X data as a list of mRNA and miRNA; the Y data set is a single data set of proteins
#' data = list(mrna = breast.TCGA$data.train$mrna, mirna = breast.TCGA$data.train$mirna)
#' # set up a full design where every block is connected
#' design = matrix(1, ncol = length(data), nrow = length(data),
#' dimnames = list(names(data), names(data)))
#' diag(design) =  0
#' design
#' # set number of component per data set
#' ncomp = c(2)
#' # set number of variables to select, per component and per data set (this is set arbitrarily)
#' list.keepX = list(mrna = rep(20, 2), mirna = rep(10,2))
#' list.keepY = c(rep(10, 2))
#'
#' TCGA.block.spls = block.spls(X = data, Y = breast.TCGA$data.train$protein,
#' ncomp = ncomp, keepX = list.keepX, keepY = list.keepY, design = design)
#' TCGA.block.spls
#' # in plotindiv we color the samples per breast subtype group but the method is unsupervised!
#' plotIndiv(TCGA.block.spls, group =  breast.TCGA$data.train$subtype, ind.names = FALSE)
#' # illustrates coefficient weights in each block
#' plotLoadings(TCGA.block.spls, ncomp = 1)
#' plotVar(TCGA.block.spls, style = 'graphics', legend = TRUE)
#' network(TCGA.block.spls)
#'
#'
#' @export block.spls
block.spls = function(X,
Y,
indY,
ncomp = 2,
keepX,
keepY,
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
    result = .mintWrapperBlock(X=X, Y=Y, indY=indY, ncomp=ncomp,
    keepX=keepX, keepY=keepY,
    design=design, scheme=scheme, mode=mode, scale=scale,
    init=init, tol=tol, max.iter=max.iter, near.zero.var=near.zero.var,
    all.outputs = all.outputs)

    # calculate weights for each dataset
    weights = .getWeights(result$variates, indY = result$indY)

    # choose the desired output from 'result'
    out=list(call = match.call(),
        X = result$A,
        indY = result$indY,
        ncomp = result$ncomp,
        mode = result$mode,
        keepX = result$keepA[-result$indY],
        keepY = result$keepA[result$indY][[1]],
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



