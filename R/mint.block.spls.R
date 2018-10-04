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
# mint.block.spls: perform a horizontal and vertical sPLS on a combination of datasets, input as a list in X
# this function is a particular setting of internal_mint.block,
# the formatting of the input is checked in internal_wrapper.mint.block, which then call 'internal_mint.block'
# ========================================================================================================

# X: a list of data sets (called 'blocks') matching on the same samples. Data in the list should be arranged in samples x variables, with samples order matching in all data sets. \code{NA}s are not allowed.
# Y: outcome
# indY: to supply if Y is missing, indicate the position of the outcome in the list X.
# study: grouping factor indicating which samples are from the same study
# ncomp: numeric vector of length the number of blocks in \code{X}. The number of components to include in the model for each block (does not necessarily need to take the same value for each block). By default set to 2 per block.
# keepX: A vector of same length as X.  Each entry keepX[i] is the number of X[[i]]-variables kept in the model.
# keepY: Only used if Y is provided. Each entry keepY[i] is the number of Y-variables kept in the model on the last components.
# design: the input design.
# scheme: the input scheme, one of "horst", "factorial" or ""centroid". Default to "centroid"
# mode: input mode, one of "canonical", "classic", "invariant" or "regression". Default to "regression"
# scale: boleean. If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
# init: intialisation of the algorithm, one of "svd" or "svd.single". Default to "svd"
# tol: Convergence stopping value.
# max.iter: integer, the maximum number of iterations.
# near.zero.var: boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations
# all.outputs: calculation of non-essential outputs (e.g. explained variance, loadings.Astar, etc)



mint.block.spls = function(X,
Y,
indY,
study,
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
    # call to 'internal_wrapper.mint.block'
    result = internal_wrapper.mint.block(X=X, Y=Y, indY=indY, study=study, ncomp=ncomp,
    keepX=keepX, keepY=keepY,
    design=design, scheme=scheme, mode=mode, scale=scale, init=init, tol=tol,
    max.iter=max.iter, near.zero.var=near.zero.var, all.outputs = all.outputs)
    
    # choose the desired output from 'result'
    out = list(
        call = match.call(),
        X = result$A,
        Y = result$A[[1]],
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
        scale = result$scale)

    class(out) = c("mint.block.spls","block.spls","sgcca")
    return(invisible(out))
    
}



