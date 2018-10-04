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
# mint.splsda: perform a vertical sPLS-DA on a combination of experiments, input as a matrix in X
# this function is a particular setting of internal_mint.block,
# the formatting of the input is checked in internal_wrapper.mint, which then call 'internal_mint.block'
# ========================================================================================================

# X: numeric matrix of predictors
# Y: a factor or a class vector for the discrete outcome
# ncomp: the number of components to include in the model. Default to 2.
# mode: input mode, one of "canonical", "classic", "invariant" or "regression". Default to "regression"
# study: grouping factor indicating which samples are from the same study
# keepX: number of \eqn{X} variables kept in the model on the last components.
# scale: boleean. If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
# tol: Convergence stopping value.
# max.iter: integer, the maximum number of iterations.
# near.zero.var: boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations
# all.outputs: calculation of non-essential outputs (e.g. explained variance, loadings.Astar, etc)


mint.splsda = function(X,
Y,
ncomp = 2,
mode = c("regression", "canonical", "invariant", "classic"),
study,
keepX = rep(ncol(X), ncomp),
scale = TRUE,
tol = 1e-06,
max.iter = 100,
near.zero.var = FALSE,
all.outputs = TRUE)
{
    
    #-- validation des arguments --#
    # most of the checks are done in 'internal_wrapper.mint'
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

    # call to 'internal_wrapper.mint'
    result = internal_wrapper.mint(X = X, Y = Y.mat, ncomp = ncomp, near.zero.var = near.zero.var, study = study, mode = mode,
    keepX = keepX, max.iter = max.iter, tol = tol, scale = scale, all.outputs = all.outputs)
    
    # choose the desired output from 'result'
    out = list(
        call = match.call(),
        X = result$A[-result$indY][[1]],
        Y = Y,
        ind.mat = result$A[result$indY][[1]],
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
        scale = result$scale,
        explained_variance = result$explained_variance)

    class(out) = c("mint.splsda","mixo_splsda","mixo_spls","DA")
    return(invisible(out))
    
}
