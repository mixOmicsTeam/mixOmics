#############################################################################################################
# Authors:
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Kim-Anh Le Cao, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD

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
# plsda: perform a PLS-DA
# this function is a particular setting of internal_mint.block, the formatting of the input is checked in internal_wrapper.mint
# ========================================================================================================

# X: numeric matrix of predictors
# Y: a factor or a class vector for the discrete outcome
# ncomp: the number of components to include in the model. Default to 2.
# scale: boleean. If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
# tol: Convergence stopping value.
# max.iter: integer, the maximum number of iterations.
# near.zero.var: boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations
# logratio: one of "none", "CLR"
# multilevel: repeated measurement. `multilevel' is passed to multilevel(design = ) in withinVariation. Y is ommited and shouldbe included in `multilevel'
# all.outputs: calculation of non-essential outputs (e.g. explained variance, loadings.Astar, etc)


plsda = function(X,
Y,
ncomp = 2,
scale = TRUE,
mode = c("regression", "canonical", "invariant", "classic"),
tol = 1e-06,
max.iter = 100,
near.zero.var = FALSE,
logratio = "none",   # one of "none", "CLR"
multilevel = NULL,
all.outputs = TRUE)
{
    
    #-- validation des arguments --#
    # most of the checks are done in the wrapper.mint.spls.hybrid function
    if (is.null(multilevel))
    {
        if (is.null(Y))
        stop("'Y' has to be something else than NULL.")
        
        if (is.null(dim(Y)))
        {
            Y = factor(Y)
        } else {
            stop("'Y' should be a factor or a class vector.")
        }
        
        if (nlevels(Y) == 1)
        stop("'Y' should be a factor with more than one level")
        
        Y.mat = unmap(Y)
        colnames(Y.mat) = levels(Y)
        
    } else {
        # we expect a vector or a 2-columns matrix in 'Y' and the repeated measurements in 'multilevel'
        multilevel = data.frame(multilevel)
        
        if ((nrow(X) != nrow(multilevel)))
        stop("unequal number of rows in 'X' and 'multilevel'.")
        
        if (ncol(multilevel) != 1)
        stop("'multilevel' should have a single column for the repeated measurements, other factors should be included in 'Y'.")
        
        if (!is.null(ncol(Y)) && !ncol(Y) %in% c(0,1,2))# multilevel 1 or 2 factors
        stop("'Y' should either be a factor, a single column data.frame containing a factor, or a 2-columns data.frame containing 2 factors.")
        multilevel = data.frame(multilevel, Y)
        multilevel[, 1] = as.numeric(factor(multilevel[, 1])) # we want numbers for the repeated measurements
        
        Y.mat = NULL
    }

    
    # call to 'internal_wrapper.mint'
    result = internal_wrapper.mint(X = X, Y = Y.mat, ncomp = ncomp, scale = scale, near.zero.var = near.zero.var, mode = mode,
    max.iter = max.iter, tol = tol, logratio = logratio, multilevel = multilevel, DA = TRUE, all.outputs=all.outputs)

    # choose the desired output from 'result'
    out = list(
        call = match.call(),
        X = result$A[-result$indY][[1]],
        Y = if (is.null(multilevel))
            {
                Y
            } else {
                result$Y.factor
            },
        ind.mat = result$A[result$indY][[1]],
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
        explained_variance = result$explained_variance,#[-result$indY],
        input.X = result$input.X,
        mat.c = result$mat.c#,
        )
    
    class(out) = c("mixo_plsda","mixo_pls","DA")
    # output if multilevel analysis
    if (!is.null(multilevel))
    {
        out$multilevel = multilevel
        class(out) = c("mixo_mlplsda",class(out))
    }

    return(invisible(out))
}

