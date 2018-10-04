#############################################################################################################
# Authors:
#   Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
#   Kim-Anh Le Cao, French National Institute for Agricultural Research and ARC Centre of Excellence in Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
#   Sebastien Dejean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
#
# created: 2009
# last modified:
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

tune.rcc =
function(X, 
         Y, 
         grid1 = seq(0.001, 1, length = 5), 
         grid2 = seq(0.001, 1, length = 5), 
         validation = c("loo", "Mfold"), 
         folds = 10,
         plot = TRUE)
{

    # validation des arguments #
	#--------------------------#
    if (length(dim(X)) != 2 || length(dim(Y)) != 2) 
        stop("'X' and/or 'Y' must be a numeric matrix.")
     
    X = as.matrix(X)
    Y = as.matrix(Y)
     
    if (!is.numeric(X) || !is.numeric(Y)) 
    stop("'X' and/or 'Y' must be a numeric matrix.")
		
    if (nrow(X) != nrow(Y)) 
    stop("unequal number of rows in 'X' and 'Y'.")
     
    validation = match.arg(validation)
    grid = expand.grid(grid1, grid2)
     
    if (validation == "loo")
    {
        M = nrow(X)
        folds = split(1:M, 1:M)
        cv.score = apply(grid, 1, function(lambda)
                                    {
                                        Mfold(X, Y, lambda[1], lambda[2], folds)
                                    })

    } else {
        nr = nrow(X)
        M = length(folds)
        if (is.null(folds) || !is.numeric(folds) || folds < 2 || folds > nr)
        {
            stop("Invalid number of folds.")
        } else {
            M = round(folds)
            folds = split(sample(1:nr), rep(1:M, length = nr))
        }
        cv.score = apply(grid, 1, function(lambda) 
                                    {
                                        Mfold(X, Y, lambda[1], lambda[2], folds)
                                    })
    }
     
    cv.score.grid = cbind(grid, cv.score)
    mat = matrix(cv.score, nrow = length(grid1), ncol = length(grid2))
     
    if (isTRUE(plot))
    image.tune.rcc(list(grid1 = grid1, grid2 = grid2, mat = mat))
    
    opt = cv.score.grid[cv.score.grid[, 3] == max(cv.score.grid[, 3]), ]
    cat("  lambda1 = ", opt[[1]], "\n", " lambda2 = ", opt[[2]], "\n",
	"CV-score = ", opt[[3]], "\n")
     
    out = list(opt.lambda1 = opt[[1]], opt.lambda2 = opt[[2]], 
	opt.score = opt[[3]], grid1 = grid1, grid2 = grid2, mat = mat)
    
    out$call = match.call()

    class(out) = "tune.rcc"
    return(invisible(out))
}


Mfold = function(X, Y, lambda1, lambda2, folds)
{
    
    xscore = NULL
    yscore = NULL
    M = length(folds)
    
    for (m in 1:M)
    {
        omit = folds[[m]]
        result = rcc(X[-omit, , drop = FALSE], Y[-omit, , drop = FALSE], ncomp = 1, lambda1, lambda2, method = "ridge")
        X[omit, ][is.na(X[omit, ])] = 0
        Y[omit, ][is.na(Y[omit, ])] = 0
        xscore = c(xscore, X[omit, , drop = FALSE] %*% result$loadings$X[, 1])
        yscore = c(yscore, Y[omit, , drop = FALSE] %*% result$loadings$Y[, 1])
    }
    
    cv.score = cor(xscore, yscore, use = "pairwise")
    return(invisible(cv.score))
}


