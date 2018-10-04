#############################################################################################################
# Authors:
#   Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
#   Kim-Anh Le Cao, French National Institute for Agricultural Research and
#   ARC Centre of Excellence ins Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
#   Leigh Coonan, Student, University of Quuensland, Australia
#   Fangzhou Yao, Student, University of Queensland, Australia
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Sebastien Dejean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
#
# created: 2009
# last modified: 08-07-2016
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


pca = function(   X,
ncomp = 2,
center = TRUE,
scale = FALSE,
max.iter = 500,
tol = 1e-09,
logratio = 'none',# one of ('none','CLR','ILR')
ilr.offset = 0.001,
V = NULL,
multilevel = NULL)
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
    
    #-- log.ratio
    choices = c('CLR', 'ILR','none')
    logratio = choices[pmatch(logratio, choices)]
    
    if (any(is.na(logratio)) || length(logratio) > 1)
    stop("'logratio' should be one of 'CLR' ,'ILR'or 'none'.", call. = FALSE)
    
    if (logratio != "none" && any(X < 0))
    stop("'X' contains negative values, you can not log-transform your data")
    
    
    #-- cheking center and scale
    if (!is.logical(center))
    {
        if (!is.numeric(center) || (length(center) != ncol(X)))
        stop("'center' should be either a logical value or a numeric vector of length equal to the number of columns of 'X'.",
        call. = FALSE)
    }
    
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
    
    #-- end checking --#
    #------------------#
    
    
    #-----------------------------#
    #-- logratio transformation --#
    
    if (is.null(V) & logratio == "ILR") # back-transformation to clr-space, will be used later to recalculate loadings etc
    V = clr.backtransfo(X)
    
    X = logratio.transfo(X = X, logratio = logratio, offset = if(logratio == "ILR") {ilr.offset} else {0})
    
    #as X may have changed
    if (ncomp > min(ncol(X), nrow(X)))
    stop("use smaller 'ncomp'", call. = FALSE)
    
    #-- logratio transformation --#
    #-----------------------------#
    
    #---------------------------------------------------------------------------#
    #-- multilevel approach ----------------------------------------------------#
    
    if (!is.null(multilevel))
    {
        # we expect a vector or a 2-columns matrix in 'Y' and the repeated measurements in 'multilevel'
        multilevel = data.frame(multilevel)
        
        if ((nrow(X) != nrow(multilevel)))
        stop("unequal number of rows in 'X' and 'multilevel'.")
        
        if (ncol(multilevel) != 1)
        stop("'multilevel' should have a single column for the repeated measurements.")
        
        multilevel[, 1] = as.numeric(factor(multilevel[, 1])) # we want numbers for the repeated measurements
        
        Xw = withinVariation(X, design = multilevel)
        X = Xw
    }
    #-- multilevel approach ----------------------------------------------------#
    #---------------------------------------------------------------------------#
    
    X = scale(X, center = center, scale = scale)
    cen = attr(X, "scaled:center")
    sc = attr(X, "scaled:scale")
    
    if (any(sc == 0))
    stop("cannot rescale a constant/zero column to unit variance.",
    call. = FALSE)
 
    is.na.X = is.na(X)
    na.X = FALSE
    if (any(is.na.X)) na.X = TRUE
    NA.X = any(is.na.X)
    
    cl = match.call()
    cl[[1]] = as.name('pca')
    result = list(call = cl, X = X, ncomp = ncomp,NA.X = NA.X,
    center = if (is.null(cen)) {FALSE} else {cen},
    scale = if (is.null(sc)) {FALSE} else {sc},
    names = list(X = X.names, sample = ind.names))
    
    
    #-- pca approach -----------------------------------------------------------#
    #---------------------------------------------------------------------------#
    
    if(logratio == 'CLR' | logratio=='none')
    {
        #-- if there are missing values use NIPALS agorithm
        if (any(is.na.X))
        {
            res = nipals(X, ncomp = ncomp, reconst = TRUE, max.iter = max.iter, tol = tol)
            result$sdev = res$eig / sqrt(max(1, nrow(X) - 1))
            names(result$sdev) = paste("PC", 1:length(result$sdev), sep = "")
            result$rotation = res$p
            dimnames(result$rotation) = list(X.names, paste("PC", 1:ncol(result$rotation), sep = ""))
            X[is.na.X] = res$rec[is.na.X]
            result$x = X %*% res$p
            dimnames(result$x) = list(ind.names, paste("PC", 1:ncol(result$x), sep = ""))
        } else {
            #-- if data is complete use singular value decomposition
            
            #-- borrowed from 'prcomp' function
            res = svd(X, nu = 0)
            
            result$sdev = res$d[1:ncomp] / sqrt(max(1, nrow(X) - 1))
            result$rotation = res$v[, 1:ncomp, drop = FALSE]
            result$x = X %*% res$v[, 1:ncomp, drop = FALSE]
        }
    } else {
        # if 'ILR', transform data and then back transform in clr space (from RobCompositions package)
        # data have been transformed above
        res = svd(X, nu = max(1, nrow(X) - 1))
        if (ncomp < ncol(X))
        {
            result$sdev = res$d[1:ncomp] / sqrt(max(1, nrow(X) - 1))  # Note: what differs with RobCompo is that they use: cumsum(eigen(cov(X))$values)/sum(eigen(cov(X))$values)
            # calculate loadings using back transformation to clr-space
            result$rotation = V %*% res$v[, 1:ncomp, drop = FALSE]
            # extract component score from the svd, multiply matrix by vector using diag, NB: this differ from our mixOmics PCA calculations
            # NB: this differ also from Filmoser paper, but ok from their code: scores are unchanged
            result$x = res$u[, 1:ncomp, drop = FALSE] %*% diag(res$d[1:ncomp, drop = FALSE])
        } else {
            result$sdev = res$d / sqrt(max(1, nrow(X) - 1))
            result$rotation = V %*% res$v
            result$x = res$u%*% diag(res$d)
        }
    }
    
    names(result$sdev) = paste("PC", 1:length(result$sdev), sep = "")
    dimnames(result$rotation) = list(X.names, paste("PC", 1:ncol(result$rotation), sep = ""))
    dimnames(result$x) = list(ind.names, paste("PC", 1:ncol(result$x), sep = ""))
    
    result$var.tot=sum(X^2 / max(1, nrow(X) - 1))# same as all res$d, or variance after nipals replacement of the missing values
    
    # to be similar to other methods, add loadings and variates as outputs
    result$loadings = list(X=result$rotation)
    result$variates = list(X=result$x)
    
    # output multilevel if needed
    if(!is.null(multilevel))
    result=c(result, list(Xw = Xw, design = multilevel))
    
    class(result) = c("pca","prcomp")
    if(!is.null(multilevel))
    class(result)=c("mlpca",class(result))
    
    #calcul explained variance
    result$explained_variance = result$sdev^2 / result$var.tot
    result$cum.var = cumsum(result$explained_variance)
    
    return(invisible(result))
}




