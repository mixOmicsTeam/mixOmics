#############################################################################################################
# Authors:
#   Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
#   Kim-Anh Le Cao, French National Institute for Agricultural Research and
#   ARC Centre of Excellence ins Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
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

nipals = function (X, ncomp = 1, reconst = FALSE, max.iter = 500, tol = 1e-09)
{
    
    #-- X matrix
    if (is.data.frame(X))
    X = as.matrix(X)
    
    if (!is.matrix(X) || is.character(X))
    stop("'X' must be a numeric matrix.", call. = FALSE)
    
    if (any(apply(X, 1, is.infinite)))
    stop("infinite values in 'X'.", call. = FALSE)
    
    nc = ncol(X)
    nr = nrow(X)
    #-- put a names on the rows and columns of X --#
    X.names = colnames(X)
    if (is.null(X.names))
    X.names = paste("V", 1:ncol(X), sep = "")
    
    ind.names = rownames(X)
    if (is.null(ind.names))
    ind.names = 1:nrow(X)
    
    #-- ncomp
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp < 1 || !is.finite(ncomp) || ncomp>min(nr,nc))
    stop("invalid value for 'ncomp'.", call. = FALSE)
    
    ncomp = round(ncomp)
    
    #-- reconst
    if (!is.logical(reconst))
    stop("'reconst' must be a logical constant (TRUE or FALSE).",
    call. = FALSE)
    
    #-- max.iter
    if (is.null(max.iter) || max.iter < 1 || !is.finite(max.iter))
    stop("invalid value for 'max.iter'.", call. = FALSE)
    
    max.iter = round(max.iter)
    
    #-- tol
    if (is.null(tol) || tol < 0 || !is.finite(tol))
    stop("invalid value for 'tol'.", call. = FALSE)
    
    #-- end checking --#
    #------------------#
    #-- pca approach -----------------------------------------------------------#
    #---------------------------------------------------------------------------#
    
    
    #-- initialisation des matrices --#
    p = matrix(nrow = nc, ncol = ncomp)
    t.mat = matrix(nrow = nr, ncol = ncomp)
    eig = vector("numeric", length = ncomp)
    nc.ones = rep(1, nc)
    nr.ones = rep(1, nr)
    is.na.X = is.na(X)
    na.X = FALSE
    if (any(is.na.X)) na.X = TRUE
    
    #-- boucle sur h --#
    for (h in 1:ncomp)
    {
        th = X[, which.max(apply(X, 2, var, na.rm = TRUE))]
        if (any(is.na(th))) th[is.na(th)] = 0
        ph.old = rep(1 / sqrt(nc), nc)
        ph.new = vector("numeric", length = nc)
        iter = 1
        diff = 1
        
        if (na.X)
        {
            X.aux = X
            X.aux[is.na.X] = 0
        }
        
        while (diff > tol & iter <= max.iter)
        {
            if (na.X)
            {
                ph.new = crossprod(X.aux, th)
                Th = drop(th) %o% nc.ones
                Th[is.na.X] = 0
                th.cross = crossprod(Th)
                ph.new = ph.new / diag(th.cross)
            } else {
                ph.new = crossprod(X, th) / drop(crossprod(th))
            }
            
            ph.new = ph.new / drop(sqrt(crossprod(ph.new)))
            
            if (na.X)
            {
                th = X.aux %*% ph.new
                P = drop(ph.new) %o% nr.ones
                P[t(is.na.X)] = 0
                ph.cross = crossprod(P)
                th = th / diag(ph.cross)
            } else {
                th = X %*% ph.new / drop(crossprod(ph.new))
            }
            
            diff = drop(sum((ph.new - ph.old)^2, na.rm = TRUE))
            ph.old = ph.new
            iter = iter + 1
        }
        
        if (iter > max.iter)
        warning(paste("Maximum number of iterations reached for comp.", h))
        
        X = X - th %*% t(ph.new)
        p[, h] = ph.new
        t.mat[, h] = th
        eig[h] = sum(th * th, na.rm = TRUE)
    }
    
    eig = sqrt(eig)
    t.mat = scale(t.mat, center = FALSE, scale = eig)
    attr(t.mat, "scaled:scale") = NULL
    result = list(eig = eig, p = p, t = t.mat)
    
    if (reconst)
    {
        X.hat = t.mat %*% diag(eig) %*% t(p)
        
        colnames(X.hat) = colnames(X)
        rownames(X.hat) = rownames(X)
        result$rec = X.hat
    }
    
    return(invisible(result))
}
