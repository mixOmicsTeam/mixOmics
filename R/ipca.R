#############################################################################################################
# Authors:
#   Kim-Anh Le Cao, French National Institute for Agricultural Research and Queensland Facility for Bioinformatics, University of Queensland, Australia
#   Fangzhou Yao, Queensland Facility for Advanced Bioinformatics, University of Queensland, Australia and Shangai University of Finance and Economics, Shanghai, P.R. China
#   Jeff Coquery, Queensland Facility for Advanced Bioinformatics, University of Queensland, Australia andSup Biotech, Paris, France
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 2009
# last modified: 13-04-2016
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

# ========================================================================================================
# The function ica.par and ica.def are borrowed from the fastICA package (see references).
# ========================================================================================================


ipca = function (X,
ncomp = 2,
mode = "deflation",
fun = "logcosh",
scale = FALSE,
w.init = NULL,
max.iter = 200,
tol = 1e-04)
{
    #-- checking general input parameters --------------------------------------#
    #---------------------------------------------------------------------------#
    
    #-- X matrix
    if (is.data.frame(X)) X = as.matrix(X)
    
    if (!is.matrix(X) || is.character(X))
    stop("'X' must be a numeric matrix.", call. = FALSE)
    
    if (any(apply(X, 1, is.infinite)))
    stop("infinite values in 'X'.", call. = FALSE)
    
    if (any(is.na(X)))
    stop("missing values in 'X'.", call. = FALSE)
    
    nc = ncol(X)
    nr = nrow(X)
    
    #-- put a names on the rows and columns of X --#
    X.names = colnames(X)
    if (is.null(X.names))
    X.names = paste("V", 1:nc, sep = "")
    
    ind.names = rownames(X)
    if (is.null(ind.names))
    ind.names = 1:nr
    
    #-- ncomp
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp < 1 || !is.finite(ncomp))
    stop("invalid value for 'ncomp'.", call. = FALSE)
    
    ncomp = round(ncomp)
    
    if (ncomp > min(nc, nr))
    stop("use smaller 'ncomp'", call. = FALSE)
    
    #-- scale
    if (!is.logical(scale))
    {
        if (!is.numeric(scale) || (length(scale) != nc))
        stop("'scale' should be either a logical value or a numeric vector of length equal to the number of columns of 'X'.",
        call. = FALSE)
    }
    
    X = scale(X, center = TRUE, scale = scale)
    sc = attr(X, "scaled:scale")
    
    if (any(sc == 0))
    stop("cannot rescale a constant/zero column to unit variance.",
    call. = FALSE)
    
    #-- mode
    choices = c("deflation", "parallel")
    mode = choices[pmatch(mode, choices)]
    
    if (is.na(mode))
    stop("'mode' should be one of 'deflation' or 'parallel'.", call. = FALSE)
    
    #-- fun
    choices = c("logcosh", "exp")
    fun = choices[pmatch(fun, choices)]
    
    if (is.na(fun))
    stop("'fun' should be one of 'logcosh' or 'exp'.", call. = FALSE)
    
    #-- w.init
    if (is.null(w.init))
    {
        w.init = matrix(1 / sqrt(ncomp), ncomp, ncomp)
    } else {
        if(!is.matrix(w.init) || length(w.init) != (ncomp^2) || !is.numeric(w.init))
        stop("'w.init' is not a numeric matrix or is the wrong size", call. = FALSE)
    }
    
    if (any(is.infinite(w.init)))
    stop("infinite values in 'w.init'.", call. = FALSE)
    
    if(sum(w.init==0,na.rm=TRUE)==length(w.init))
    stop("'w.init' has to be a non-zero matrix", call. = FALSE)
    
    if(any(is.na(w.init)))
    stop("'w.init' has to be a numeric matrix matrix with non-NA values", call. = FALSE)
    
    
    #-- max.iter
    if (is.null(max.iter) || !is.numeric(max.iter) || max.iter < 1 || !is.finite(max.iter))
    stop("invalid value for 'max.iter'.", call. = FALSE)
    
    max.iter = round(max.iter)
    
    #-- tol
    if (is.null(tol) || !is.numeric(tol) || tol < 0 || !is.finite(tol))
    stop("invalid value for 'tol'.", call. = FALSE)
    
    
    #-- end checking --#
    #------------------#
    
    
    #-- ipca approach ----------------------------------------------------------#
    #---------------------------------------------------------------------------#
    
    V = svd(X)$v
    V = scale(V, center = TRUE, scale = TRUE)
    V = t(V)[1:ncomp, , drop = FALSE]
    
    if (mode == "deflation")
    {
        W = ica.def(V, ncomp = ncomp, tol = tol, fun = fun, alpha = 1,
        max.iter = max.iter, verbose = FALSE, w.init = w.init)
    }else if (mode == "parallel") {
        W = ica.par(V, ncomp = ncomp, tol = tol, fun = fun, alpha = 1,
        max.iter = max.iter, verbose = FALSE, w.init = w.init)
    }
    
    #-- independent loadings --#
    S = matrix(W %*% V, nrow = ncomp)
    
    #-- order independent loadings by kurtosis --#
    kurt = apply(S, 1, function(x) { n = length(x)
        x = x - mean(x)
        n * sum(x^4) / (sum(x^2)^2) - 3 } )
    ord = order(kurt, decreasing = TRUE)
    kurt = kurt[ord]
    S = S[ord, , drop = FALSE]
    norm = apply(S, 1, function(x) { crossprod(x) })
    S = t(sweep(S, 1, norm, "/"))
    
    #-- independent PCs / force orthonormality --#
    ipc = matrix(nrow = nr, ncol = ncomp)
    ipc[, 1] = X %*% S[, 1]
    ipc[, 1] = ipc[, 1] / sqrt(crossprod(ipc[, 1]))
    
    if (ncomp > 1)
    {
        for (h in 2:ncomp)
        {
            ipc[, h] = lsfit(y = X %*% S[, h], ipc[, 1:(h - 1)], intercept = FALSE)$res
            ipc[, h] = ipc[, h] / sqrt(crossprod(ipc[, h]))
        }
    }
    
    
    #-- output -----------------------------------------------------------------#
    #---------------------------------------------------------------------------#
    dimnames(S) = list(X.names, paste("IPC", 1:ncol(S), sep = ""))
    dimnames(ipc) = list(ind.names, paste("IPC", 1:ncol(ipc), sep = ""))
    
    cl = match.call()
    cl[[1]] = as.name('ipca')
    
    result = list(call = cl, X = X, ncomp = ncomp, x = ipc,
    loadings=list(X=S),rotation=S,variates=list(X=ipc), kurtosis = kurt, unmixing = t(W),
    mixing = t(t(W) %*% solve(W %*% t(W))),
    names = list(var = X.names, sample = ind.names))
    
    class(result) = c("ipca", "pca")
    
    #calcul explained variance
    explX=explained_variance(X,result$variates$X,ncomp)
    result$explained_variance=explX
    
    return(invisible(result))
}
