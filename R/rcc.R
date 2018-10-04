#############################################################################################################
# Authors:
#   Sebastien Dejean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
#   Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
#   Francois Bartolo, Institut National des Sciences Appliquees et Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 2009
# last modified: 03-03-2016
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


rcc =
function(   X,
            Y,
            ncomp = 2,
            method = "ridge", #c("ridge", "shrinkage")
            lambda1 = 0,
            lambda2 = 0)
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
    
    #-- data set names --#
    data.names = c(deparse(substitute(X)), deparse(substitute(Y)))
    
    #-- method
    choices = c("ridge", "shrinkage")
    method = choices[pmatch(method, choices)]
    if (is.na(method))
    stop("'method' should be one of 'ridge' or 'shrinkage'.",
    call. = FALSE)
    
    #-- X matrix
    if (is.data.frame(X)) X = as.matrix(X)
    
    if (!is.matrix(X) || is.character(X))
    stop("'X' must be a numeric matrix.", call. = FALSE)
    
    if (any(apply(X, 1, is.infinite)))
    stop("infinite values in 'X'.", call. = FALSE)
    
    if (method == "shrinkage"){
        if (any(is.na(X)))
        stop("missing values in 'X' matrix. NAs not are allowed if method = 'shrinkage'.", call. = FALSE)
    }
    #-- Y matrix
    if (is.data.frame(Y)) Y = as.matrix(Y)
    
    if (!is.matrix(Y) || is.character(Y))
    stop("'Y' must be a numeric matrix.", call. = FALSE)
    
    if (any(apply(Y, 1, is.infinite)))
    stop("infinite values in 'Y'.", call. = FALSE)
    
    if (method == "shrinkage")
    if (any(is.na(Y)))
    stop("missing values in 'Y' matrix. NAs not are allowed if method = 'shrinkage'.", call. = FALSE)
    
    #-- equal number of rows in X and Y
    if ((n = nrow(X)) != nrow(Y))
    stop("unequal number of rows in 'X' and 'Y'.", call. = FALSE)
    
    p = ncol(X)
    q = ncol(Y)
    
    #-- put a names on the columns of X and Y --#
    X.names = colnames(X)
    if (is.null(X.names)) X.names = paste("X", 1:ncol(X), sep = "")
    
    Y.names = colnames(Y)
    if (is.null(Y.names)) Y.names = paste("Y", 1:ncol(Y), sep = "")
    
    #-- put a names on the samples --#
    ind.names = dimnames(X)[[1]]
    if (is.null(ind.names)) ind.names = dimnames(Y)[[1]]
    if (is.null(ind.names)) ind.names = 1:n
    
    #-- ncomp
    if (is.null(ncomp) || ncomp < 1 || !is.finite(ncomp))
    stop("invalid value for 'ncomp'.", call. = FALSE)
    
    ncomp = round(ncomp)
    
    if (ncomp > min(p, q))
    stop("'comp' must be smaller or equal than ", min(p, q), ".",
    call. = FALSE)
    
    #-- lambda1
    if (!is.finite(lambda1) || is.null(lambda1))
    stop("invalid value for 'lambda1'.", call. = FALSE)
    
    if(lambda1 < 0)
    stop("'lambda1' must be a non-negative value.", call. = FALSE)
    
    #-- lambda2
    if (!is.finite(lambda2) || is.null(lambda2))
    stop("invalid value for 'lambda2'.", call. = FALSE)
    
    if(lambda2 < 0)
    stop("'lambda2' must be a non-negative value.", call. = FALSE)
    
    #-- end checking --#
    #------------------#
    
    
    #-- rcc approach -----------------------------------------------------------#
    #---------------------------------------------------------------------------#
    
    #-- covariance matrices regularization --#
    if (method == "ridge") {
        Cxx = var(X, na.rm = TRUE, use = "pairwise") + diag(lambda1, ncol(X))
        Cyy = var(Y, na.rm = TRUE, use = "pairwise") + diag(lambda2, ncol(Y))
        Cxy = cov(X, Y, use = "pairwise")
    }
    else { # if method == 'shrinkage'
        Cxx = cov.shrink(X, verbose = FALSE)
        Cyy = cov.shrink(Y, verbose = FALSE)
        
        lambda.x = attr(Cxx, "lambda")
        lambda.y = attr(Cyy, "lambda")
        
        sc.x = sqrt(var.shrink(X, verbose = FALSE))
        sc.y = sqrt(var.shrink(Y, verbose = FALSE))
        
        w = rep(1/n, n)
        xs = wt.scale(X, w, center = TRUE, scale = TRUE)
        ys = wt.scale(Y, w, center = TRUE, scale = TRUE)
        
        #-- bias correction factor
        h1 = n / (n - 1)
        
        #-- unbiased empirical estimator
        Cxy = h1 * crossprod(sweep(sweep(xs, 1, sqrt((1 - lambda.x) * w), "*"), 2, sc.x, "*"),
        sweep(sweep(ys, 1, sqrt((1 - lambda.y) * w), "*"), 2, sc.y, "*"))
    }
    
    #-- calculation of the canonical correlations and canonical variables --#
    Cxx.fac = chol(Cxx)
    Cyy.fac = chol(Cyy)
    Cxx.fac.inv = solve(Cxx.fac)
    Cyy.fac.inv = solve(Cyy.fac)
    mat = t(Cxx.fac.inv) %*% Cxy %*% Cyy.fac.inv
    
    if (p >= q) {
        result = svd(mat, nu = ncomp, nv = ncomp)
        cor = result$d
        xcoef = Cxx.fac.inv %*% result$u
        ycoef = Cyy.fac.inv %*% result$v
    }
    else {
        result = svd(t(mat), nu = ncomp, nv = ncomp)
        cor = result$d
        xcoef = Cxx.fac.inv %*% result$v
        ycoef = Cyy.fac.inv %*% result$u
    }
    
    #-- output -----------------------------------------------------------------#
    #---------------------------------------------------------------------------#
    names(cor) = 1:length(cor)
    X.aux = scale(X, center = TRUE, scale = FALSE)
    Y.aux = scale(Y, center = TRUE, scale = FALSE)
    X.aux[is.na(X.aux)] = 0
    Y.aux[is.na(Y.aux)] = 0
    
    U = X.aux %*% xcoef
    V = Y.aux %*% ycoef
    
    cl = match.call()
    cl[[1]] = as.name('rcc')
    
    if (method == "ridge"){
    lambda = c("lambda1" = lambda1, "lambda2" = lambda2)
    } else { # if method == 'shrinkage')
      lambda = c("lambda1" = lambda.x, "lambda2" = lambda.y)
    }
    
    result = list(call = cl,
    X = X,
    Y = Y,
    ncomp = ncomp,
    method = method,
    cor = cor,
    loadings = list(X = xcoef, Y = ycoef),
    variates = list(X = U, Y = V),
    names = list(sample = ind.names, colnames = list(X=colnames(X),Y=colnames(Y)), blocks = c("X","Y"),#list(X = X.names, Y = Y.names, indiv = ind.names,
    data = data.names),
    lambda = lambda)
    
    #calcul explained variance
    explX=explained_variance(result$X,result$variates$X,ncomp)
    explY=explained_variance(result$Y,result$variates$Y,ncomp)
    result$explained_variance=list(X=explX,Y=explY)
    
    
    class(result) = "rcc"
    return(invisible(result))
}
