#############################################################################################################
# Authors:
#  Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
#  Kim-Anh Le Cao, French National Institute for Agricultural Research and ARC Centre of Excellence ins Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
#  Fangzhou Yao, Queensland Facility for Advanced Bioinformatics, University of Queensland, Australia
#
# created: 2009
# last modified: 2011
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


##==========================SPARSE PCA =========================##

spca <- 
function(X, 
         ncomp = 2,
         center = TRUE, 
         scale = TRUE,
         keepX = rep(ncol(X), ncomp),
         max.iter = 500, 
         tol = 1e-06,
         logratio = 'none',# one of ('none','CLR')
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

    #-- ncomp
    if (is.null(ncomp))
    ncomp = min(nrow(X),ncol(X))

    ncomp = round(ncomp)

    if ( !is.numeric(ncomp) || ncomp < 1 || !is.finite(ncomp))
    stop("invalid value for 'ncomp'.", call. = FALSE)

    if (ncomp > min(ncol(X), nrow(X)))
    stop("use smaller 'ncomp'", call. = FALSE)

    #-- keepX
    if (length(keepX) != ncomp)
    stop("length of 'keepX' must be equal to ", ncomp, ".")
    if (any(keepX > ncol(X)))
    stop("each component of 'keepX' must be lower or equal than ", ncol(X), ".")

    #-- log.ratio
    choices = c('CLR','none')
    logratio = choices[pmatch(logratio, choices)]

    if (any(is.na(logratio)) || length(logratio) > 1)
    stop("'logratio' should be one of 'CLR'or 'none'.", call. = FALSE)

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
    X = logratio.transfo(X = X, logratio = logratio, offset = 0)#if(logratio == "ILR") {ilr.offset} else {0})
    
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


    #--scaling the data--#
    X=scale(X,center=center,scale=scale)
    cen = attr(X, "scaled:center")
    sc = attr(X, "scaled:scale")
    if (any(sc == 0))
        stop("cannot rescale a constant/zero column to unit variance.")


    #--initialization--#
    X=as.matrix(X)
    X.temp=as.matrix(X)
    n=nrow(X)
    p=ncol(X)
    
    # put a names on the rows and columns
    X.names = dimnames(X)[[2]]
    if (is.null(X.names)) X.names = paste("X", 1:p, sep = "")
    colnames(X) = X.names
    
    ind.names = dimnames(X)[[1]]
    if (is.null(ind.names)) ind.names = 1:nrow(X)
    rownames(X) = ind.names
    
    vect.varX=vector(length=ncomp)
    names(vect.varX) = paste("PC", 1:ncomp, sep = "")#c(1:ncomp)

    vect.iter=vector(length=ncomp)
    names(vect.iter) = paste("PC", 1:ncomp, sep = "")#c(1:ncomp)

    vect.keepX = keepX
    names(vect.keepX) = paste("PC", 1:ncomp, sep = "")#c(1:ncomp)

# KA: to add if biplot function (but to be fixed!)
    #sdev = vector(length = ncomp)

    mat.u=matrix(nrow=n, ncol=ncomp)
    mat.v=matrix(nrow=p, ncol=ncomp)
    colnames(mat.u)=paste("PC", 1:ncomp, sep = "")#c(1:ncomp)
    colnames(mat.v)=paste("PC", 1:ncomp, sep = "")#c(1:ncomp)
    rownames(mat.v)=colnames(X)

    #--loop on h--#
    for(h in 1:ncomp){
       
       #--computing the SVD--#
       svd.X = svd(X.temp, nu = 1, nv = 1)
       u = svd.X$u[,1]
       loadings = svd.X$v[,1]#svd.X$d[1]*svd.X$v[,1]
       v.stab = FALSE
       u.stab = FALSE
       iter = 0

       #--computing nx(degree of sparsity)--#
       nx = p-keepX[h]
       #vect.keepX[h] = keepX[h]
       
       u.old = u
       loadings.old = loadings
       #--iterations on v and u--#
       repeat{
           
            iter = iter +1
            
            loadings = t(X.temp) %*% u
            
            #--penalisation on loading vectors--#
            if (nx != 0) {
                absa = abs(loadings)
                if(any(rank(absa, ties.method = "max") <= nx)) {
                    loadings = ifelse(rank(absa, ties.method = "max") <= nx, 0, sign(loadings) * (absa - max(absa[rank(absa, ties.method = "max") <= nx])))
                }
            }
            loadings = loadings / drop(sqrt(crossprod(loadings)))
            
            u = as.vector(X.temp %*% loadings)
            #u = u/sqrt(drop(crossprod(u))) # no normalisation on purpose: to get the same $x as in pca when no keepX.
       
           #--checking convergence--#
           if(crossprod(u-u.old)<tol){break}
           if(crossprod(loadings-loadings.old)<tol){break}
           
           if (iter >= max.iter)
           {
               warning(paste("Maximum number of iterations reached for the component", h),call. = FALSE)
               break
           }
           
           u.old = u
           loadings.old = loadings
       
       }



        #v.final = v.new/sqrt(drop(crossprod(v.new)))
             
       #--deflation of data--#
       c = crossprod(X.temp, u) / drop(crossprod(u))
       X.temp= X.temp - u %*% t(c)  #svd.X$d[1] * svd.X$u[,1] %*% t(svd.X$v[,1])
       
       
       vect.iter[h] = iter
       mat.v[,h] = loadings
       mat.u[,h] = u
       
       #--calculating adjusted variances explained--#
       X.var = X %*% mat.v[,1:h]%*%solve(t(mat.v[,1:h])%*%mat.v[,1:h])%*%t(mat.v[,1:h])
       vect.varX[h] = sum(X.var^2)

# KA: to add if biplot function (but to be fixed!)
       #sdev[h] = sqrt(svd.X$d[1])
       
        
    }#fin h
    
    rownames(mat.u) = ind.names

   cl = match.call()
		cl[[1]] = as.name('spca')

    result = (list(call = cl, X = X,
		   ncomp = ncomp,	
                   #sdev = sdev,  # KA: to add if biplot function (but to be fixed!)
                   #center = center, # KA: to add if biplot function (but to be fixed!)
                   #scale = scale,   # KA: to add if biplot function (but to be fixed!)
                   varX = vect.varX/sum(X^2),
                   keepX = vect.keepX,
                   iter = vect.iter,
                   rotation = mat.v,
                   x = mat.u,
                   names = list(X = X.names, sample = ind.names),
                   loadings = list(X = mat.v),
                   variates = list(X = mat.u)
              ))
			  
    class(result) = c("spca", "prcomp", "pca")
    
    #calcul explained variance
    explX=explained_variance(X,result$variates$X,ncomp)
    result$explained_variance=explX
    
    
    return(invisible(result))
}
