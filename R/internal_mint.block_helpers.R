################################################################################
# Author :
#   Florian Rohart,
#   Kim-Anh Le Cao,
#   Benoit Gautier,
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
# Internal helpers functions to run "mixOmics" and "internal_mint.block"
# =============================================================================

# Some of these functions have been borrowed from the RGCCA package,
#   as indicated below

# --------------------------------------
# get.weights: used in (MCV.)(block.)(s)pls(da).R
# --------------------------------------
#' Get block weights for in analyses
#'
#' weights saved in objects and used by \code{internal.predict.DA}
#'
#' @param variates Named list of variates from internal_wrapper.mint.block
#' @param indY Integer, the Y matrix index from internal_wrapper.mint.block
#'
#' @return named vector of average correlations b/w block components and Y components
#' @noRd
get.weights = function(variates, indY)
{
    ## for all X blocks, calculate correlation matrix of block components and Y components
    ## and keep the intra-component correlations only
    block_correlation_with_Y <- lapply(variates[-indY], function(x) {
        diag(cor(x, variates[[indY]]))
        })
    block_correlation_with_Y <- data.frame(t(data.frame(block_correlation_with_Y)))
    return(block_correlation_with_Y)
    
}


# --------------------------------------
# study_split: used in 'internal_mint.block.R' and 'predict.mint.block.pls.R'
# --------------------------------------
study_split = function(data, study)
{
    #data should be a matrix
    if(!is(data, "matrix"))
    data = as.matrix(data)
    
    M = length(levels(study))
    
    #---------------------- split data
    if(M>1)
    {
        data.list.study = split.data.frame(data,study)
    } else {
        data.list.study = list(data)
        names(data.list.study) = levels(study)
    }

    result = data.list.study
    return(invisible(result))
}


# --------------------------------------
# soft_thresholding: used in sparsity (below)
# --------------------------------------
# x: vector
# nx: number of entries to put to zero

soft_thresholding_L1 = function(x,nx)
{
    #selection on a (loadings.X). modified on 19/02/15 to make sure that a!=0
    if (nx!=0)
    {
        absa = abs(x)
        if (any(rank(absa, ties.method = "max") <= nx))
        {
            x = ifelse(rank(absa, ties.method = "max") <= nx, 0,
                        sign(x) * (absa - max(absa[rank(absa,
                        ties.method = "max") <= nx])))
        }
    }
    
    x
}

# -----------------------------------------------------------------------------
# soft.threshold() - soft-thresholds a vector such that
# the L1-norm constraint is satisfied.
# -----------------------------------------------------------------------------
soft.threshold = function (x, sumabs = 1)
return(soft(x, BinarySearch(x, sumabs)))

BinarySearch = function(argu,sumabs)
{
    if (norm2(argu)==0 || sum(abs(argu/norm2(argu)))<=sumabs)
    return(0)
    
    lam1 = 0
    lam2 = max(abs(argu))-1e-5
    iter = 1
    while (iter < 500)
    {
        su = soft(argu,(lam1+lam2)/2)
        if (sum(abs(su/norm2(su)))<sumabs)
        {
            lam2 = (lam1+lam2)/2
        } else {
            lam1 = (lam1+lam2)/2
        }
        if ((lam2-lam1)<1e-10)
        return((lam1+lam2)/2)
        
        iter = iter+1
    }
    warning("Didn't quite converge")
    return((lam1+lam2)/2)
}

soft = function(x,d) return(sign(x)*pmax(0, abs(x)-d))

norm2 = function(vec)
{
    a = sqrt(sum(vec^2))
    if (a == 0)
    a = .05
    
    return(a)
}


# --------------------------------------
# sparsity function: used in 'internal_mint.block.R'
# --------------------------------------
sparsity=function(loadings.A, keepA, penalty=NULL)
{
    
    if (!is.null(keepA)) {
        nx = length(loadings.A) - keepA
        loadings.A = soft_thresholding_L1(loadings.A, nx = nx)
    } else if (!is.null(penalty)) {
        loadings.A = soft.threshold(loadings.A, penalty)
    }

    return(loadings.A)
}



# --------------------------------------
# scaling with or without bias: used in mean_centering_per_study (below)
# --------------------------------------
scale.function_old=function(temp, scale = TRUE)
# problem: divide by n instead of n-#NA
{
    meanX = colMeans(temp, na.rm = TRUE)
    data.list.study.scale_i = t(t(temp) - meanX)
    if (scale)
    {
        sqrt.sdX = sqrt(colSums(data.list.study.scale_i^2, na.rm = TRUE) /
        (nrow(temp) - 1))
        data.list.study.scale_i = t(t(data.list.study.scale_i) / sqrt.sdX)
    } else {
        sqrt.sdX = NULL
    }
    
    #is.na.data = is.na(data.list.study.scale_i)
    #if (sum(is.na.data) > 0)
    #data.list.study.scale_i[is.na.data] = 0
    
    out = list(data_scale=data.list.study.scale_i, meanX=meanX,
    sqrt.sdX=sqrt.sdX)
    return(out)
}

# --------------------------------------
# scaling, using colSds from library(matrixStats),
# used in mean_centering_per_study (below)
# --------------------------------------
scale.function=function(temp, scale = TRUE)
{
    meanX = colMeans(temp, na.rm = TRUE)
    
    if (scale)
    {
        sqrt.sdX = colSds(temp,  na.rm=TRUE)
        # first possiblity: scale(), too long
        # second possibility: matrix approach: transpose is too long
        # data.list.study.scale_i = t( (t(temp)-meanX) / sqrt.sdX)
        
        # third possibility (ell.equal=>TRUE)
        data.list.study.scale_i = temp-rep(1, nrow(temp)) %*% t(meanX)
        data.list.study.scale_i = data.list.study.scale_i /
            rep(1, nrow(temp)) %*% t(sqrt.sdX)
        
        
        ind = which(sqrt.sdX == 0) # scaling can creates NA
        if(length(ind) >0)
        data.list.study.scale_i[,ind] = 0
        
    } else {
        sqrt.sdX = NULL
        # data.list.study.scale_i = t( (t(temp)-meanX)) # too long bc t()
        data.list.study.scale_i = temp-rep(1, nrow(temp)) %*% t(meanX)
    }
    
    #is.na.data = is.na(data.list.study.scale_i)
    #if (sum(is.na.data) > 0)
    #data.list.study.scale_i[is.na.data] = 0
    
    out = list(data_scale=data.list.study.scale_i, meanX=meanX,
    sqrt.sdX=sqrt.sdX)
    return(out)
}


# --------------------------------------
# Mean centering/scaling per study: used in 'internal_mint.block.R'
# --------------------------------------
mean_centering_per_study=function(data, study, scale)
{
    
    M = length(levels(study))   # number of groups
    # split the data
    data.list.study = study_split(data, study)

    # center and scale data per group, and concatene the data
    res = lapply(data.list.study, scale.function, scale = scale)
    
    meanX = lapply(res, function(x){x[[2]]})
    sqrt.sdX = lapply(res, function(x){x[[3]]})
    rownames.study = lapply(res, function(x){rownames(x[[1]])})

    #rename rows and cols of concatenated centered (and/or scaled) data
    #colnames(concat.data) = colnames(data) #already done
    
    #sort the samples as in the original X
    if(M>1) # otherwise already same order
    {
        concat.data = do.call("rbind", lapply(res,function(x){x[[1]]}))
        indice.match = match(rownames(data),rownames(concat.data))
        concat.data = concat.data[indice.match, ,drop=FALSE]
    } else{
        concat.data = res[[1]][[1]]
    }
    if (M > 1)
    {
        for (m in seq_len(M))
        {
            attr(concat.data,paste0("means:", levels(study)[m])) = meanX[[m]]
            if(scale)
            {
                attr(concat.data,paste0("sigma:", levels(study)[m])) =
                sqrt.sdX[[m]]
            } else {
                attr(concat.data,paste0("sigma:", levels(study)[m])) = NULL
            }
        }
    } else {
        attr(concat.data,"scaled:center") = meanX[[1]]
        if (scale)
        {
            attr(concat.data,"scaled:scale") = sqrt.sdX[[1]]
        } else {
            attr(concat.data,"scaled:scale") = NULL
        }
    }
    
    return(list(concat.data=concat.data, rownames.study=rownames.study))
}


# --------------------------------------
# l2.norm: used in 'internal_mint.block.R'
# --------------------------------------
l2.norm=function(x)
{
    if (!is.vector(x))
    stop("x has to be a vector")
    
    out = x / drop(sqrt(crossprod(x)))
}

# ---------------------------------------------------
# tau.estimate() - Estimation of tau accoring to Strimmer formula
# ---------------------------------------------------
#used in 'internal_mint.block.R'
tau.estimate = function (x)
{
    if (is.matrix(x) == TRUE && is.numeric(x) == FALSE)
    stop("The data matrix must be numeric!")
    
    p = NCOL(x)
    n = NROW(x)
    #covm = cov(x)
    corm = cor(x)
    xs = scale(x, center = TRUE, scale = TRUE)
    xs2 = xs^2
    v = (n/((n - 1)^3)) * (crossprod(xs2) - 1/n * (crossprod(xs))^2)
    diag(v) = 0
    m = matrix(rep(apply(xs2, 2, mean), p), p, p)
    I = diag(NCOL(x))
    d = (corm - I)^2
    tau = (sum(v))/sum(d)
    tau = max(min(tau, 1), 0)
    return(tau)
}


################################################################################
# Functions acquired from RGCCA R-library
################################################################################
# ------------------------------------------------------------------------------
# cov2() - Compute biased and unbiased covariance and variance estimates
# ------------------------------------------------------------------------------
# used in 'internal_mint.block.R'
cov2 = function (x, y = NULL, bias = FALSE) {
    n = NROW(x)
    if (is.null(y)) {
        x = as.matrix(x)
        if (bias) {
            C = ((n - 1)/n) * cov(x, use = "pairwise.complete.obs")
        } else {
            C = cov(x, use = "pairwise.complete.obs")
        }
    } else {
        if (bias) {
            C = ((n - 1)/n) * cov(x, y, use = "pairwise.complete.obs")
        } else {
            C = cov(x, y, use = "pairwise.complete.obs")
        }
    }
    return(C)
}


# ------------------------------------------------------------------------------
# miscrossprod() - Compute cross-product between vectors x and y
# ------------------------------------------------------------------------------
# used in 'internal_mint.block.R'
miscrossprod = function (x, y) {
    d.p = sum(drop(x) * drop(y), na.rm = TRUE)
    #d.p = as.vector(d.p)/norm2(d.p)     ## change made
    return(d.p)
}


# ------------------------------------------------------------------------------
# deflation()
# ------------------------------------------------------------------------------
# used in defl.select (below)
deflation = function(X, y, misdata.q, is.na.A.q, ind.NA){
    # Computation of the residual matrix R
    # Computation of the vector p.
    
    #is.na.tX <- is.na(t(X))
    #save(list=ls(),file="temp3.Rdata")
    if (misdata.q)
    {
        #is.na.tX = t(is.na.A.q)
        #p = apply(t(X),1,miscrossprod,y)/as.vector(crossprod(y))
        
        #variates.A[, q] =  apply(A[[q]], 1, miscrossprod, loadings.A[[q]])
        #A.temp = replace(t(X), is.na.tX, 0) # replace NA in A[[q]] by 0
        loadings.A.temp = crossprod(X, y)
        #temp = drop(y) %o% rep(1, ncol(A.temp.q))
        #temp[is.na.A.q] = 0
        # we only want the diagonal, which is the norm of each column of temp
        #loadings.A.norm = crossprod(temp)
        #p = variates.A.temp / diag(loadings.A.norm)
        
        #d.loadings.A.norm = apply(temp,2, crossprod)
        #only calculating the ones where there's a NA
        d.loadings.A.norm = rep(crossprod(y), ncol(X))
        #ind.NA = which(apply(is.na.A.q, 2, sum) == 1)
        
        
        if(length(ind.NA)>0)
        {
            temp = drop(y) %o% rep(1, length(ind.NA))
            # should be n*p, but we limit it to n* where there's NA
            temp[is.na.A.q[,ind.NA,drop=FALSE]] = 0
            d.loadings.A.norm[ind.NA] = apply(temp,2, crossprod)
        }
        
        p = loadings.A.temp / d.loadings.A.norm
        # we can have 0/0, so we put 0
        a = is.na(p)
        if (any(a))
        p[a] = 0
        
    } else {
        p <- crossprod(X,y) / as.vector(crossprod(y))
    }
    
    R <- X - tcrossprod(y,p)
    return(list(p=p,R=R))
}


# ------------------------------------------------------------------------------
# defl.select() - computes residual matrices
# ------------------------------------------------------------------------------
# used in 'internal_mint.block.R'
defl.select = function(yy, rr, nncomp, nn, nbloc, indY = NULL,
mode = "canonical", aa = NULL, misdata, is.na.A, ind.NA) {
    ### Start: Add new parameter for estimation classic mode
    #save(list=ls(),file="temp2.Rdata")
    resdefl = pdefl = vector("list",length=nbloc)
    for (q in seq_len(nbloc)) {
        # for each block we create missing data parameters to be passed to
        #  deflation()
        if(misdata[q])
        {
            is.na.A.q = is.na.A[[q]]
        } else {
            is.na.A.q = NULL
        }
        ### Start: insertion of new deflations
        # (See La regression PLS Theorie et pratique p204 (Chap 11))
        if ( nn <= nncomp[q] ) {
            if ((mode == "canonical") || (q != indY)) {
                #deflation of each block independently from the others,
                # except indY
                defltmp = deflation(rr[[q]], yy[ , q], misdata[q],
                is.na.A.q, ind.NA[[q]])
                #save(list=ls(),file="temp.Rdata")
                resdefl[[q]] = defltmp$R
                pdefl[[q]]   = defltmp$p
            } else if (mode == "classic") {
                resdefl[[q]] = Reduce("+", lapply(c(seq_len(nbloc))[-q], function(x)
                {rr[[q]] - yy[ ,x]%*%t(aa[[q]])}))/(nbloc-1)
                pdefl[[q]]   =  rep(0,NCOL(rr[[q]]))
            } else if (mode == "invariant") { #no deflation
                resdefl[[q]] = rr[[q]]
                pdefl[[q]]   =  rep(0,NCOL(rr[[q]]))
            } else if (mode == "regression") {
                resdefl[[q]] = Reduce("+", lapply(c(seq_len(nbloc))[-q], function(x)
                {deflation(rr[[q]],yy[, x], misdata[q], is.na.A.q,
                    ind.NA[[q]])$R}))/(nbloc-1)
                pdefl[[q]]   =  rep(0,NCOL(rr[[q]]))
            }
            ### End: insertion of new deflations
        } else {
            resdefl[[q]] = rr[[q]]
            pdefl[[q]]   =  rep(0,NCOL(rr[[q]]))
        }
    }
    names(resdefl) = names(pdefl) = names(rr)
    
    return(list(resdefl=resdefl,pdefl=pdefl))
}

# ------------------------------------------------------------------------------
# initsvd() - performs SVD on matrix X
# ------------------------------------------------------------------------------
# used in 'internal_mint.block.R'
initsvd = function (X) {
    n = NROW(X)
    p = NCOL(X)
    
    if(p>3) #use svds
    {
        ifelse(n >= p, return(svds(X, k=1, nu = 1, nv = 1)$v),
        return(svds(X, k=1, nu = 1, nv = 1)$u))
        
    } else {
        ifelse(n >= p, return(svd(X, nu = 0, nv = 1)$v),
        return(svd(X, nu = 1, nv = 0)$u))
        
    }
}

# ------------------------------------------------------------------------------
# init svd
# ------------------------------------------------------------------------------
initialisation_by_svd = function(A, indY = NULL, misdata, is.na.A = NULL,
init = "svd")
{
    
    J = length(A)
    loadings.A = vector("list",length=J)
    
    if (init == "svd")
    {
        
        # same step with or without NA, as they are already replaced by 0
        M = lapply(c(seq_len(J))[-indY], function(x){crossprod(A[[x]], A[[indY]])})
        #ssvd faster with svds, only if more than 3 column, otherwise break down
        svd.M = lapply(M, function(x){if(ncol(x)>3)
            {svds(x, k=1, nu = 1, nv = 1)} else {svd(x, nu = 1, nv = 1)}})
        
        loadings.A[c(seq_len(J))[-indY]] = lapply(seq_len(length(M)), function(x)
        {svd.M[[x]]$u})
        loadings.A[[indY]] = svd.M[[1]]$v
        
    } else if (init=="svd.single") {
        
        alpha =  lapply(seq_len(J), function(y){initsvd(A[[y]])})

        for (j in seq_len(J))
        {
            if (nrow(A[[j]]) >= ncol(A[[j]]))
            {
                loadings.A[[j]] = alpha[[j]]
            } else {
                alpha[[j]] = drop(1/sqrt( t(alpha[[j]]) %*% A[[j]] %*%
                (t(A[[j]]) %*% alpha[[j]]))) * alpha[[j]]
                
                loadings.A[[j]] = crossprod(A[[j]],alpha[[j]])
            }
        }
    } else {
        stop("init should be either 'svd' or 'svd.single'.")
    }

    return(loadings.A)

}

