#############################################################################################################
# Authors:
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Benoit Gautier, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 2013
# last modified: 05-10-2017
#
# Copyright (C) 2013
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


wrapper.sgcca = function(
X,
design = 1 - diag(length(X)),
penalty = NULL,
ncomp = 1,
keepX,
scheme = "horst",
mode = "canonical",
scale = TRUE,
init = "svd.single",
tol = .Machine$double.eps,
max.iter = 1000,
near.zero.var = FALSE,
all.outputs = TRUE
){
    
    
    check=Check.entry.sgcca(X = X, design = design ,ncomp = ncomp , scheme = scheme , scale = scale,
    init = init , tol = tol, mode = mode, max.iter = max.iter,near.zero.var = near.zero.var,keepX = keepX)
    
    
    A = check$A
    design = check$design
    ncomp = check$ncomp
    init = check$init
    scheme = check$scheme
    near.zero.var = check$near.zero.var
    keepA = check$keepA
    nzv.A = check$nzv.A
    
    keepAA = vector("list", length = max(ncomp)) # one keepA per comp
    names(keepAA) = paste0("comp",1:max(ncomp))
    for(comp in 1:max(ncomp)) # keepA[[block]] [1:ncomp]
    keepAA[[comp]] = lapply(keepA, function(x) x[comp])
    
    keepA = lapply(keepAA, expand.grid)
   
    result.sgcca = internal_mint.block(A = A, design = design, tau = NULL,
    ncomp = ncomp,
    scheme = scheme, scale = scale,
    init = init, tol = tol,
    keepA = keepA,
    max.iter = max.iter,
    study = factor(rep(1,nrow(A[[1]]))),#mint.sgcca not coded yet
    mode = mode,penalty = penalty,
    all.outputs = all.outputs
    )
    
   
    out = list(
    call = match.call(),
    X = result.sgcca$A,
    variates = result.sgcca$variates,
    loadings = result.sgcca$loadings,
    loadings.star = result.sgcca$loadings.star,
    design = result.sgcca$design,
    penalty = penalty,
    scheme = result.sgcca$scheme,
    ncomp = result.sgcca$ncomp,
    crit = result.sgcca$crit,
    AVE = result.sgcca$AVE,
    names = result.sgcca$names,#names = list(indiv = rownames(X[[1]]), var = sapply(X, colnames)),
    init = result.sgcca$init,
    tol = result.sgcca$tol,
    iter = result.sgcca$iter,
    max.iter = result.sgcca$max.iter,
    nzv = result.sgcca$nzv,
    scale = result.sgcca$scale,
    design = result.sgcca$design,
    scheme = result.sgcca$scheme,
    explained_variance = result.sgcca$explained_variance
    )
    
    class(out) = 'sgcca'
    return(invisible(out))
}

