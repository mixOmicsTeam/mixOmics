################################################################################
# Authors:
#   Amrit Singh,
#   Florian Rohart,
#   Kim-Anh Le Cao,
#
# created: 2015
# last modified: 19-08-2016
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



################################################
#
## 2) cim_diablo
#
################################################

# This function is a small wrapper of cim.
# For more customisation, please use cim

cimDiablo = function(object,
color = NULL,
color.Y,
color.blocks,
comp = NULL,
margins = c(2, 15),
legend.position="topright",
transpose = FALSE,
row.names = TRUE,
col.names = TRUE,
size.legend=1.5)
{
    
    # check input object
    if (!any(class(object) == "block.splsda"))
    stop("cimDiablo is only available for 'block.splsda' objects")

    if (length(object$X) <= 1)
    stop("This function is only available when there are more than 3 blocks")
    # so 2 blocks in X + the outcome Y
    
    ncomp = min(object$ncomp)
    #-- comp
    if(is.null(comp))
    {comp = 1:ncomp}
    if (length(comp) > 1) {
        comp = unique(comp)
        if (!is.numeric(comp) || any(comp < 1))
        stop("invalid vector for 'comp'.", call. = FALSE)
        if (any(comp > ncomp))
        stop("the elements of 'comp' must be smaller or equal than ", ncomp,
        ".", call. = FALSE)
    }
    
    if (length(comp) == 1) {
        if (is.null(comp) || !is.numeric(comp) || comp <= 0 || comp > ncomp)
        stop("invalid value for 'comp'.", call. = FALSE)
        comp=c(comp,comp)
    }
    
    comp = round(comp)
    
    # color
    if(missing(color.Y))
    {
        color.Y = color.mixo(1:nlevels(object$Y))
    } else {
        if(length(color.Y) != nlevels(object$Y))
        stop("'color.Y' needs to be of length ", nlevels(object$Y))
    }
    
    
    if(missing(color.blocks))
    {
        color.blocks = brewer.pal(n = 12, name = 'Paired')[seq(2, 12, by = 2)]
    } else {
        if(length(color.blocks) != length(object$X))
        stop("'color.blocks' needs to be of length ", length(object$X))
    }
    
    X = object$X
    Y = object$Y

    #need to reorder variates and loadings to put 'Y' in last
    indY = object$indY
    object$variates = c(object$variates[-indY], object$variates[indY])
    object$loadings = c(object$loadings[-indY], object$loadings[indY])
    
    #reducing loadings for ncomp
    object$loadings = lapply(object$loadings, function(x){x[, comp,
        drop=FALSE]})
    
    keepA = lapply(object$loadings, function(i) apply(abs(i), 1, sum) > 0)
    XDatList = mapply(function(x, y){
        x[, y]
    }, x=X, y=keepA[-length(keepA)], SIMPLIFY=FALSE)
    XDat = do.call(cbind, XDatList)
    XDat[which(XDat > 2)] = 2
    XDat[which(XDat < -2)] = -2
    
    #dark = brewer.pal(n = 12, name = 'Paired')[seq(2, 12, by = 2)]
    VarLabels = factor(rep(names(X), lapply(keepA[-length(keepA)], sum)),
    levels = names(X))#[order(names(X))])
    
    ## Plot heatmap
    opar = par()[! names(par()) %in% c("cin", "cra", "csi", "cxy", "din",
        "page")]
    par(mfrow=c(1,1))
    cim(XDat,transpose= transpose, color = color,
    row.names = row.names, col.names = col.names,
    col.sideColors = color.blocks[as.numeric(VarLabels)],
    row.sideColors = color.Y[as.numeric(Y)], margins = margins)
    
    if(!transpose)
    {
        legend(legend.position,
        c("Rows", c(levels(Y)[order(levels(Y))], "", "Columns", names(X))),
        col = c(1, color.Y, 1, 1,
        color.blocks[1:nlevels(VarLabels)][match(levels(VarLabels), names(X))]),
        pch = c(NA, rep(19, nlevels(Y)), NA, NA, rep(19, nlevels(VarLabels))),
        bty="n",
        cex = size.legend,
        text.font = c(2, rep(1, nlevels(Y)), NA, 2, rep(1, nlevels(VarLabels))))

    } else { # if transpose == TRUE, rows and columns must be switched
        legend(legend.position,
        c("Rows", names(X), "", "Columns", c(levels(Y)[order(levels(Y))])),
        col = c(1, color.blocks[1:nlevels(VarLabels)][match(levels(VarLabels),
        names(X))], 1, 1, color.Y),
        pch = c(NA, rep(19, nlevels(VarLabels)), NA, NA, rep(19, nlevels(Y))),
        bty="n",
        cex = size.legend,
        text.font = c(2, rep(1, nlevels(VarLabels)), NA, 2, rep(1, nlevels(Y))))
        
    }
    par(opar)
    
    return(invisible(XDat))
}

