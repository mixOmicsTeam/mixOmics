#############################################################################################################
# Author :
#   Amrit Singh. University of British Columbia, Vancouver, Canada.
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 04-2015
# last modified: 25-08-2016
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
# functions
#     1) plotIndiv_diablo; 2 sub-functions; splotMatPlot() and panel.ellipses
#     2) circosPlot_diablo
#     3) heatmap_diablo
#     4) enrichPathwayNetwork_diablo
#
# ========================================================================================================


plot.sgccda = plotDiablo = function(x,
ncomp = 1,
legend = TRUE,
legend.ncol,
...)
{
    
    object=x
    #need to reorder variates and loadings to put 'Y' in last
    opar = par()[! names(par()) %in% c("cin", "cra", "csi", "cxy", "din", "page")]

    indY=object$indY
    object$variates=c(object$variates[-indY],object$variates[indY])
    object$loadings=c(object$loadings[-indY],object$loadings[indY])
    
    VarX = do.call(cbind, lapply(object$variates, function(i) i[, ncomp]))
    datNames = colnames(VarX)
    
    if(ncol(VarX)<=2)
    stop("This function is only available when there are more than 3 blocks") # so 2 blocks + the outcome Y
    
    # check input parameters
    Y=object$Y
    
    if (length(ncomp) != 1 | ncomp > min(object$ncomp))
    stop(paste0("'ncomp' must be a numeric value lower than ", min(object$ncomp),", which is min(object$ncomp)"))
    # end check parameters
    
    if(missing(legend.ncol))
    legend.ncol = min(5, nlevels(Y))

    numberOfCols = ncol(VarX)-1
    numberOfRows = numberOfCols #- 1
    
    mat = matrix(0, nrow = numberOfRows, ncol = numberOfRows)
    for(i in 1:nrow(mat))
    {
        for(j in 1:ncol(mat))
        mat[i,j] = paste(i,j, sep="_")
    }
    plotType = list(cor=mat[lower.tri(mat)], scatter=mat[upper.tri(mat)],
    lab=diag(mat))#,
    #bar=paste(1:(numberOfRows-1), numberOfCols, sep="_"),
    #stackedbar=paste(numberOfRows, numberOfCols, sep="_"))
    
    par(mfrow = c(numberOfRows+1, numberOfCols), mar = rep.int(1/2, 4), oma = c(2,2,2,2))
    layout(matrix(c(1:(numberOfCols)^2, rep((numberOfCols)^2+1,numberOfCols)),numberOfRows+1,numberOfCols, byrow=TRUE),
    heights = c(rep(1,numberOfRows), 0.25 * floor(nlevels(Y)/legend.ncol)))
    for(i in 1:numberOfRows)
    {
        for(j in 1:numberOfCols)
        {
            ptype = unlist(lapply(plotType, function(x)
            {
                intersect(paste(i,j,sep="_"), x)
            }))
            splotMatPlot(x=VarX[, j], y=VarX[, i], datNames, Y, ptype)
            
            if(i == 1 & j %in% seq(2, numberOfRows, 1))
            Axis(side = 3, x=VarX[, i])
            
            if(j == numberOfRows & i %in% seq(1, numberOfRows-1, 1))
            Axis(side = 4, x=VarX[, i])
        }
    }
    #add legend
    plot(1:3,1:3,type="n",axes=FALSE,xlab="",ylab="")
    if(legend)
    legend("center",legend=levels(Y), col = color.mixo(1:nlevels(Y)), pch = 19, ncol = legend.ncol, cex = 1.5)
    
    par(opar)
}


splotMatPlot = function(x, y, datNames, Y, ptype)
{
    if(names(ptype) == "cor")
    {
        plot(1, type = "n", axes = FALSE)
        r = round(cor(x, y), 2)
        text(1, 1, labels=r, cex = 0.6/strwidth(abs(r))*abs(r))
        box()
    }
    if(names(ptype) == "scatter")
    panel.ellipses(x=x, y=y, Y = Y)
    
    if(names(ptype) == "lab")
    {
        plot(1, type = "n", axes = FALSE)
        ind = as.numeric(unlist(lapply(strsplit(ptype, "_"), unique)))
        text(x=1, y=1, labels=datNames[ind], cex = 2)
        box()
    }
    if(FALSE)
    {
    if(names(ptype) == "bar")
    {
        Y2 = factor(as.character(Y), levels = groupOrder)
        boxplot(x ~ Y2, horizontal=TRUE, axes = FALSE, ylim = c(min(x)-3, max(x)),
        col= color.mixo(match(levels(Y2), levels(Y))))
        axis(4, at=1:nlevels(Y2), labels=levels(Y2))
    }
    if(names(ptype) == "stackedbar")
    {
        Y2 = factor(as.character(Y), levels = groupOrder)
        bars = table(Y2)
        barplot(bars, col= color.mixo(match(levels(Y2), levels(Y))),
        axes = FALSE)
        axis(4, at=seq(0,max(bars),length.out=5), labels=seq(0,max(bars),length.out=5))
    }
    }
}

panel.ellipses = function(x, y, Y = Y, pch = par("pch"), col.lm = "red", axes = FALSE, ...)
{
    ind.gp = matrice = cdg = variance = list()
    for(i in 1:nlevels(Y))
    ind.gp[[i]] = which(as.numeric(Y)==i)
    
    matrice = lapply(ind.gp, function(z){matrix(c(x[z], y[z]), ncol = 2)})
    cdg = lapply(matrice, colMeans)
    variance = lapply(matrice, var)
    
    #library(ellipse)
    coord.ellipse = lapply(1:nlevels(Y), function(x){ellipse(variance[[x]], centre = cdg[[x]], ellipse.level = 0.95)})
    max.ellipse = sapply(coord.ellipse, function(x){apply(x, 2, max)})
    min.ellipse = sapply(coord.ellipse, function(x){apply(x, 2, min)})
    ind.names = names(Y)
    cex = 0.5
    
    plot(x, y, xlab = "X.label", ylab = "Y.label", col=color.mixo(as.numeric(Y)), pch=20, axes=axes,
    xlim = c(min(x, min.ellipse[1, ]), max(x, max.ellipse[1, ])), ylim = c(min(y, min.ellipse[2, ]), max(y, max.ellipse[2, ])))
    #text(x, y, ind.names, col = col, cex = cex)
    box()
    for (z in 1:nlevels(Y))
    points(coord.ellipse[[z]], type = "l", col = color.mixo(c(1:nlevels(Y))[z]))
    
}

