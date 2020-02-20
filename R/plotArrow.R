############################################################################################################
# Authors:
#   Francois Bartolo, Institut National des Sciences Appliquees et Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
#   Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 2015
# last modified: 24-05-2016
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


#----------------------------------------------------------------------------------------------------------#
#-- Includes plotArrow for PLS, sPLS, rCC,  rGCCA, sGCCA, sGCCDA --#
#----------------------------------------------------------------------------------------------------------#

plotArrow <-
function(object,
comp = NULL,
abline = FALSE,
xlim = NULL,
ylim = NULL,
group = NULL,
col,
cex,
pch,
title = NULL,
plot.arrows = TRUE,
legend = FALSE,
X.label = NULL,
Y.label = NULL,
ind.names = FALSE,
position.names = 'centroid'
)
{
    
    class.object = class(object)
    object.pls=c("mixo_pls", "mixo_plsda", "mixo_spls", "smixo_plsda", "rcc")
    object.blocks=c("sgcca", "sgccda", "rgcca")
    
    if (! any(class.object %in% c(object.pls,object.blocks)))
    stop( " 'plotArrow' is only implemented for the following objects: pls, plsda, spls, splsda, rcc, sgcca, sgccda, rgcca", call.=FALSE)
    
    
    ### Start: Validation of arguments
    ncomp = object$ncomp
    
    if ((!position.names %in% c("centroid", "start", "end")))
    stop("'position.names' must be one of 'centroid', 'start' , 'end' .", call. = FALSE)
    
    
    if (any(class.object %in% object.blocks))
    {
        
        blocks = object$names$blocks
        
        #if (class.object[1] == "sgccda")
        # blocks = blocks[-object$indY]
        object$variates = object$variates[names(object$variates) %in% blocks]
        
        if (any(object$ncomp[blocks] == 1))
        stop(paste("The number of components for one selected block '", paste(blocks, collapse = " - "),"' is 1. The number of components must be superior or equal to 2."), call. = FALSE)
        
        ncomp = object$ncomp[blocks]
    }
    
    #-- arrows
    if (!is.logical(plot.arrows))
    stop("'plot.arrows' must be a logical constant (TRUE or FALSE).", call. = FALSE)
    
    #-- xlim,ylim
    lim.X = FALSE
    if (!is.null(xlim))
    {
        if (!is.numeric(xlim) || length(xlim) != 2)
        stop("'xlim' must be a vector of length 2.",call. = FALSE)
        
        lim.X = TRUE
    }
    
    lim.Y = FALSE
    if (!is.null(ylim))
    {
        if (!is.numeric(ylim) || length(ylim) != 2)
        stop("'ylim' must be a vector  of length 2.",call. = FALSE)
        
        lim.Y = TRUE
    }
    
    
    #-- comp
    if (is.null(comp))
    comp=c(1:2)
    if (length(comp) != 2)
    stop("'comp' must be a numeric vector of length 2.", call. = FALSE)
    
    if (!is.numeric(comp))
    stop("Invalid vector for 'comp'.")
    
    if (any(ncomp < max(comp)))
    stop("Each element of 'comp' must be smaller or equal than ", max(object$ncomp), ".", call. = FALSE)
    
    comp1 = round(comp[1])
    comp2 = round(comp[2])
    
    
    if (any(class.object %in% object.pls))
    blocks=c("X","Y");object$variates = object$variates[names(object$variates) %in% blocks]
    
    if (is.null(X.label))
        X.label = sprintf("Dimension %s", comp1)

    if (is.null(Y.label))
        Y.label = sprintf("Dimension %s", comp2)
    
    
    
    #-- ind.names
    display.names.start = FALSE
    display.centroid = FALSE
    display.names.end = FALSE
    
    if (isTRUE(ind.names))
    {
        ind.names = object$names$sample
        
        if ( position.names=="centroid")
        {
            display.centroid = TRUE
        } else if (position.names=="start") {
            display.names.start =TRUE
        } else if (position.names=="end") {
            display.names.end = TRUE
        }
        
    } else if (length(ind.names) > 1 ) {
        if (length(ind.names) != length(object$names$sample))
        stop("'ind.names' must be a character vector of length ", length(object$names$sample), " or boolean  .")
        
        if ( position.names=="centroid"){
            display.centroid = TRUE
        } else if (position.names=="start"){
            display.names.start =TRUE
        } else if (position.names=="end") {
            display.names.end = TRUE
        }
    }
    
    #-- Start: Retrieve variates from object
    x = y =list()
    
    x = lapply(object$variates, function(x){x[, comp1, drop = FALSE]})
    y = lapply(object$variates, function(x){x[, comp2, drop = FALSE]})
    
    
    
    
    #-- End: Retrieve variates from object
    
    
    
    #-- Define group
    missing.group = FALSE
    if (is.null(group) & any(class.object %in% c("DA")))
    {
        group = object$Y#factor(map(object$ind.mat), labels = object$names$colnames$Y)
    } else if (!is.null(group)) {
        missing.group = TRUE
        if (!is.factor(group))
        group = as.factor(group)

        object$ind.mat = unmap(group)
        
        if (length(group) != length(x[[1]]))
        stop("Length of 'group' should be of length ", length(x[[1]]), ", the sample size of your data")
    } else {
        group = factor(rep("No group", length(x[[1]])))
        object$ind.mat = unmap(group)
    }
    
    #-- col.per.group argument
    if (nlevels(group) < 10)
    {   #only 10 colors in color.mixo
        col.per.group = color.mixo(1:nlevels(group))
    } else {
        #use color.jet
        col.per.group = color.jet(nlevels(group))
    }
    
    levels.color = vector(, length(x[[1]]))
    if (length(col.per.group) != length(x[[1]]))
    {
        for (i in 1 : nlevels(group))
        levels.color[group == levels(group)[i]] = col.per.group[i]
    } else {
        levels.color = col.per.group
    }
    
    
    #-- col argument
    missing.col = FALSE
    if (!missing(col))
    {
        if (length(col) > length(x[[1]]))
        stop("Length of 'col' should be of length inferior or equal to ", length(x[[1]]),".")
        
        col = factor(rep(col, ceiling(length(x[[1]])/length(col)))[1 : length(x[[1]])])
        if (!missing.group)
        {
            group = col
            levels.color = col
            col.per.group = levels(col)
            object$ind.mat = unmap(group)
        }
        missing.col = TRUE
    } else {
        col = levels.color
    }
    
    #-- cex argument
    if (missing(cex))
    {
        cex = rep(1, length(x[[1]]))
    } else {
        if (length(cex) == 1)
        {
            cex = rep(cex, length(x[[1]]))
        } else if (length(cex) > length(x[[1]])) {
            stop("Length of 'cex' should be of length inferior or equal to ", length(x[[1]]),".")
        } else if (length(cex) == length(unique(group))){
            cex = cex[as.factor(group)]
        }else {
            cex = rep(cex, ceiling(length(x[[1]])/length(cex)))[1 : length(x[[1]])]
        }
    }
    
    #-- pch argument
    if (missing(pch))
    {
        if (missing.col)
        {
            pch = as.numeric(col)
        } else {
            pch = as.numeric(group)
        }
    } else {
        if (length(pch) == 1)
        {
            pch = rep(pch, length(x[[1]]))
        } else if (length(pch) > length(x[[1]])){
            stop("Length of 'pch' should be of length inferior or equal to ", length(group),".")
        } else if (length(pch) == length(unique(group))){
            pch = pch[as.factor(group)]
        } else {
            pch = rep(pch, ceiling(length(x[[1]])/length(pch)))[1 : length(x[[1]])]
        }
    }
    
    
    arrow1 = function(x0, y0, x1, y1, length = 0.12, angle = 15, color)
    {
        
        d0 <- sqrt((x0 - x1)^2 + (y0 - y1)^2)
        
        if (d0 < 1e-07)
        return(invisible())
        
        arrows(x0, y0, x1, y1, angle = angle, length = max(c(length, length * cex)),
        col = color)
    }  #fin function arrow
    
    
    #-- Start: data set
    df = list()
    for (i in 1 : length(x))
    df[[i]] = data.frame(x = x[[i]], y = y[[i]], group = group)
    
    df = data.frame(do.call(rbind, df), "Block" = paste0("Block: ", unlist(lapply(1 : length(df), function(z){rep(blocks[z], nrow(df[[z]]))}))))
    
    
    names(df)[1:2] = c("x", "y")
    
    if (display.names.start ||display.names.end || display.centroid)
    df$names = rep(ind.names, length(x))
    
    df$pch = pch; df$cex = cex; df$col = as.character(col)
    
    
    
    opar <- par()[! names(par()) %in% c("cin", "cra", "csi", "cxy", "din", "page")]
    #-- Define layout
    if (legend)
    {
        par(mai=c( 1.360000, 1.093333, 1.093333,(max(strwidth(levels(group),"inches")))+0.6),xpd=TRUE)
    } else {
        par(mar=c(5,4,4,2))
    }
    
    plot(df[, "x" ], df[, "y" ],
    type = "n", xlab = X.label, ylab = Y.label, main = title,
    xlim = xlim, ylim = ylim)
    
    #-- arrows --#
    
    for (j in 1 : length(df[ df$Block %in% paste0("Block: ", blocks[1]), "x"]))
    {
        if (length(blocks)==2)
        {
            x0=df[df$Block %in% paste0("Block: ", blocks[1]), "x"][j]
            y0=df[df$Block %in% paste0("Block: ", blocks[1]), "y"][j]
            x1=df[df$Block %in% paste0("Block: ", blocks[2]), "x"][j]
            y1=df[df$Block %in% paste0("Block: ", blocks[2]), "y"][j]
            d0 <- sqrt((x0 - x1)^2 + (y0 - y1)^2)
            if (d0 < 1e-07)
            return(invisible())
            
            if (plot.arrows)
            arrows(x0, y0, x1, y1, col=df[, "col"][j],angle = 15, length = max(c(0.12, 0.12 * df$cex[j])),xpd=FALSE)
            
            if (display.centroid)
            {
                x2=mean(c(x0,x1))
                y2=mean(c(y0,y1))
                text(x2, y2, df[df$Block %in% paste0("Block: ", blocks[1]), "names"][j],col = df[df$Block %in% paste0("Block: ", blocks[1]), "col"][j], cex = df$cex[j],xpd=FALSE)
            }
        } else if (length(blocks)>2) {
            if (display.names.start)
            {
                display.names.start=FALSE
                display.centroid=TRUE
            }
            x0=y0=0
            x=y=c()
            for (i in 1 : length(blocks))
            {
                x0=x0+df[df$Block %in% paste0("Block: ", blocks[i]), "x"][j]
                y0=y0+df[df$Block %in% paste0("Block: ", blocks[i]), "y"][j]
                x[i]=df[df$Block %in% paste0("Block: ", blocks[i]), "x"][j]
                y[i]=df[df$Block %in% paste0("Block: ", blocks[i]), "y"][j]
            }
            x0=x0/length(blocks)
            y0=y0/length(blocks)
            
            if (display.centroid)
            {
                text(x0, y0, df[df$Block %in% paste0("Block: ", blocks[1]), "names"][j],col = df[df$Block %in% paste0("Block: ", blocks[1]), "col"][j], cex = df$cex[j],xpd=FALSE)
            } else{
                points(cbind(x0,y0),pch=8,cex=df$cex[j],col=df[df$Block %in% paste0("Block: ", blocks[1]), "col"][j],xpd=FALSE)
            }
            if (plot.arrows)
            {
                for (i in 1 : length(blocks))
                arrows(x0, y0, x[i], y[i], col=df[, "col"][j],angle = 15, length = max(c(0.12, 0.12 * df$cex[j])),xpd=FALSE)
            }
        }
    }
    
    names.end.blocks=TRUE
    if (display.names.start)
    {
        print.names=paste0("Block: ", blocks[1])
        print.points=paste0("Block: ", blocks[2])
    } else if (display.names.end) {
        if (length(blocks)==2)
        {
            print.names=paste0("Block: ", blocks[2])
            print.points=paste0("Block: ", blocks[1])
        } else {
            names.end.blocks=FALSE
            print.names=paste0("Block: ", blocks[1:length(blocks)])
        }
    } else {
        print.points=paste0("Block: ", blocks[1:length(blocks)])
    }
    
    
    
    #-- color samples according to col
    for (i in unique(col))
    {
        if (display.names.end || display.names.start)
        text(x = df[df$Block %in% print.names, "x"],
        y = df[df$Block %in% print.names, "y"],
        labels = df[df$Block %in% print.names, "names"],
        col = df[df$Block %in% print.names, ]$col, cex = df[df$Block %in% print.names, ]$cex,xpd=FALSE)
        if (names.end.blocks)
        points(x = df[df$Block %in% print.points, "x"],
        y = df[df$Block %in% print.points , "y"],
        col = df[df$Block %in% print.points, ]$col, cex = df[df$Block %in% print.points, ]$cex, pch = df[df$Block %in% print.points, ]$pch,xpd=FALSE)
    }
    
    
    #-- Abline
    if (abline)
    abline(v = 0, h = 0, lty = 1,xpd=FALSE)
    
    
    pch.legend=NULL
    if (missing.col)
    {
        for (i in 1:nlevels(factor(col)))
        pch.legend=c(pch.legend,df[df$col == levels(factor(col))[i], ]$pch)
    } else {
        for (i in 1:nlevels(group))
        pch.legend=c(pch.legend,df[df$group == levels(group)[i], ]$pch)
    }
    
    if (legend && (any(class.object %in% c("sgccda", "DA"))||missing.group))
        legend(par()$usr[2]+0.1,par()$usr[4] - (par()$usr[4]-par()$usr[3])/2, col = col.per.group, legend = levels(group), pch = if (display.names.end || display.names.start) {16} else {unique(pch.legend)}, title = 'Legend', cex = 0.8)
    
    opar["usr"]=par()["usr"]
    
    par(opar)
    
    return(invisible(df))
    
}
