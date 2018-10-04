#############################################################################################################
# Authors:
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Benoit Gautier, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Francois Bartolo, Institut National des Sciences Appliquees et Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
#   Kim-Anh Le Cao, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 16-03-2016
# last modified: 24-08-2016
#
# Copyright (C) 2016
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
#-- Includes plotIndiv rGCCA, sGCCA, sGCCDA --#
#----------------------------------------------------------------------------------------------------------#

plotIndiv.sgcca = 
plotIndiv.rgcca = 
#plotIndiv.sgccda =   #because sgcca

function(object, 
comp = NULL, 
blocks = NULL, # to choose which block data to plot, when using GCCA module
ind.names = TRUE, 
group, # factor indicating the group membership for each sample, useful for ellipse plots. Coded as default for the -da methods, but needs to be input for the unsupervised methods (PCA, IPCA...)
col.per.group, 
style = "ggplot2", # can choose between graphics, 3d, lattice or ggplot2
ellipse = FALSE, 
ellipse.level = 0.95, 
centroid = FALSE, 
star = FALSE, 
title = NULL, 
subtitle, 
legend = FALSE, 
X.label = NULL, 
Y.label = NULL, 
Z.label = NULL, 
abline = FALSE, 
xlim = NULL, 
ylim = NULL, 
col, 
cex, 
pch, 
pch.levels,
alpha = 0.2,
axes.box = "box", 
layout = NULL, 
size.title = rel(2), 
size.subtitle = rel(1.5), 
size.xlabel = rel(1), 
size.ylabel = rel(1), 
size.axis = rel(0.8), 
size.legend = rel(1), 
size.legend.title = rel(1.1), 
legend.title = "Legend",
legend.title.pch = "Legend",
legend.position = "right",
point.lwd = 1, 
...
)
{
    plot_parameters = list(size.title = size.title, size.subtitle = size.subtitle, size.xlabel = size.xlabel, size.ylabel = size.ylabel,
    size.axis = size.axis, size.legend = size.legend, size.legend.title = size.legend.title, legend.title = legend.title,
    legend.title.pch = legend.title.pch, legend.position = legend.position, point.lwd = point.lwd, alpha = alpha)

    if(any(class(object)%in%c("mint.block.pls", "mint.block.spls", "mint.block.plsda", "mint.block.splsda")))
    stop("No plotIndiv for the following functions at this stage: mint.block.pls, mint.block.spls, mint.block.plsda, mint.block.splsda.")
    
    #-- rep.space
    rep.space = "multi" # rep.space is not used afterwards, put to "multi" to plot all blocks
    
    
    if (is.null(blocks))
    {
        blocks = names(object$X)#names$blocks
    } else if (is.numeric(blocks) & min(blocks) > 0 &  max(blocks) <=  length(object$names$blocks)) {
        blocks = object$names$blocks[blocks]
    } else if (is.character(blocks)) {
        if (!any(blocks %in% object$names$blocks))
        stop("One element of 'blocks' does not match with the names of the blocks")
    } else {
        stop("Incorrect value for 'blocks'", call. = FALSE)
    }
    #object$variates = object$variates[names(object$variates) %in% blocks] # reduce the variate to the 'blocks' we are looking at
    object$variates = object$variates[match(blocks, names(object$variates))] # reduce the variate to the 'blocks' we are looking at
    
    if (any(object$ncomp[blocks] ==  1))
    stop(paste("The number of components for one selected block '", paste(blocks, collapse = " - "), "' is 1. The number of components must be superior or equal to 2."), call. = FALSE)
    ncomp = object$ncomp[blocks]

    
    if(length(blocks)!= length(unique(blocks)))
    stop("Duplicate in 'blocks' not allowed")

    if(!missing(subtitle))
    {
        if(length(subtitle)!= length(blocks) | length(subtitle)!= length(unique(subtitle)))
        stop("'subtitle' indicates the subtitle of the plot for each 'blocks'; it needs to be the same length as 'blocks' and duplicate are not allowed.")
    }



    #-- check inputs
    check = check.input.plotIndiv(object = object, comp = comp, blocks = blocks, ind.names = ind.names, 
    style = style, ellipse = ellipse, ellipse.level = ellipse.level, centroid = centroid, 
    star = star, legend = legend, X.label = X.label, Y.label = Y.label, Z.label = Z.label, abline = abline, 
    xlim = xlim, ylim = ylim, alpha = alpha, axes.box = axes.box, plot_parameters = plot_parameters)
    #-- retrieve outputs from the checks
    axes.box = check$axes.box
    comp = check$comp
    xlim = check$xlim
    ylim = check$ylim
    ind.names = check$ind.names
    display.names = check$display.names


    #-- get the variates
    variate = internal_getVariatesAndLabels(object, comp, blocks = blocks, style = style, X.label = X.label, Y.label = Y.label, Z.label = Z.label, rep.space = rep.space)
    #-- retrieve outputs
    x = variate$x
    y = variate$y
    z = variate$z
    X.label = variate$X.label
    Y.label = variate$Y.label
    Z.label = variate$Z.label

    n = nrow(object$X[[1]])

    # create data frame df that contains (almost) all the ploting information
    out = shape.input.plotIndiv(object = object, n = n, blocks = blocks, x = x, y = y, z = z, ind.names = ind.names, group, col.per.group = col.per.group, 
    style = style, study = "global", ellipse = ellipse, ellipse.level = ellipse.level, 
    centroid = centroid, star = star, title = title, xlim = xlim, ylim = ylim, 
    col = col, cex = cex, pch = pch, pch.levels = pch.levels, display.names = display.names, plot_parameters = plot_parameters)
    #-- retrieve outputs
    df = out$df
    df.ellipse = out$df.ellipse
    col.per.group = out$col.per.group
    title = out$title
    display.names = out$display.names
    xlim = out$xlim
    ylim = out$ylim
    #missing.col = out$missing.col
    ellipse = out$ellipse
    centroid = out$centroid
    star = out$star
    plot_parameters = out$plot_parameters


    # change the levels of df.final$Block to "subtitle"
    if(!missing(subtitle))
    {
        df$Block = factor(df$Block, labels = subtitle)
        if(ellipse)
        df.ellipse$Block = factor(df.ellipse$Block, labels = subtitle)
    }
    
    #call plot module (ggplot2, lattice, graphics, 3d)
    res = internal_graphicModule(df = df, centroid = centroid, col.per.group = col.per.group, title = title, X.label = X.label,
    Y.label = Y.label, Z.label = Z.label, xlim = xlim, ylim = ylim, class.object = class(object), display.names = display.names, legend = legend,
    abline = abline, star = star, ellipse = ellipse, df.ellipse = df.ellipse, style = style, layout = layout,
    #missing.col = missing.col,
    axes.box = axes.box, plot_parameters = plot_parameters, alpha = alpha)

    return(invisible(list(df = df, df.ellipse = df.ellipse, graph = res)))


}
