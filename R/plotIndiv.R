#############################################################################################################
# Authors:
#   Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Benoit Gautier, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Francois Bartolo, Institut National des Sciences Appliquees et Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD

#
# created: 2009
# last modified: 24-08-2016
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


plotIndiv  =
function(object, ...) UseMethod("plotIndiv")


# --------------------------------------------------------------------------------------
# Internal helpers functions to run "plotIndiv" functions
# --------------------------------------------------------------------------------------


check.input.plotIndiv = function(object,
comp = NULL,
blocks = NULL, # to choose which block data to plot, when using GCCA module
ind.names = TRUE,
style = "ggplot2", # can choose between graphics, 3d, lattice or ggplot2
#study = "global",
ellipse = FALSE,
ellipse.level = 0.95,
centroid = FALSE,
star = FALSE,
legend = FALSE,
X.label = NULL,
Y.label = NULL,
Z.label = NULL,
abline = FALSE,
xlim = NULL,
ylim = NULL,
alpha = 0.2,
axes.box = "box",
plot_parameters)
{
    
    
    # --------------------------------------------------------------------------------------
    #           independent from class.object
    # --------------------------------------------------------------------------------------
    
    ### Start: Validation of arguments
    ncomp = object$ncomp
    
    #-- style
    if (!style %in% c("ggplot2", "lattice", "graphics", "3d"))
    stop("'style' must be one of 'ggplot2', 'lattice', 'graphics' or '3d' .", call. = FALSE)
    
    #-- axes.box
    choices = c("box", "bbox", "both")
    axes.box = choices[pmatch(axes.box, choices)]
    
    if (is.na(axes.box))
    stop("'axes.box' should be one of 'box', 'bbox' or 'both'.", call. = FALSE)
    
    #-- ellipse.level
    if ((ellipse.level > 1) | (ellipse.level < 0))
    stop("The value taken by 'ellipse.level' must be between 0 and 1")
    
    #-- legend
    if (length(legend) !=  1 || !is.logical(legend))
    stop("'legend' must be a logical value.", call. = FALSE)
    
    #-- alpha correlation
    if (!is.numeric(alpha) | (alpha > 1) | (alpha < 0))
    stop("The value taken by 'alpha' must be between 0 and 1", call. = FALSE)
    
    
    #-- comp
    if (is.null(comp))
    {
        if (style == "3d")
        {
            comp = c(1:3)
        } else {
            comp = c(1:2)
        }
    }
    if (length(comp) !=  2 && !(style == "3d"))
    {
        stop("'comp' must be a numeric vector of length 2.", call. = FALSE)
    } else if (length(comp) !=  3 && (style == "3d")) {
        stop("'comp' must be a numeric vector of length 3.", call. = FALSE)
    }
    
    if (!is.numeric(comp))
    stop("Invalid vector for 'comp'.")
    
    if (any(ncomp < max(comp)))
    stop(paste0("The number of components of the object to be plotted (ncomp = ", max(object$ncomp), ") is smaller than 'comp' (", paste(comp, collapse = ", "), "). Please increase ncomp  or decrease 'comp'"), call. = FALSE)
    
    comp1 = round(comp[1])
    comp2 = round(comp[2])
    if (style == "3d")
    {
        comp3 = round(comp[3])
    } else {
        comp3 = NULL
    }
    
    
    #ellipse
    if (!is.logical(ellipse))
    stop("'ellipse' must be either TRUE or FALSE", call. = FALSE)
    
    #centroid
    if (!is.logical(centroid))
    stop("'centroid' must be either TRUE or FALSE", call. = FALSE)
    
    #star
    if (!is.logical(star))
    stop("'star' must be either TRUE or FALSE", call. = FALSE)
    
    #legend
    if (!is.logical(legend))
    stop("'legend' must be either TRUE or FALSE", call. = FALSE)
    
    # abline
    if (!is.logical(abline))
    stop("'abline' must be either TRUE or FALSE", call. = FALSE)
    
    #X.label, Y.label, Z.label
    if (!is.null(X.label))
    {
        if (length(X.label)!= 1 | !is.vector(X.label))
        stop("'X.label' must be a vector of length 1", call. = FALSE)
    }
    
    if (!is.null(Y.label))
    {
        if (length(Y.label)!= 1 | !is.vector(Y.label))
        stop("'Y.label' must be a vector of length 1", call. = FALSE)
    }
    
    if (!is.null(Z.label))
    {
        if (style!= "3d")
        warning("'Z.label' is not used as style!= '3d'")
        if (length(Z.label)!= 1 | !is.vector(Z.label))
        stop("'Z.label' must be a vector of length 1", call. = FALSE)
    }
    
    # plot_parameters
    #plot_parameters = list(size.title = size.title, size.subtitle = size.subtitle, size.xlabel = size.xlabel, size.ylabel = size.ylabel, size.axis = size.axis,
    #size.legend = size.legend, size.legend.title = size.legend.title, legend.position = legend.position)
    size.title = plot_parameters$size.title
    size.subtitle = plot_parameters$size.subtitle
    size.xlabel = plot_parameters$size.xlabel
    size.ylabel = plot_parameters$size.ylabel
    size.axis = plot_parameters$size.axis
    size.legend = plot_parameters$size.legend
    size.legend.title = plot_parameters$size.legend.title
    legend.title = plot_parameters$legend.title
    legend.title.pch = plot_parameters$legend.title.pch
    legend.position = plot_parameters$legend.position
    point.lwd = plot_parameters$point.lwd
    
    if (!is.numeric(size.title) || length(size.title)>1 || size.title<0)
    stop("'size.title' needs to be a non negative number")
    
    if (!is.numeric(size.subtitle) || length(size.subtitle)>1 || size.subtitle<0)
    stop("'size.subtitle' needs to be a non negative number")
    
    if (!is.numeric(size.xlabel) || length(size.xlabel)>1 || size.xlabel<0)
    stop("'size.xlabel' needs to be a non negative number")
    
    if (!is.numeric(size.ylabel) || length(size.ylabel)>1 || size.ylabel<0)
    stop("'size.ylabel' needs to be a non negative number")
    
    if (!is.numeric(size.axis) || length(size.axis)>1 || size.axis<0)
    stop("'size.axis' needs to be a non negative number")
    
    if (!is.numeric(size.legend) || length(size.legend)>1 || size.legend<0)
    stop("'size.legend' needs to be a non negative number")
    
    if (!is.numeric(size.legend.title) || length(size.legend.title)>1 || size.legend.title<0)
    stop("'size.legend.title' needs to be a non negative number")
    
    if (length(legend.position)>1 || !legend.position%in%c("bottom", "left", "right", "top"))
    stop('"legend.position" needs to be one of "bottom", "left", "right" or "top"')
    
    if (length(legend.title)>1)
    stop("'legend.title' needs to be a single value (length 1)")

    if (length(legend.title.pch)>1)
    stop("'legend.title.pch' needs to be a single value (length 1)")

    if (!is.numeric(point.lwd) || length(point.lwd)>1 || point.lwd<0)
    stop("'point.lwd' needs to be a non negative number")
    
    if (is.logical(ind.names) & isTRUE(ind.names))
    ind.names = object$names$sample
    if (length(ind.names) > 1)
    {
        if (length(ind.names) !=  length(object$names$sample))
        stop("'ind.names' must be a character vector of length ", length(object$names$sample), " or a boolean atomic vector.")
    }
    
    
    display.names = FALSE
    if (length(ind.names) == length(object$names$sample))
    display.names = TRUE
    
    # --------------------------------------------------------------------------------------
    #           need blocks
    # --------------------------------------------------------------------------------------
    #-- xlim, ylim
    if(style %in% c("lattice", "graphics"))
    {
        if (!is.null(xlim))
        {
            if (length(blocks) == 1) # a single graph is plotted, xlim needs to be a vector of length 2
            {
                if (!is.numeric(xlim) || length(xlim) !=  2)
                stop("'xlim' must be a vector of length 2.", call. = FALSE)
                
                xlim = list(xlim)
                
            } else { # multiple graphs are plotted, xlim needs to be a list of vectors of length 2
                
                if (!is.list(xlim) || length(xlim) !=  length(blocks) || length(unlist(xlim)) !=  2 * length(blocks))
                stop("'xlim' must be a list of ", length(blocks), " vectors of length 2.", call. = FALSE)
            }
        }
        
        if (!is.null(ylim))
        {
            if (length(blocks) == 1) # a single graph is plotted, ylim needs to be a vector of length 2
            {
                if (!is.numeric(ylim) || length(ylim) !=  2)
                stop("'ylim' must be a vector of length 2.", call. = FALSE)
                
                ylim = list(ylim)
                
            } else { # multiple graphs are plotted, ylim needs to be a list of vectors of length 2
                
                if (!is.list(ylim) || length(ylim) !=  length(blocks) || length(unlist(ylim)) !=  2 * length(blocks))
                stop("'ylim' must be a list of ", length(blocks), " vectors of length 2.", call. = FALSE)
            }
        }
    }else if(style =="ggplot2") { #xlim/ylim needs to be a vector: same limits for all graphs
        if (!is.null(xlim))
        {
                if (!is.numeric(xlim) || length(xlim) !=  2)
                stop("'xlim' must be a vector of length 2.", call. = FALSE)
        }
        if (!is.null(ylim))
        {
                if (!is.numeric(ylim) || length(ylim) !=  2)
                stop("'ylim' must be a vector of length 2.", call. = FALSE)
        }
        
    }# for style = 3d, no xlim, ylim used
    
    out = list(axes.box = axes.box, comp = c(comp1, comp2, comp3), xlim = xlim, ylim = ylim, ind.names = ind.names, display.names = display.names)
    
}


internal_getVariatesAndLabels = function(object, comp, blocks.init, blocks, rep.space, style, X.label, Y.label, Z.label)
{
    
    class.object = class(object)
    object.pls = c("mixo_pls", "mixo_spls", "mixo_mlspls", "rcc")
    object.blocks = c("sgcca", "rgcca")
    object.mint = c("mint.pls", "mint.spls", "mint.plsda", "mint.splsda")
    
    #-- Start: Retrieve variates from object
    x = y = z = list()
    if (any(class.object %in%  c(object.pls, object.blocks, object.mint)))
    {
        
        
        x = lapply(object$variates, function(x){x[, comp[1], drop = FALSE]})
        y = lapply(object$variates, function(x){x[, comp[2], drop = FALSE]})
        if (style == "3d")
        {
            z = lapply(object$variates, function(x){x[, comp[3], drop = FALSE]})
        } else {
            z = NULL
        }
        
        
        #-- Variance explained for pls and block.pls obejcts
        
        # we display explained variance when only 1block is plotted and the object is not MINT
        #        if (!any(class.object%in%object.mint) & length(blocks) == 1 && rep.space !=  "XY-variate")
        if (length(blocks) == 1 && rep.space !=  "XY-variate")
        {
            if(!is.null(object$explained_variance))
            {
                if (style == "3d")
                {
                    inf = 100*round(object$explained_variance[[blocks]][c(comp[1], comp[2], comp[3])], 2)#c((object$sdev[comp[1]])^2/object$var.tot, (object$sdev[comp[2]])^2/object$var.tot, (object$sdev[comp[3]]^2)/object$var.tot)
                } else {
                    if (any(class.object%in%object.mint))
                    {
                        if (blocks%in%c("X", "Y")) # means that study == "global"
                        {
                            inf = 100*round(object$explained_variance[[blocks]]$"all data"[c(comp[1], comp[2])], 2)
                        } else {
                            inf = 100*round(object$explained_variance[[blocks.init]][[blocks]][c(comp[1], comp[2])], 2)# c((object$sdev[comp[1]])^2/object$var.tot, (object$sdev[comp[2]])^2/object$var.tot)
                            
                        }
                    } else {
                        inf = 100*round(object$explained_variance[[blocks]][c(comp[1], comp[2])], 2)
                        
                    }
                }
            } else {
                inf = rep("NC", 3)
            }
            
            # for future development: if a label with explained variance for each blocks :
            #   inf = lapply(object$explained_variance, function(x){round(x[c(comp[1], comp[2])], 2)})
            #   lapply(inf, function(x){paste0("variate ", comp[2], ": ", x[comp[2]], "% expl. var")})}
            
            
            if (is.null(X.label))
            {
                percentage = paste0(inf[1], "% expl. var")
                if (rep.space == "multi")
                X.label = paste0("variate ", comp[1], ": ", percentage)
                if (rep.space == "X-variate")
                X.label = paste0("X-variate ", comp[1], ": ", percentage)
                if (rep.space == "Y-variate")
                X.label = paste0("Y-variate ", comp[1], ": ", percentage)
                if (rep.space == "XY-variate")
                X.label = paste0("XY-variate ", comp[1], ": ", percentage)
            }
            if (is.null(Y.label))
            {
                percentage = paste0(inf[2], "% expl. var")
                if (rep.space == "multi")
                Y.label = paste0("variate ", comp[2], ": ", percentage)
                if (rep.space == "X-variate")
                Y.label = paste0("X-variate ", comp[2], ": ", percentage)
                if (rep.space == "Y-variate")
                Y.label = paste0("Y-variate ", comp[2], ": ", percentage)
                if (rep.space == "XY-variate")
                Y.label = paste0("XY-variate ", comp[2], ": ", percentage)
            }
            if (is.null(Z.label)&&style == "3d")
            {
                percentage = paste0(inf[3], "% expl. var")
                if (rep.space == "multi")
                Z.label = paste0("variate ", comp[3], ": ", percentage)
                if (rep.space == "X-variate")
                Z.label = paste0("X-variate ", comp[3], ": ", percentage)
                if (rep.space == "Y-variate")
                Z.label = paste0("Y-variate ", comp[3], ": ", percentage)
                if (rep.space == "XY-variate")
                Z.label = paste0("XY-variate ", comp[3], ": ", percentage)
            }
        }
        
        
        
        if (is.null(X.label))
        {
            if (rep.space == "multi")
            X.label = paste0("variate ", comp[1])
            if (rep.space == "X-variate")
            X.label = paste("X-variate", comp[1])
            if (rep.space == "Y-variate")
            X.label = paste("Y-variate", comp[1])
            if (rep.space == "XY-variate")
            X.label = paste("XY-variate", comp[1])
        }
        if (is.null(Y.label))
        {
            if (rep.space == "multi")
            Y.label = paste("variate", comp[2])
            if (rep.space == "X-variate")
            Y.label = paste("X-variate", comp[2])
            if (rep.space == "Y-variate")
            Y.label = paste("Y-variate", comp[2])
            if (rep.space == "XY-variate")
            Y.label = paste("XY-variate", comp[2])
        }
        if (is.null(Z.label)&&style == "3d")
        {
            if (rep.space == "multi")
            Z.label = paste("variate", comp[3])
            if (rep.space == "X-variate")
            Z.label = paste("X-variate", comp[3])
            if (rep.space == "Y-variate")
            Z.label = paste("Y-variate", comp[3])
            if (rep.space == "XY-variate")
            Z.label = paste("XY-variate", comp[3])
        }
        
    }
    #-- End: Retrieve variates from object
    out = list(x = x, y = y, z = z, X.label = X.label, Y.label = Y.label, Z.label = Z.label)
    
    
}



shape.input.plotIndiv = function(object,
n, #number of total samples
blocks = NULL, # to choose which block data to plot, when using GCCA module
x, y, z,
ind.names = TRUE,
group, # factor indicating the group membership for each sample, useful for ellipse plots. Coded as default for the -da methods, but needs to be input for the unsupervised methods (PCA, IPCA...)
col.per.group,
style = "ggplot2", # can choose between graphics, 3d, lattice or ggplot2
study = "global",
ellipse = FALSE,
ellipse.level = 0.95,
centroid = FALSE,
star = FALSE,
title = NULL,
xlim = NULL,
ylim = NULL,
col,
cex,
pch,
pch.levels,
display.names,
plot_parameters)
{
    
    class.object = class(object)
    object.mint = c("mint.pls", "mint.spls", "mint.plsda", "mint.splsda")
    
    size.title = plot_parameters$size.title
    size.subtitle = plot_parameters$size.subtitle
    size.xlabel = plot_parameters$size.xlabel
    size.ylabel = plot_parameters$size.ylabel
    size.axis = plot_parameters$size.axis
    size.legend = plot_parameters$size.legend
    size.legend.title = plot_parameters$size.legend.title
    legend.title = plot_parameters$legend.title
    legend.title.pch = plot_parameters$legend.title.pch
    legend.position = plot_parameters$legend.position
    point.lwd = plot_parameters$point.lwd

    # --------------------------------------------------------------------------------------
    #           need class.object whether it's DA
    # --------------------------------------------------------------------------------------
    
    #-- Define group
    missing.group = TRUE
    if (missing(group) & is(object, "DA"))
    {
        group = object$Y#factor(map(object$ind.mat), labels = object$names$colnames$Y)
        object$ind.mat = unmap(group) # added in v6 cause $ind.mat is the scaled (if scale = TRUE) version of ind.mat( = unmap(Y))
        missing.group = FALSE #not user defined
    } else if (!missing(group)) {
        missing.group = FALSE
        if (!is.factor(group))
        group = as.factor(group)

        object$ind.mat = unmap(group)
        
        if (length(group) !=  n)
        stop("Length of 'group' should be of length ", n, ", the sample size of your data")
    } else {
        if (star || centroid || ellipse)
        warning('star , ellipse and centroid work only if !group == NULL')
        star = centroid = ellipse = FALSE
        group = factor(rep("No group", n))
        object$ind.mat = unmap(group)
    }
    
    
    
    # --------------------------------------------------------------------------------------
    #           independent from class.object
    # --------------------------------------------------------------------------------------
    
    #at this stage, we have a 'group' - user defined or DA, or by default 1 single group
    
    # col and col.per.group
    if(!missing.group) # group is user defined or DA; we require a col.per.group input, if only a 'col' input: we use it as col.per.group
    {
        if(missing(col.per.group) & !missing(col))  # we use col as a col.per.group
        {
            if(length(col) !=  nlevels(group))
            stop("Length of 'col' should be of length ", nlevels(group), " (the number of groups).")
            col.per.group = col

        } else if (missing(col.per.group) & missing(col)) { # we create a col.per.group
            
            if (nlevels(group) < 10)
            {
                #only 10 colors in color.mixo
                col.per.group = color.mixo(1:nlevels(group))
            } else {
                #use color.jet
                col.per.group = color.jet(nlevels(group))
            }

            
        } else if (!missing(col.per.group) & !missing(col)) { # we ignore 'col'
            warning("'col' is ignored as 'group' has been set.")

        } else if (!missing(col.per.group) & missing(col)) {# all good
            
        }
        
        if(length(col.per.group) !=  nlevels(group))
        stop("Length of 'col.per.group' should be of length ", nlevels(group), " (the number of groups).")


        # make a vector of one color per sample from col.per.group
        levels.color = vector(, n)
        if (length(col.per.group) !=  n)
        {
            for (i in 1 : nlevels(group))
            levels.color[group == levels(group)[i]] = col.per.group[i]
        } else {
            levels.color = col.per.group
        }
        
    } else { #missing group, we require a 'col' of length n (or repeated to length n) and not a 'col.per.group'
        # col creates a group argument, which creates a col.per.group (levels from 'col')
        if(!missing(col.per.group))
        warning("'col.per.group' is ignored as 'group' has not been set.")
        
        if(!missing(col))
        {
            if (length(col) > n)
            stop("Length of 'col' should be of length inferior or equal to ", n, ".")

            col = rep(col, ceiling(n/length(col)))[1 : n]
            group = factor(col)
            levels.color = col
            col.per.group = levels(group)
            object$ind.mat = unmap(group)
        } else { # no group, no color => a single color by default
            col.per.group = color.mixo(1)
            levels.color = rep(col.per.group, n)
            col = levels.color
        }

    }
    # 'group' and 'col' are always the same factor, but different values
    # here we have a group, a col.per.group (length nlevels(group)) and a levels.color (length n)
    

    #-- cex argument
    if (missing(cex))
    {
        if (style == "ggplot2")
        {
            cex = rep(3, n)
            cex = cex[as.factor(group)]
        } else {
            cex = rep(1, n)
            cex = cex[as.factor(group)]
        }
    } else {
        if (length(cex) == 1)
        {
            cex = rep(cex, n)
            cex = cex[as.factor(group)]
        } else if (length(cex) == length(unique(group))) {
            cex = cex[as.factor(group)]
        } else {
            if(length(unique(group))>1)
            {
                stop("Length of 'cex' should be either 1 or ",length(unique(group)), ".")
            }else{ # one group
                stop("Length of 'cex' should be 1.")
            }
            #cex = rep(cex, ceiling(n/length(cex)))[1 : n]
        }
    }

    if (ellipse)
    {
        # removing calculation for classes with only one sample
        
        nlevels.class = 1 : ncol(object$ind.mat)
        ind.unique = which(apply(object$ind.mat, 2, sum) == 1)
        if(length(ind.unique) > 0)
        {
            nlevels.class = nlevels.class[-ind.unique]
        }
       
        #-- Start: Computation ellipse
        min.ellipse = max.ellipse = xlim.min = xlim.max = ylim.min = ylim.max = list()
        ind.gp = matrice = cdg = variance = list()
        ind.gp = lapply(1 : ncol(object$ind.mat), function(x){which(object$ind.mat[, x] == 1)})
        matrice = lapply(1 : length(x), function(z1) {lapply(ind.gp, function(z2){matrix(c(x[[z1]][z2], y[[z1]][z2]), ncol = 2)})})
        cdg = lapply(1 : length(x), function(z){ lapply(matrice[[z]], colMeans)})
        variance = lapply(1 : length(x), function(z){lapply(matrice[[z]], var)})
        coord.ellipse = lapply(1 : length(x), function(z1){ lapply(nlevels.class, function(z2){ellipse(variance[[z1]][[z2]],
            centre = cdg[[z1]][[z2]],
            level = ellipse.level)})})
        max.ellipse = lapply(1 : length(x), function(z1) {sapply(coord.ellipse[[z1]], function(z2){apply(z2, 2, max)})})
        min.ellipse = lapply(1 : length(x), function(z1) {sapply(coord.ellipse[[z1]], function(z2){apply(z2, 2, min)})})
        #-- End: Computation ellipse
        
        if (is.null(xlim))
        xlim = lapply(1 : length(x), function(z) {c(min(x[[z]], min.ellipse[[z]][1, ]), max(x[[z]], max.ellipse[[z]][1, ]))})
        if (is.null(ylim))
        ylim = lapply(1 : length(x), function(z) {c(min(y[[z]], min.ellipse[[z]][2, ]), max(y[[z]], max.ellipse[[z]][2, ]))})
        if(style == "ggplot2") # no lists, a single vector of two values is expected
        {
            temp = matrix(unlist(xlim),ncol=2,byrow=TRUE)
            xlim = c(min(temp[,1]),max(temp[,2]))
            temp = matrix(unlist(ylim),ncol=2,byrow=TRUE)
            ylim = c(min(temp[,1]),max(temp[,2]))
        }
        
    }
    # no need for xlim and ylim as ggplot2, lattice and graphics are good without by default
    #else {
    #  if (is.null(xlim))
    #    xlim = lapply(1 : length(x), function(z) {c(min(x[[z]]), max(x[[z]]))})
    #    if (is.null(ylim))
    #    ylim = lapply(1 : length(x), function(z) {c(min(y[[z]]), max(y[[z]]))})
    #}
    
    
    # --------------------------------------------------------------------------------------
    #           not independent from class.object: for the title of the plot, either "PlotIndiv" or "block:.."
    # --------------------------------------------------------------------------------------
    
    #-- pch argument
    if (missing(pch) & !any(class.object%in%object.mint))
    {
      
        if (style == "3d")
        {
            pch = unlist(lapply(1 : length(length(levels(group))), function(x){rep(c("sphere", "tetra", "cube", "octa", "icosa", "dodeca")[x], length(group == x))}))
            
        } else {
            pch = as.numeric(group)
        }
        pch.levels = pch

    }else if (any(class.object%in%object.mint)) {
        if (missing(pch))
        {
            # a pch per study, forced
            pch = as.numeric(object$study)
        } else {
            if (length(pch)!= length(object$study))
            stop("'pch' needs to be of length 'object$study' as each of 'pch' represents a specific study", call. = FALSE)
        }
        pch.levels = pch

    } else {
        if (style == "3d")
        {
            if (!all(unlist(pch) %in% c("sphere", "tetra", "cube", "octa", "icosa", "dodeca")))
            stop("pch' must be a simple character or character vector from {'sphere', 'tetra', 'cube', 'octa', 'icosa', 'dodeca'}.",
            call. = FALSE)
        }
        
        if(!missing(pch.levels))
        {
            if(length(pch.levels) != length(pch))
            stop("'pch.levels' needs to be a vector of the same length as 'pch': ", length(pch))
        } else {
            pch.levels = pch
        }
        
        if (length(pch) == 1)
        {
            pch = rep(pch, n)
            pch.levels = rep(pch.levels, n)

        } else if (length(pch) > n) {
            stop("Length of 'pch' should be of length inferior or equal to ", n, ".")
        } else if (length(pch) == length(unique(group)) & length(pch)!=n ) { # prevent from reordering pch when 1 group per sample (length(pch)=length(group)=n)
            pch = pch[as.factor(group)]
            pch.levels = pch.levels[as.factor(group)]

        } else {
            pch = rep(pch, ceiling(n/length(pch)))[1 : n]
            pch.levels = rep(pch.levels, ceiling(n/length(pch.levels)))[1 : n]
        }
        
        # if pch is given and ind.names is TRUE, pch takes over
        if(display.names)
        warning("'ind.names' is set to FALSE as 'pch' overrides it")

        display.names = FALSE
      
    }
    
    
    
    # constructing data.frame df
    if (any(study == "global"))# | length(study) == 1)
    {
        #-- Start: data set
        df = list()
        if (style == "3d")
        {
            for (i in 1 : length(x))
            df[[i]] = data.frame(x = x[[i]], y = y[[i]], z = z[[i]], group = group)
        } else {
            for (i in 1 : length(x))
            df[[i]] = data.frame(x = x[[i]], y = y[[i]], group = group)
        }
        
        title.save = title # to use for ellipse
        #if (any(class.object %in% c("ipca", "sipca", "pca", "spca", "prcomp", "splsda", "plsda")) &
        if(length(blocks) == 1 & !any(class(object)%in%c(object.mint, "sgcca", "rgcca"))) # add blocks == 1 to allow "multi" with plsda
        {
            if (is.null(title))
            {
                df = data.frame(do.call(rbind, df), "Block" = "PlotIndiv")
                if (style %in%c("graphics"))
                title = "PlotIndiv" # to add title to graphics
                
            } else {
                df = data.frame(do.call(rbind, df), "Block" = title)
                if (style %in%c("ggplot2", "lattice"))
                title = NULL # to avoid double title
                
            }
            
            # no subtitle with these objects
            if(size.title != rel(2)) # rel(2) is the default
            size.subtitle = size.title


            df$Block = as.factor(df$Block)
        } else {
            df = data.frame(do.call(rbind, df), "Block" = paste0("Block: ", unlist(lapply(1 : length(df), function(z){rep(blocks[z], nrow(df[[z]]))}))))
            df$Block = factor(df$Block, levels = unique(df$Block))
        }
        
        if (style == "3d")
        {
            names(df)[1:3] = c("x", "y", "z")
        } else {
            names(df)[1:2] = c("x", "y")
        }
        
        if (display.names)
        df$names = rep(ind.names, length(x))
        
        df$pch = pch; df$pch.levels = pch.levels
        df$cex = cex
        #df$col.per.group = levels.color#[group] #FR: don't understand what is that changing as levels.color is already group?
        df$col = levels.color#as.character(col)
        
        
        if (centroid == TRUE || star == TRUE)
        {
            df = cbind(df, rep(0, nrow(df)))
            n = ncol(df)
            df = cbind(df, rep(0, nrow(df)))
            for (i in 1:nlevels(group))
            {
                if (length(x)>1)
                {
                    for (k in 1 : length(x))
                    {
                        x0 = mean(df[df$group == levels(group)[i] & df$Block %in% paste0("Block: ", blocks[k]), "x"])
                        y0 = mean(df[df$group == levels(group)[i] & df$Block %in% paste0("Block: ", blocks[k]) , "y"])
                        df[df$group == levels(group)[i] & df$Block %in% paste0("Block: ", blocks[k]), n] = x0
                        df[df$group == levels(group)[i] & df$Block %in% paste0("Block: ", blocks[k]), n+1] = y0
                        names(df)[c(ncol(df)-1, ncol(df))] = c("x0", "y0")
                    }
                } else {
                    x0 = mean(df[df$group == levels(group)[i] , "x"])
                    y0 = mean(df[df$group == levels(group)[i]  , "y"])
                    df[df$group == levels(group)[i] , n] = x0
                    df[df$group == levels(group)[i] , n+1] = y0
                    names(df)[c(ncol(df)-1, ncol(df))] = c("x0", "y0")
                }
                
            }
        }
        
        if (ellipse == TRUE)
        {
            df.ellipse = data.frame(do.call("rbind", lapply(1 : length(x), function(k){do.call("cbind", coord.ellipse[[k]])})), "Block" = paste0("Block: ", rep(blocks, each = 100)))
            if(length(ind.unique) > 0)
            {
                names(df.ellipse)[1 : (2*length(nlevels.class))] = paste0("Col", c(1 : (2*nlevels(group)))[-c(2*(ind.unique-1)+1,2*ind.unique)])
            } else {
                names(df.ellipse)[1 : (2*nlevels(group))] = paste0("Col", c(1 : (2*nlevels(group))))

            }
            df.ellipse$ellipse.level = ellipse.level
        } else {
            df.ellipse = NULL
        }
        
        if(ellipse == TRUE && length(blocks) == 1 & !any(class(object)%in%c(object.mint, "sgcca", "rgcca"))) # add blocks == 1 to allow "multi" with plsda
        #if (ellipse == TRUE && any(class.object %in% c("ipca", "sipca", "pca", "spca", "prcomp", "splsda", "plsda"))& length(blocks) == 1& !any(class(object)%in%object.mint))
        {
            if (is.null(title.save))
            {
                df.ellipse$Block = "PlotIndiv"
            } else {
                df.ellipse$Block = title.save
            }
        }
        
        
        
        pch.legend = NULL
        for (i in 1:nlevels(group))
        {
            pch.legend = c(pch.legend, df[df$group == levels(group)[i], ]$pch)
        }
        df$pch.legend = pch.legend
        
    } else {
        
        #mint object
        #display.names = FALSE # so far ggplot and lattice require a unique vector of names. when the code changes, we can use ind.names (list)
        group.mint = split(group, object$study)[study]
        group = as.factor(unlist(group.mint))
        pch = as.vector(unlist(split(pch, object$study)[study]))
        cex = as.vector(unlist(split(cex, object$study)[study]))
        pch.levels = as.vector(unlist(split(pch.levels, object$study)[study]))
        col.per.group.mint = as.vector(unlist(split(levels.color, object$study)[study]))
        #col = as.vector(unlist(split(col, object$study)[study]))
        
        
        
        #-- Start: data set
        df = list()
        if (style == "3d")
        {
            for (i in 1 : length(x))
            {
                df[[i]] = data.frame(x = x[[i]], y = y[[i]], z = z[[i]], group = group.mint[[i]])
            }
        } else {
            for (i in 1 : length(x))
            {
                df[[i]] = data.frame(x = x[[i]], y = y[[i]], group = group.mint[[i]])
            }
        }
        
        df = data.frame(do.call(rbind, df), "Block" = paste0("Study: ", unlist(lapply(1 : length(df), function(z){rep(blocks[z], nrow(df[[z]]))}))))
        df$Block = factor(df$Block, levels = unique(df$Block))
        
        #print(df)
        
        if (style == "3d")
        {
            names(df)[1:3] = c("x", "y", "z")
        } else {
            names(df)[1:2] = c("x", "y")
        }
        
        # no names for MINT object
        #if (display.names)
        #df$names = rep(ind.names, length(x))
        
        df$pch = pch; df$pch.levels = pch.levels
        df$cex = cex
        #df$col.per.group = col.per.group.mint;
        df$col = col.per.group.mint#as.character(col)
        
        
        pch.legend = NULL
        for (i in 1:nlevels(group))
        pch.legend = c(pch.legend, df[df$group == levels(group)[i], ]$pch)
        
        df$pch.legend = pch.legend
        df.ellipse = NULL # no ellipse so far
    }
    
    if (any(study == "global"))
    study = levels(object$study)
    
    # match study with names of the study
    study.ind = match(study, levels(object$study))
    
    #print(df)
    plot_parameters = list(size.title = size.title, size.subtitle = size.subtitle, size.xlabel = size.xlabel, size.ylabel = size.ylabel,
    size.axis = size.axis, size.legend = size.legend, size.legend.title = size.legend.title, legend.title = legend.title,
    legend.title.pch = legend.title.pch, legend.position = legend.position, point.lwd = point.lwd)
    
    out = list(df = df, study.ind = study.ind, df.ellipse = df.ellipse, col.per.group = col.per.group, title = title, display.names = display.names, xlim = xlim, ylim = ylim, ellipse = ellipse, centroid = centroid, star = star, plot_parameters = plot_parameters)
}


# --------------------------------------------------------------------------------------
# Internal helpers functions to run some plots functions
# --------------------------------------------------------------------------------------

#-- Function to display an error message (used for the parameters var.names, cex, col, pch and font)
stop.message = function(argument, data){
    if (length(data) == 1) {
        count.data = sapply(data, length)
    } else {
        count.data = paste(paste(sapply(data[-length(data)], length), collapse =  ", "), length(data[[length(data)]]), sep = " and ")
    }
    stop(argument, " must be either a vector of length ", length(data),
    " or a list of ", length(data), " vector components of length ", count.data, " respectively.",call.= FALSE)
}
