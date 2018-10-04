#############################################################################################################
# Authors:
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 16-03-2016
# last modified: 25-08-2016
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
#-- Includes plotIndiv for the MINT module --#
#----------------------------------------------------------------------------------------------------------#

plotIndiv.mint.pls      =
plotIndiv.mint.spls     =
plotIndiv.mint.plsda    =
plotIndiv.mint.splsda   =

function(object, 
comp = NULL, 
study = "global",
rep.space = NULL,
group, # factor indicating the group membership for each sample, useful for ellipse plots. Coded as default for the -da methods, but needs to be input for the unsupervised methods (PCA, IPCA...)
col.per.group, 
style = "ggplot2", # can choose between graphics, lattice or ggplot2
ellipse = FALSE,
ellipse.level = 0.95, 
centroid = FALSE,
star = FALSE,
title = NULL,
subtitle, 
legend = FALSE,
X.label = NULL, 
Y.label = NULL, 
abline = FALSE, 
xlim = NULL, 
ylim = NULL, 
col, 
cex, 
pch, 
layout = NULL, 
size.title = rel(2), 
size.subtitle = rel(1.5), 
size.xlabel = rel(1), 
size.ylabel = rel(1), 
size.axis = rel(0.8), 
size.legend = rel(1), 
size.legend.title = rel(1.1), 
legend.title = "Legend",
legend.position = "right",
point.lwd = 1, 
...
)
{
    plot_parameters = list(size.title = size.title, size.subtitle = size.subtitle, size.xlabel = size.xlabel, size.ylabel = size.ylabel, size.axis = size.axis, 
    size.legend = size.legend, size.legend.title = size.legend.title, legend.title = legend.title,
    legend.position = legend.position, point.lwd = point.lwd)


    if (any(class(object)%in%c("mint.block.pls", "mint.block.spls", "mint.block.plsda", "mint.block.splsda")))
    stop("No plotIndiv for the following functions at this stage: mint.block.pls, mint.block.spls, mint.block.plsda, mint.block.splsda.")
    

    #-- rep.space
    if (is.null(rep.space))#"splsda", "plsda", "mlsplsda")))
    rep.space = "X-variate"

    rep.space  =  match.arg(rep.space, c("XY-variate", "X-variate", "Y-variate", "multi"))

    ind.names = FALSE
    # --------------------------------------------------------------------------------------
    #           need study
    # --------------------------------------------------------------------------------------
    
    # check study
    #study needs to be either: from levels(object$study), numbers from 1:nlevels(study) or "global"
    if (any(!study%in%c(levels(object$study), "global" , "all.partial")))
    stop("'study' must be one of 'object$study', 'global' or 'all.partial', see help file.")
    
    if (length(study)!=length(unique(study)))
    stop("Duplicate in 'study' not allowed")
    
    if (any(study != "global"))
    {
        if (ellipse == TRUE)
        stop("'ellipse' must be FALSE when study is different from 'global'")

        if (star == TRUE)
        stop("'star' must be FALSE when study is different from 'global'")
    }


    #LOOP ON STUDY, to get a plot with every single one, could be a mixed of numbers and "global", only if there is both "global" and something else.
    
    object.init = object
    study.init = unique(study)
    
    # replace "all.partial" by all levels of object$study
    ind.all.partial = which(study.init == "all.partial")
    if (length(ind.all.partial) > 0)
    {
        if (ind.all.partial > 1 & ind.all.partial < length(study.init))
        {
            # there are things before and after "all.partial"
            study.init = c(study.init[1:(ind.all.partial-1)], levels(object$study), study.init[(ind.all.partial+1) : length(study.init)])
        } else if (ind.all.partial == 1 & ind.all.partial < length(study.init)) {
            # there are only things after "all.partial"
            study.init = c(levels(object$study), study.init[(ind.all.partial+1) : length(study.init)])
        } else if (ind.all.partial > 1 & ind.all.partial == length(study.init)) {
            # there are things only before "all.partial"
            study.init = c(study.init[1:(ind.all.partial-1)], levels(object$study))
        } else if (ind.all.partial == 1 & ind.all.partial == length(study.init)) {
            # there's only "all.partial"
            study.init = levels(object$study)
        }
    }
    
    study.init = unique(study.init) #once again cause we added studies if "all.partial"
    
    if (!missing(subtitle))
    {
        if (length(subtitle)!=length(study.init)| length(subtitle)!=length(unique(subtitle)))
        stop("'subtitle' indicates the subtitle of the plot for each study and it needs to be the same length as 'study' (", length(study.init),") and duplicate are not allowed. 'study' includes: ", paste(study.init, collapse = ", "))
    }
    
    df.final = data.frame()
    
    indice.all = grep("global", study.init) # can go faster before and after "global"
    if (length(indice.all)>0)
    {
        study.list = list()
        i = 1
        if (indice.all>1)
        {
            study.list[[1]] = study.init[1:(indice.all-1)]
            i = i+1
        }
        
        study.list[[i]] = study.init[indice.all]
        
        if (indice.all<length(study.init))
        study.list[[i+1]] = study.init[-(1:indice.all)]
    } else {
        study.list = list(study.init)
    }
    

    # the following loop consider subset of studies all together, up until "global", and subset of studies after "global"
    for (length.study in 1 : length(study.list))
    {
        object = object.init #reinitialise $variates
        study = study.list[[length.study]]
        
        #-- define 'blocks'
        if (any(study == "global"))
        {
            # can plot both X and Y when one study or when study="global"
            # same as class.object==pls
            
            if (rep.space == "multi")
            {
                blocks = c("X", "Y")
                object$variates = object$variates[names(object$variates) %in% blocks]
            }
            
            if (rep.space == "X-variate")
            {
                object$variates = object$variates["X"]
                blocks = "X"
            }
            
            if (rep.space == "Y-variate")
            {
                object$variates = object$variates["Y"]
                blocks = "Y"
            }
            
            if (rep.space == "XY-variate")
            {
                object$variates$XYvariates = (object$variates$X + object$variates$Y)/2
                object$variates = object$variates["XYvariates"]
                blocks = "XY combined"
            }
            
        } else if (length(study) == 1) {
            # can plot only X, Y or XY variate when more than one study
            # can plot both X and Y when one study or when study="global"
            
            blocks = c("X", "Y")
            
            if (rep.space == "X-variate")
            blocks = "X"
            
            if (rep.space == "Y-variate")
            blocks = "Y"
            
            #extract variates for each "blocks" for "study"
            object$variates = lapply(object$variates.partial, function(x){x[[study]]})[names(object$variates) %in% blocks]
            
            #if XY-variate, combine the previous variates (relative to "blocks" and "study")
            if (rep.space == "XY-variate")
            {
                object$variates$XYvariates = (object$variates$X + object$variates$Y)/2
                object$variates = object$variates["XYvariates"]
                blocks = "XY combined"
            }
            blocks.init = blocks #save for "internal_getVariatesAndLabels"
            blocks = study
            
            
        } else { #length(study)>1
            
            blocks = c("X", "Y")

            if (rep.space == "multi")
            {
                rep.space = "X-variate"
                warning("More than one study is plotted, 'rep.space' is set to 'X-variate'. Alternatively, you can input 'Y-variate'")
            }
            
            if (rep.space == "X-variate")
            blocks = "X"
            
            if (rep.space == "Y-variate")
            blocks = "Y"
            

            #extract variates for each "blocks" for "study"
            object$variates = lapply(object$variates.partial, function(x)
            {
                out = lapply(study, function(y){x[[y]]})
                names(out) = study
                out
            })[names(object$variates) %in% blocks]#[[1]]
            
            #if XY-variate, combine the previous variates (relative to "blocks" and "study")
            if (rep.space == "XY-variate")
            {
                for (i in 1:length(object$variates$X))
                object$variates$XYvariates[[i]] = (object$variates$X[[i]]+object$variates$Y[[i]])/2

                names(object$variates$XYvariates) = names(object$variates$X)
                object$variates = object$variates[["XYvariates"]]
            } else {
                object$variates = object$variates[[1]] # get rid of the $X or $Y
            }
            
            # blocks becomes study, so each study is plotted
            blocks = study
            object$names$sample = lapply(object$variates, rownames)
            ellipse = FALSE
            star = FALSE
            centroid = FALSE
        }

        #-- check inputs
        # check style as we do not do 3d at the moment:
        if (!style %in% c("ggplot2", "lattice", "graphics"))
        stop("'style' must be one of 'ggplot2', 'lattice' or 'graphics'.", call. = FALSE)
        
        check = check.input.plotIndiv(object = object, comp = comp, blocks = blocks, ind.names = ind.names, 
        style = style, ellipse = ellipse, ellipse.level = ellipse.level, centroid = centroid, 
        star = star, legend = legend, X.label = X.label, Y.label = Y.label, abline = abline, 
        xlim = xlim, ylim = ylim, plot_parameters = plot_parameters)
        #-- retrieve some outputs from the checks
        comp = check$comp
        xlim = check$xlim
        ylim = check$ylim
        ind.names = check$ind.names
        display.names = FALSE#check$display.names
        

        #-- get the variates
        variate = internal_getVariatesAndLabels(object, comp, blocks.init = blocks.init, blocks = blocks, rep.space = rep.space,
        style = style, X.label = X.label, Y.label = Y.label, Z.label = NULL)
        #-- retrieve outputs
        x = variate$x
        y = variate$y
        z = variate$z
        X.label = variate$X.label #only the last one of the loop is used
        Y.label = variate$Y.label #only the last one of the loop is used

        n = nrow(object$X)

        # create data frame df that contains (almost) all the ploting information
        out = shape.input.plotIndiv(object = object, n = n, blocks = blocks, x = x, y = y, z = z, ind.names = ind.names, group = group,
        col.per.group = col.per.group, style = style, study = study, ellipse = ellipse, ellipse.level = ellipse.level,
        centroid = centroid, star = star, title = title, xlim = xlim, ylim = ylim, 
        col = col, cex = cex, pch = pch, display.names = display.names, plot_parameters = plot_parameters)
        #-- retrieve outputs
        df = out$df
        df.ellipse = out$df.ellipse
        col.per.group = out$col.per.group
        title = out$title
        display.names = out$display.names
        xlim = out$xlim
        ylim = out$ylim
        #missing.col = out$missing.col
        plot_parameters = out$plot_parameters

        #save(list=ls(),file="temp.Rdata")
        # concatenate results
        df.final = rbind(df.final, df)
    }
    # add study information on df.final, for pch legend
    study.levels = study.init[which(!study.init == "global")]
    if (any(study.init == "global"))
    study.levels = levels(object$study)
    
    
    # change the levels of df.final$Block to "subtitle"
    if (!missing(subtitle))
    {
        df.final$Block = factor(df.final$Block, labels = subtitle)

        if(ellipse)
        df.ellipse$Block = factor(df.ellipse$Block, labels = subtitle)
    }
    df = df.final
        
    if (style == "ggplot2")
    style = "ggplot2-MINT"
        
    #call plot module (ggplot2, lattice, graphics, 3d)
    res = internal_graphicModule(df = df, centroid = centroid, col.per.group = col.per.group, title = title,
    X.label = X.label, Y.label = Y.label, xlim = xlim, ylim = ylim, class.object = class(object),
    display.names = display.names, legend = legend, abline = abline,
    star = star, ellipse = ellipse, df.ellipse = df.ellipse, style = style, layout = layout,
    #missing.col = missing.col,
    #for ggplot2-MINT
    study.levels = study.levels, plot_parameters = plot_parameters
    )

    return(invisible(list(df = df, graph = res)))
}
