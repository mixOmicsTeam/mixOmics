#############################################################################################################
# Author :
#   Kim-Anh Le Cao, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Benoit Gautier, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 19-04-2016
# last modified: 24-05-2016
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
#-- Includes plotIndiv for PLS-DA, SPLS-DA, sGCCDA --#
#----------------------------------------------------------------------------------------------------------#


plotLoadings.mixo_plsda =
#plotLoadings.mlplsda =    # because plsda too
plotLoadings.mixo_splsda =
#plotLoadings.mlsplsda =   # because splsda too
plotLoadings.sgccda =


function(object,
contrib = NULL,  # choose between 'max" or "min", NULL does not color the barplot
method = "mean", # choose between 'mean" or "median"
block, #single value, for sgccda object
comp = 1,
plot = TRUE,
show.ties = TRUE,
col.ties = "white",
ndisplay = NULL,
size.name = 0.7,
size.legend = 0.8,
name.var = NULL,
name.var.complete = FALSE,
title = NULL,
subtitle,
size.title = rel(1.8),
size.subtitle = rel(1.4),
legend = TRUE,
legend.color = NULL,
legend.title = 'Outcome',
layout = NULL,
border = NA,
xlim = NULL,
...
) {
    
    # -- input checks
    check = check.input.plotLoadings(object = object, block = block, subtitle = subtitle, size.name = size.name, size.legend = size.legend,
    title = title, col = NULL, contrib = contrib, name.var = name.var, xlim = xlim)
    
    size.name = check$size.name
    size.legend = check$size.legend
    block = check$block
    xlim = check$xlim

    # contrib
    # --
    
    # if contrib is NULL, then we switch to the classical plotLoadings (without contribution/colors)
    if(is.null(contrib))
    {
        if(plot)
        {
            plotLoadings.mixo_pls(object = object, block = block, comp = comp, ndisplay = ndisplay,
            size.name = size.name,
            name.var = name.var,
            name.var.complete = name.var.complete,
            title = title,
            subtitle = subtitle,
            xlim = xlim,
            layout = layout,
            size.title = size.title,
            size.subtitle = size.subtitle,
            border = TRUE,
            col = "white")
        } else {
            stop("'contrib' is NULL and 'plot' is FALSE => no results to show", call. = FALSE)
        }
        # stop the script without error message
        # blankMsg <- sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" "))
        # stop(simpleError(blankMsg))
    } else {
        
        # -- layout
        res = layout.plotLoadings(layout = layout, plot = plot, legend = legend, block = block)
        reset.mfrow = res$reset.mfrow
        opar = res$opar
        omar = par("mar") #reset mar at the end
        
        # method
        # ----
        if (length(method) !=1 || !method %in% c("mean","median"))
        {
            method = "median"
            warning("'method' should be either 'mean' or 'median', set to 'median' by default")
        }
        
        if (length(block) == 1 & !is.null(name.var))
        name.var = list(name.var = name.var)
        
        for (i in 1 : length(block))
        {
            res = get.loadings.ndisplay(object = object, comp = comp, block = block[i], name.var = name.var[[i]], name.var.complete = name.var.complete, ndisplay = ndisplay)
            X = res$X
            names.block = res$names.block
            colnames.X = res$colnames.X
            name.selected.var = res$name.selected.var
            value.selected.var = res$value.selected.var
            
            Y = object$Y #v6: all $Y are factors for DA methods
            
            #legend.color
            #-----
            if (!is.null(legend.color) & (length(legend.color) != nlevels(Y)))
            {
                warning('legend.color must be the same length than the number of group, by default set to default colors')
                legend.color = color.mixo(1:10)  # by default set to the colors in color.mixo (10 colors)
            }
            if (is.null(legend.color))
            legend.color = color.mixo(1:10)[1:nlevels(Y)] # by default set to the colors in color.mixo (10 colors)
            
            if (col.ties%in%legend.color[1:nlevels(Y)])
            stop("'col.ties' should not be in 'legend.color'")
            
            #  determine the colors/groups matching max contribution
            df = get.contrib.df(Y = Y, X = X, method = method, contrib = contrib, value.selected.var = value.selected.var, colnames.X = colnames.X, name.selected.var = name.selected.var, legend.color = legend.color, col.ties = col.ties)
            
            # when working with sparse counts in particular and using the median to measure contribution
            # ties to determine the contribution of a variable may happen, in that case remove them, otherwise they are showns as blank
            if (show.ties == FALSE)
            {
                df = df[!df$color %in% col.ties, ]
                colnames.X = rownames(df)
            }

            # display barplot with names of variables
            colnames.X <- .trim_long_names(colnames.X) ## issue 45
            if (plot) # condition if all we need is the contribution stats
            {
                if (!is.null(title) & length(block) > 1)
                {
                    par(mar = c(4, max(7, max(sapply(colnames.X, nchar),na.rm = TRUE)/3), 6, 2))
                } else {
                    par(mar = c(4, max(7, max(sapply(colnames.X, nchar), na.rm = TRUE)/3), 4, 2))
                }
                
                barplot(df$importance, horiz = TRUE, las = 1, col = df$color, axisnames = TRUE, names.arg = colnames.X, #names.arg = row.names(df),
                cex.names = size.name, cex.axis = 0.7, beside = TRUE, border = border, xlim = xlim[i, ])
                
                if ( length(block) == 1 & is.null(title) )
                {
                    title(paste0('Contribution on comp ', comp), line=0, cex.main = size.title)
                } else if (length(block) == 1) {
                    title(paste(title), line=1, cex.main = size.title)
                } else if ((length(block) > 1 & missing(subtitle))) {
                    title(paste0('Contribution on comp ', comp, "\nBlock '", names.block,"'"), line=0, cex.main = size.subtitle)
                } else if (length(block) > 1 & !missing(subtitle)) {
                    title(paste(subtitle[i]), line=1, cex.main = size.subtitle)
                }
                
                if (legend)
                {
                    par(mar = c(5, 0, 4, 3) + 0.1)
                    plot(1,1, type = "n", axes = FALSE, ann = FALSE)
                    legend(0.8, 1.1, col = legend.color[1:nlevels(Y)], legend = levels(Y), pch = 19,
                    title = paste(legend.title),
                    cex = size.legend)
                }
            } # end if plot
        }
        
        if(plot) # overall title and reset par if needed
        {
            # legend
            if (length(block) > 1 & !is.null(title))
            title(title, outer=TRUE, line = -2, cex.main = size.title)
            
            if (reset.mfrow)
            par(opar)#par(mfrow = omfrow)

            par(mar = omar) #reset mar
        }
        
        # return the contribution matrix
        return(invisible(df))
    }# end contrib missing
}
