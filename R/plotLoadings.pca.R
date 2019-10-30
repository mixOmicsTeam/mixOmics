#############################################################################################################
# Author :
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Kim-Anh Le Cao, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 15-04-2016
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
#-- Includes plotLoadings for PLS, sPLS, rCC, PCA, sPCA, IPCA, sIPCA, rGCCA, sGCCA --#
#----------------------------------------------------------------------------------------------------------#


plotLoadings.pca =

function(object,
comp = 1,
col = NULL,
ndisplay = NULL,
size.name = 0.7,
name.var = NULL,
name.var.complete = FALSE, #name.var.complete
title = NULL,
size.title = rel(2),
layout = NULL,
border = NA,
xlim = NULL,
...)
{
    
    # -- input checks
    object$names$blocks = "X"
    check = check.input.plotLoadings(object = object, block = "X", size.name = size.name, title = title, col = col, name.var = name.var, xlim = xlim)
    
    col = check$col
    size.name = check$size.name
    block = 1
    object$X = list(X=object$X)
    xlim = check$xlim


    # -- layout
    res = layout.plotLoadings(layout = layout, plot = TRUE, legend = FALSE, block = block)
    reset.mfrow = res$reset.mfrow
    opar = res$opar
    omar = par("mar") #reset mar at the end

    res = get.loadings.ndisplay(object = object, comp = comp, block = block, name.var = name.var, name.var.complete = name.var.complete, ndisplay = ndisplay)
    X = res$X
    names.block = res$names.block
    colnames.X = res$colnames.X
    value.selected.var = res$value.selected.var

    df = data.frame(importance = value.selected.var) # contribution of the loading
    
    # barplot with contributions
    colnames.X <- .trim_long_names(colnames.X) ## issue 45
    par(mar = c(4, max(7, max(sapply(colnames.X, nchar),na.rm = TRUE)/3), 4, 2))

    barplot(df$importance, horiz = TRUE, las = 1, col = col, axisnames = TRUE, names.arg = colnames.X, #names.arg = row.names(df),
    cex.names = size.name, cex.axis = 0.7, beside = TRUE, border = border, xlim = xlim)
    
    if (is.null(title))
    {
        title(paste0('Loadings on comp ', comp), cex.main = size.title)
    } else {
        title(paste(title), cex.main = size.title)
    }
    
    if (reset.mfrow)
    par(opar)#par(mfrow = c(1,1))

    par(mar = omar) #reset mar

    # return the contribution matrix
    return(invisible(df))
}
