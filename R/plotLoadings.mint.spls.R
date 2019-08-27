#############################################################################################################
# Author :
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Kim-Anh Le Cao, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
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
#-- Includes plotLoadings for mint.pls and mint.spls --#
#----------------------------------------------------------------------------------------------------------#


plotLoadings.mint.pls    =
plotLoadings.mint.spls   =

function(object, 
study = "global",
comp = 1,
col = NULL,
ndisplay = NULL,
size.name = 0.7,
name.var = NULL,
name.var.complete = FALSE,
title = NULL,
subtitle,
size.title = rel(1.8),
size.subtitle = rel(1.4),
layout = NULL,
border = NA,
xlim = NULL,
...
) {
    
    # what I want is to modify the input and call plotLoadings.pls and plotLoadings.splsda where blocks are now studies
    # do not forget to change object$names$block in levels(object$study) and it should work, see you tomorrow
    
    if(any(study == "global"))
    {
        # if study == "global" then we plot the results on the concatenated data, thus direct call to plotLoadings.plsda
        plotLoadings.mixo_pls(object = object, block = "X", comp = comp, ndisplay = ndisplay,
        size.name = size.name,
        name.var = name.var,
        name.var.complete = name.var.complete,
        title = title,
        subtitle = subtitle,
        layout = layout,
        size.title = size.title,
        size.subtitle = size.subtitle,
        border = border,
        xlim = xlim,
        col = col)
        
    } else {
        # if study != "global" then we plot the results on each study

        # -- input checks
        check = check.input.plotLoadings(object = object, block = "X", title = title, col = col, size.name = size.name, name.var = name.var)

        col = check$col
        size.name = check$size.name
        block = check$block # "X"

        #study needs to be either: from levels(object$study), numbers from 1:nlevels(study) or "global"
        if (any(!study%in%c(levels(object$study), "global" , "all.partial")))
        stop("'study' must from one of 'object$study', 'global' or 'all.partial', see help file.")

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
        study = study.init

        if (!isNULL(subtitle))
        {
            if (length(subtitle)!=length(study))
            stop("'subtitle' indicates the subtitle of the plot for each study and it needs to be the same length as 'study' (", length(study),"), which includes: ", paste(study, collapse = ", "))
        }


        # swap block for study
        block = study
        
        # check xlim, has to be a matrix with number of rows=number of studies, or a vector of two values
        if(length(study) == 1 & !is.null(xlim))
        {
            if(length(xlim) !=2)
            stop("'xlim' must be a vector of length 2")
            
            xlim = matrix(xlim, nrow = 1)
        }
        
        
        if(length(study)>1 & !is.null(xlim))
        {
            if(is.matrix(xlim) && ( !nrow(xlim) %in%c(1, length(study))  | ncol(xlim) != 2 ))
            stop("'xlim' must be a matrix with ",length(study)," rows (length(study)) and 2 columns")
            
            if(is.vector(xlim))
            {
                if(length(xlim) !=2)
                stop("'xlim' must be a matrix with ",length(study)," rows (length(study)) and 2 columns")
                
                xlim = matrix(xlim, nrow = 1)
            }
            
            if(nrow(xlim) != length(study)) # we complete xlim to have one xlim per block
            xlim = matrix(rep(xlim, length(study)), nrow = length(study), byrow=TRUE)
        }
        
        
        # -- layout
        res = layout.plotLoadings(layout = layout, plot = TRUE, legend = FALSE, block = block)
        reset.mfrow = res$reset.mfrow
        opar = res$opar
        omar = par("mar") #reset mar at the end
        
        # get the selected variables on the concatenated data
        res = get.loadings.ndisplay(object = object, comp = comp, block = "X", name.var = name.var, name.var.complete = name.var.complete, ndisplay = ndisplay)
        X = res$X
        colnames.X = res$colnames.X
        name.selected.var = res$name.selected.var
        value.selected.var = res$value.selected.var


        # swap loadings partial for loadings
        object$loadings.global = object$loadings
        object$loadings = object$loadings.partial$X
        object$names$block = levels(object$study)
        
        df.final = list()
        for (i in 1 : length(block))
        {
            value.selected.var = object$loadings.partial$X [[block[i]]][, comp] [name.selected.var]

            df = data.frame(importance = value.selected.var, color = col, stringsAsFactors = FALSE) # contribution of the loading
           
            #display barplot with names of variables
            #added condition if all we need is the contribution stats
            if (!is.null(title) & length(block) > 1)
            {
                par(mar = c(4, max(7, max(sapply(colnames.X, nchar),na.rm = TRUE)/2), 6, 2))
            } else {
                par(mar = c(4, max(7, max(sapply(colnames.X, nchar),na.rm = TRUE)/2), 4, 2))
            }

            mp = barplot(df$importance, horiz = TRUE, las = 1, col = df$color, axisnames = TRUE, names.arg = colnames.X, #names.arg = row.names(df),
            cex.names = size.name, cex.axis = 0.7, beside = TRUE, border = border, xlim = xlim)
            
            if ( length(block) == 1 & is.null(title) )
            {
                title(paste0('Loadings on comp ', comp), line=1, cex.main = size.title)
            } else if (length(block) == 1) {
                title(paste(title), line=0, cex.main = size.title)
            } else if ((length(block) > 1 & missing(subtitle))) {
                title(paste0('Loadings on comp ', comp, "\nStudy '", block[i],"'"), line=0, cex.main = size.subtitle)
            } else if (length(block) > 1 & !missing(subtitle)) {
                title(paste(subtitle[i]), line=0, cex.main = size.subtitle)
            }
            
            df.final[[i]] = df
        }
        names(df.final) = block
        
        if (length(block) > 1 & !is.null(title))
        title(title, outer=TRUE, line = -2, cex.main = size.title)
        
        if (reset.mfrow)
        par(opar)#par(mfrow = omfrow)
        
        par(mar = omar) #reset mar
        
        return(invisible(df.final))

    }
}
