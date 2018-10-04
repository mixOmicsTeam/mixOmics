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


plotLoadings  =
function(object, ...) UseMethod("plotLoadings")


# --------------------------------------------------------------------------------------
# Internal helpers functions to run "plotLoadings" functions
# --------------------------------------------------------------------------------------


check.input.plotLoadings = function(object, block, study, subtitle, size.name, size.legend, title, col, contrib, name.var, xlim)
{
    
    if (is.null(object$loadings))
    stop("'plotLoadings' should be used on object for which object$loadings is present.")
    
    # block
    # --
    if (missing(block))
    {
        if (all(class(object) != "DA"))
        {
            block = object$names$blocks
        } else  if (any(class(object) %in% c("mixo_plsda", "mixo_splsda"))) {
            block = "X"
        } else {
            if (!is.null(object$indY))
            {
                block = object$names$blocks[-object$indY]
            } else {
                block = object$names$blocks
            }
        }
    }
    
    if (any(class(object) %in% c("mixo_plsda", "mixo_splsda")) & (!all(block %in% c(1,"X")) | length(block) > 1 ))
    stop("'block' can only be 'X' or '1' for plsda and splsda object")
    
    if (any(class(object) %in% c("mixo_plsda", "mixo_splsda","pca")))
    {
        object$indY = 2
    } else if (any(class(object) %in% c("mixo_pls", "mixo_spls"))) {
        object$indY = 3 # we don't want to remove anything in that case, and 3 is higher than the number of blocks which is 2
    }
    
    if(all(class(object) != "DA"))
    object$indY = length(object$names$blocks)+1  # we don't want to remove anything in that case, and 3 is higher than the number of blocks which is 2
    
    if(is.numeric(block))
    {
        if(any(block>length(object$names$blocks[-object$indY])))
        stop("'block' needs to be lower than the number of blocks in the fitted model, which is ",length(object$names$blocks)-1)
        
    }else if(is.character(block) & any(is.na(match(block,object$names$blocks[-object$indY])))) {
        stop("Incorrect value for 'block', 'block' should be among the blocks used in your object: ", paste(object$names$blocks[-object$indY],collapse=", "), call. = FALSE)
    }


    if (!missing(subtitle))
    {
        if (length(subtitle)!=length(block))
        stop("'subtitle' indicates the subtitle of the plot for each block and it needs to be the same length as 'block'.")
    }

    if(!missing(study))
    {
    #study needs to be either: from levels(object$study), numbers from 1:nlevels(study) or "global"
        if (any(!study%in%c(levels(object$study), "global")))
        stop("'study' must be one of 'object$study' or 'all'.")

        if (length(study)!=length(unique(study)))
        stop("Duplicate in 'study' not allowed")
    }

    # cex
    # --
    if (size.name <= 0)
    size.name = 0.7
    
    if (!missing(size.legend))
    {
        if(size.legend <= 0)
        size.legend = 0.8
    } else {
        size.legend = NULL
    }
    
    # contrib
    # --
    if(!missing(contrib))
    {
        if(length(contrib) > 1 | !all(contrib %in% c("min", "max")))
        stop("'contrib' must be either 'min' or 'max'")
        
    }
    
    # xlim
    #---
    if(!missing(xlim))
    {
        # check xlim, has to be a matrix with number of rows=number of blocks, or a vector of two values
        if(length(block) == 1 & !is.null(xlim))
        {
            if(length(xlim) !=2)
            stop("'xlim' must be a vector of length 2")
            
            xlim = matrix(xlim, nrow = 1)
        }

        if(length(block)>1 & !is.null(xlim))
        {
            if(is.matrix(xlim) && ( !nrow(xlim) %in%c(1, length(block))  | ncol(xlim) != 2 ))
            stop("'xlim' must be a matrix with ",length(block)," rows (length(block)) and 2 columns")
            
            if(is.vector(xlim))
            {
                if(length(xlim) !=2)
                stop("'xlim' must be a matrix with ",length(block)," rows (length(block)) and 2 columns")

                xlim = matrix(xlim, nrow = 1)
            }
            
            if(nrow(xlim) != length(block)) # we complete xlim to have one xlim per block
            xlim = matrix(rep(xlim, length(block)), nrow = length(block), byrow=TRUE)
        }
        
    } else {
        xlim = NULL
    }
    
    #names.var
    #-----
    if(!is.null(name.var))
    {
        if (length(block) >1 && length(block) != length(name.var))
        stop("'names' has to be a list of length the number of block to plot: ", length(block))
        
        if (length(block) > 1)
        {
            for (block_i in block)
            {
                if(length(name.var[[block_i]])!= nrow(object$loadings[[block_i]]))
                stop("For block '", block_i,"', 'name.var' should be a vector of length ", nrow(object$loadings[[block_i]]))
            }
        } else {
            if(length(name.var)!= nrow(object$loadings[[block]]))
            stop("For block '", block,"', 'name.var' should be a vector of length ", nrow(object$loadings[[block]]))

        }
    }
    #title
    #-----
    if (!is.null(title) & !is.character(title))
    warning('title needs to be of type character')

    #col
    #-----
    if (!is.null(col) & (length(col) !=  1))
    {
        warning('col must be the of length 1, by default set to default colors')
        col = color.mixo(1)  # by default set to the colors in color.mixo (10 colors)
    }
    if (is.null(col))
    col = color.mixo(1) # by default set to the colors in color.mixo (10 colors)


    return(list(col = col, size.name = size.name, size.legend = size.legend, block = block, xlim = xlim))
}

layout.plotLoadings = function(layout, plot, legend, block)
{
    # layout
    # --
    if(plot == TRUE)
    {
        opar = par(no.readonly = TRUE)
        reset.mfrow = FALSE # if set to TRUE, the algorithm ends up with  par(mfrow=reset.mfrow)
        nResp = length(block) + length(block) * legend  #number of blocks *2 if legend is plotted

        if (is.null(layout))
        {
            # check if there are enough plots in mfrow
            omfrow = par("mfrow")
            available.plots = prod(omfrow)
            if (available.plots<nResp) # if not enough plots available, we create our new plot
            {

                if (legend)
                {
                    nRows = min(c(2, ceiling(nResp / 4)))
                    nCols = min(c(4, nResp))
                    layout(matrix(1 : (nCols * nRows), nRows, nCols, byrow=TRUE),rep(c(0.7,0.7 -0.4*legend),nCols/(1+legend)))
                } else {
                    nRows = min(c(3, ceiling(nResp/3)))
                    nCols = min(c(3, ceiling(nResp / nRows)))
                    
                    layout(matrix(1 : (nCols * nRows), nRows, nCols, byrow=TRUE))

                }
                if (nRows * nCols < nResp)
                devAskNewPage(TRUE)
                
                reset.mfrow=TRUE # we changed mfrow to suits our needs, so we reset it at the end
            }
        } else {
            if (length(layout) != 2 || !is.numeric(layout) || any(is.na(layout)))
            stop("'layout' must be a numeric vector of length 2.")
            
            nRows = layout[1]
            nCols = layout[2]
            par(mfrow = layout)
            
            if (nRows * nCols < nResp)
            devAskNewPage(TRUE)
        }
        
    } else {
        reset.mfrow = FALSE
        opar = NULL
    }

    return(list(reset.mfrow = reset.mfrow, opar = opar))

}


get.loadings.ndisplay = function(object, comp, block, name.var, name.var.complete, ndisplay)
{
    ##selectvar
    selected.var = selectVar(object, comp = comp, block = block) # gives name and values of the blocks in 'block'
    name.selected.var = selected.var[[1]]$name
    value.selected.var = selected.var[[1]]$value
    
    # ndisplay
    # ------
    # if null set by default to all variables from selectVar
    if (is.null(ndisplay))
    {
        ndisplay.temp = length(name.selected.var)
    } else if (ndisplay > length(name.selected.var)) {
        message("'ndisplay' value is larger than the number of selected variables! It has been reseted to ", length(name.selected.var), " for block ", block)
        ndisplay.temp = length(name.selected.var)
    } else {
        ndisplay.temp = ndisplay
    }
    
    name.selected.var = name.selected.var[1:ndisplay.temp]
    value.selected.var = value.selected.var[1:ndisplay.temp,]
    
    #comp
    # ----
    if (any(class(object) %in% c("mixo_pls","mixo_spls", "rcc")))# cause pls methods just have 1 ncomp, block approaches have different ncomp per block
    {
        ncomp = object$ncomp
        object$X = list(X = object$X, Y = object$Y) # so that the data is in object$X, either it's a pls or block approach
    } else {
        ncomp = object$ncomp[block]
    }
    
    if (any(max(comp) > ncomp))
    stop(paste("Argument 'comp' should be less or equal to ", ncomp))
    
    names.block = as.character(names(selected.var)[1]) #it should be one block and ncomp, so we take the first one
    
    X = object$X[names.block][[1]]
    
    #name.var
    ind.match = match(name.selected.var, colnames(X)) # look at the position of the selected variables in the original data X
    if(!is.null(name.var))
    {
        if(length(name.var)!= ncol(X))
        stop("For block '", names.block,"', 'name.var' should be a vector of length ", ncol(X))
        
        colnames.X = as.character(name.var[ind.match]) # get the
    }else{
        colnames.X = as.character(colnames(X))[ind.match]
    }
    X = X[, name.selected.var, drop = FALSE] #reduce the problem to ndisplay
    
    #completing colnames.X by the original names of the variables when missing
    if (name.var.complete == TRUE)
    {
        ind = which(colnames.X == "")
        if (length(ind) > 0)
        colnames.X[ind] = colnames(X)[ind]
    }
    
    
    return(list(X = X, names.block = names.block, colnames.X = colnames.X, name.selected.var = name.selected.var, value.selected.var = value.selected.var))
}



get.contrib.df = function(Y, X, method, contrib, value.selected.var, colnames.X, name.selected.var, legend.color, col.ties)
{
    # Start: Initialisation
    which.comp = method.group = list()
    which.contrib = data.frame(matrix(FALSE, ncol = nlevels(Y) + 2, nrow = length(colnames.X),
    dimnames = list(name.selected.var, c(paste0("Contrib.", levels(Y)), "Contrib", "GroupContrib"))))
    # End: Initialisation
    
    # calculate the max.method per group for each variable, and identifies which group has the max max.method
    for(k in 1:ncol(X))
    {
        method.group[[k]] = tapply(X[, k], Y, method, na.rm=TRUE) #method is either mean or median
        # determine which group has the highest mean/median
        which.contrib[k, 1:nlevels(Y)] = (method.group[[k]]) == get(contrib)((method.group[[k]])) # contrib is either min or max
    }
    
    # we also add an output column indicating the group that is max
    # if ties, we set the color to white
    which.contrib$color = apply(which.contrib, 1, function(x)
    {
        if (length(which(x)) > 1)
        {
            return(col.ties)
        } else { # otherwise we use legend color provided
            return(legend.color[1 : nlevels(Y)][which(x)])
        }
    })
    
    which.contrib$GroupContrib = apply(which.contrib[, 1:(nlevels(Y))], 1, function(x)
    {
        if (length(which(x)) > 1)
        {
            return("tie")
        } else {
            return(levels(Y)[which(x)])
        }
    })
    
    method.group = do.call(rbind, method.group)
    df = data.frame(method.group, which.contrib, importance = value.selected.var)
    return(df)
}
