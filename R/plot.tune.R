#############################################################################################################
# Authors:
#   Kim-Anh Le Cao, French National Institute for Agricultural Research and ARC Centre of Excellence ins Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Leigh Coonan, Queensland Faculty for Advanced Bioinformatics, Australia
#
# created: 20-08-2016
# last modified: 25-08-2016
#
# Copyright (C) 2010
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


plot.tune.splsda = plot.tune.spls = #plot.spca <- plot.ipca <- plot.sipca <-
function(x, optimal = TRUE, sd = TRUE, legend.position = "topright", col, ...)
{
    
    if (!is.logical(optimal))
    stop("'optimal' must be logical.", call. = FALSE)


    error <- x$error.rate
    if(sd & !is.null(x$error.rate.sd))
    {
        error.rate.sd = x$error.rate.sd
        ylim = range(c(error + error.rate.sd), c(error - error.rate.sd))
    } else {
        error.rate.sd = NULL
        ylim = range(error)
    }

    select.keepX <- x$choice.keepX[colnames(error)]
    comp.tuned = length(select.keepX)
    
    legend=NULL
    measure = x$measure
    
    if (length(select.keepX) < 10)
    {
        #only 10 colors in color.mixo
        if(missing(col))
        col = color.mixo(1:comp.tuned)
    } else {
        #use color.jet
        if(missing(col))
        col = color.jet(comp.tuned)
    }
    if(length(col) != comp.tuned)
    stop("'col' should be a vector of length ", comp.tuned,".")

    if(measure == "overall")
    {
         ylab = "Classification error rate"
    } else if (measure == "BER")
    {
        ylab = "Balanced error rate"
    } else if (measure == "MSE"){
        ylab = "MSE"
    }else if (measure == "MAE"){
        ylab = "MAE"
    }else if (measure == "Bias"){
        ylab = "Bias"
    }else if (measure == "R2"){
        ylab = "R2"
    }else if (measure == "AUC"){
        ylab = "AUC"
    }

    matplot(rownames(error),error, type = "l", axes = TRUE, lwd = 2, lty = 1, log = "x",
    xlab = "Number of selected features", ylab = ylab,
    col = col, ylim = ylim)
    
    if(optimal)
    {
        for(i in 1:comp.tuned)
        {
            # store coordinates of chosen keepX
            index = which(rownames(error) == select.keepX[i])
            # choseen keepX:
            points(rownames(error)[index], error[index,i], col = col[i], lwd=2, cex=3, pch = 18)
        }
    }
    
    if(!is.null(error.rate.sd))
    {
        for(j in 1:ncol(error))
        plot_error_bar(x = as.numeric(names(error[, j])), y =error[, j] , uiw=error.rate.sd[, j], add=TRUE, col = rep(col[j],each=nrow(error)))#, ...)
    }



    if(length(x$choice.keepX) == 1) #only first comp tuned
    {
        legend = "comp1"
    } else if(length(x$choice.keepX) == comp.tuned) # all components have been tuned
    {
        legend = c("comp1", paste("comp1 to", colnames(error)[-1]))
    } else { #first component was not tuned
        legend = paste("comp1 to", colnames(error))
    }

    legend(legend.position, lty = 1, lwd = 2, horiz = FALSE, col = col,
    legend = legend)
    
}



plot.tune.block.splsda =
function(x, sd = TRUE, col, ...)
{
    
    # R check
    error.sd=NULL
    
    error <- x$error.rate
    if(sd & !is.null(x$error.rate.sd))
    {
        error.rate.sd = x$error.rate.sd
        ylim = range(c(error + error.rate.sd), c(error - error.rate.sd))
    } else {
        error.rate.sd = NULL
        ylim = range(error)
    }
    select.keepX <- x$choice.keepX
    comp.tuned = length(select.keepX[[1]])
    
    if (length(select.keepX) < 10)
    {
        #only 10 colors in color.mixo
        if(missing(col))
        col = color.mixo(1:comp.tuned)
    } else {
        #use color.jet
        if(missing(col))
        col = color.jet(comp.tuned)
    }
    if(length(col) != comp.tuned)
    stop("'col' should be a vector of length ", comp.tuned,".")
    
    legend=NULL
    measure = x$measure


    if(measure == "overall")
    {
        ylab = "Classification error rate"
    } else if (measure == "BER")
    {
        ylab = "Balanced error rate"
    }

    if(FALSE)
    {
    # not ordered graph
    
    # creating one dataframe with all the comp
    error.plot = data.frame(comp = rep(colnames(error), each = nrow(error)), names = do.call("rbind", as.list(rownames(error))), error = do.call("rbind", as.list(error)), error.sd = do.call("rbind", as.list(error.rate.sd)), color = rep(col, each = nrow(error)))
    
    #    p = ggplot(error.plot, aes(x=reorder(names, -error), y=error)) +
    p = ggplot(error.plot, aes(x=names, y=error)) +
    geom_bar(stat="identity", fill = error.plot$color)
    if(sd) p = p + geom_errorbar(aes(ymin=error-error.sd, ymax = error+error.sd), width=0.2)
    
    p= p +
    ylab(ylab)+
    xlab("Number of selected features for each block")+
    coord_flip()+
    facet_grid(~comp,scales='free')
    p
    }
    
    pp=list()
    for(comp in 1:comp.tuned)
    {
        # order error per comp
        so = sort(error[,comp], index.return=TRUE, decreasing = TRUE)
        
        error.ordered = so$x
        error.sd.ordered = error.rate.sd[so$ix,comp]
        
        error.plot = data.frame (names = names(error.ordered), error = error.ordered, error.sd = error.sd.ordered, color = col[comp])
        
        ## ggplot
        p = ggplot(error.plot, aes(x=reorder(names, -error), y=error)) +
        geom_bar(stat="identity", fill = error.plot$color)
        if(sd) p = p + geom_errorbar(aes(ymin=error-error.sd, ymax = error+error.sd), width=0.2)
        
        p= p +
        ylab(ylab)+
        xlab("Number of selected features for each block")+
        ggtitle(colnames(error)[comp])+
        coord_flip()
        
        
        if(comp==1)
        p1=p
        if(comp==2)
        p2=p
        #+theme(axis.text.x = element_text(angle = 90, hjust = 1))
        
        pp[[comp]] = p#assign(paste0("p", colnames(error)[comp]), p)
        
    }
    
    do.call("grid.arrange", c(pp, nrow=ceiling(comp.tuned/3)))
    
    
}
