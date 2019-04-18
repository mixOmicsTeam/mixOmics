###############################################################################
# Authors:
#   Kim-Anh Le Cao
#   Florian Rohart
#   Leigh Coonan
#
# created: 20-08-2016
# last modified: 29-01-2019
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
###############################################################################


plot.tune.splsda = plot.tune.spls = #plot.spca <- plot.ipca <- plot.sipca <-
function(x, optimal = TRUE, sd = TRUE, col, ...)
{
    # to satisfy R CMD check that doesn't recognise x, y and group (in aes)
    y = Comp = lwr = upr = NULL

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
        col = color.mixo(seq_len(comp.tuned))
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
    
    #legend
    names.comp = substr(colnames(error),5,10) # remove "comp" from the name
    if(length(x$choice.keepX) == 1){
         #only first comp tuned
         legend = "1"
    } else if(length(x$choice.keepX) == comp.tuned) {
        # all components have been tuned
        legend = c("1", paste("1 to", names.comp[-1]))
    } else {
        #first components were not tuned
        legend = paste("1 to", names.comp)
    }

    
    # creating data.frame with all the information
    df = melt(error)
    colnames(df) = c("x","Comp","y")
    df$Comp = factor(df$Comp, labels=legend)

    p = ggplot(df, aes(x = x, y = y, color = Comp)) +
    labs(x = "Number of selected features", y = ylab) +
    theme_bw() +
    geom_line()+ geom_point()
    p = p+ scale_x_continuous(trans='log10') +
    scale_color_manual(values = col)

    # error bar
    if(!is.null(error.rate.sd))
    {
        dferror = melt(error.rate.sd)
        df$lwr = df$y - dferror$value
        df$upr = df$y + dferror$value
        
        #adding the error bar to the plot
        p = p + geom_errorbar(data=df,aes(ymin=lwr, ymax=upr))
    }
    
    if(optimal)
    {
        index = NULL
        for(i in seq_len(comp.tuned))
            index = c(index, which(df$x == select.keepX[i] & df$Comp == levels(df$Comp)[i]))
            
        # adding the choseen keepX to the graph
        p=p + geom_point(data=df[index,],size=7, shape = 18)
        p = p + guides(color = guide_legend(override.aes =
        list(size=0.7,stroke=1)))
    }
    
    p
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
        col = color.mixo(seq_len(comp.tuned))
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
    for(comp in seq_len(comp.tuned))
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
