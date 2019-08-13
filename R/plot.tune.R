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

#' Plot for model performance
#' 
#' Function to plot performance criteria, such as classification error rate or
#' balanced error rate on a tune.splsda result.
#' 
#' \code{plot.tune.splsda} plots the classification error rate or the balanced
#' error rate from x$error.rate, for each component of the model. A lozenge
#' highlights the optimal number of variables on each component.
#' 
#' \code{plot.tune.block.splsda} plots the classification error rate or the
#' balanced error rate from x$error.rate, for each component of the model. The
#' error rate is ordered by increasing value, the yaxis shows the optimal
#' combination of keepX at the top.
#' (e.g. `'keepX on block 1'_'keepX on block 2'_'keepX on block 3'``)
#' 
#' @aliases plot.tune.block.splsda plot.tune.splsda plot.tune
#' @param x an \code{tune.splsda} object.
#' @param optimal If TRUE, highlights the optimal keepX per component
#' @param sd If 'nrepeat' was used in the call to 'tune.splsda', error bar
#' shows the standard deviation if sd=TRUE
#' @param col character (or symbol) color to be used, possibly vector. One
#' color per component.
#' @param \dots Further arguments sent to \code{\link{xyplot}} function.
#' @return none
#' @author Kim-Anh Lê Cao, Florian Rohart, Francois Bartolo.
#' @seealso \code{\link{tune.mint.splsda}}, \code{\link{tune.splsda}}
#' \code{\link{tune.block.splsda}} and http://www.mixOmics.org for more
#' details.
#' @keywords regression multivariate hplot
#' @examples
#' 
#' \dontrun{
#' ## validation for objects of class 'splsda'
#' 
#' data(breast.tumors)
#' X = breast.tumors$gene.exp
#' Y = as.factor(breast.tumors$sample$treatment)
#' out = tune.splsda(X, Y, ncomp = 3, nrepeat = 5, logratio = "none",
#' test.keepX = c(5, 10, 15), folds = 10, dist = "max.dist",
#' progressBar = TRUE)
#' 
#' 
#' plot(out)
#' plot(out, sd=FALSE)
#' 
#' 
#' \dontrun{
#' ## validation for objects of class 'mint.splsda'
#' 
#' data(stemcells)
#' data = stemcells$gene
#' type.id = stemcells$celltype
#' exp = stemcells$study
#' 
#' out = tune(method="mint.splsda", X=data,Y=type.id, ncomp=2, study=exp, test.keepX=seq(1,10,1))
#' out$choice.keepX
#' 
#' plot(out)
#' 
#' 
#' 
#' ## validation for objects of class 'mint.splsda'
#' 
#' data("breast.TCGA")
#' # this is the X data as a list of mRNA and miRNA; the Y data set is a single data set of proteins
#' data = list(mrna = breast.TCGA$data.train$mrna, mirna = breast.TCGA$data.train$mirna,
#' protein = breast.TCGA$data.train$protein)
#' # set up a full design where every block is connected
#' # could also consider other weights, see our mixOmics manuscript
#' design = matrix(1, ncol = length(data), nrow = length(data),
#' dimnames = list(names(data), names(data)))
#' diag(design) =  0
#' design
#' # set number of component per data set
#' ncomp = 5
#' 
#' # Tuning the first two components
#' # -------------
#' 
#' # definition of the keepX value to be tested for each block mRNA miRNA and protein
#' # names of test.keepX must match the names of 'data'
#' test.keepX = list(mrna = seq(10,40,20), mirna = seq(10,30,10), protein = seq(1,10,5))
#' 
#' # the following may take some time to run, note that for through tuning
#' # nrepeat should be > 1
#' tune = tune.block.splsda(X = data, Y = breast.TCGA$data.train$subtype,
#' ncomp = ncomp, test.keepX = test.keepX, design = design, nrepeat = 3)
#' 
#' tune$choice.ncomp
#' tune$choice.keepX
#' 
#' plot(tune)
#' }
#' }
#' 

#' @importFrom reshape2 melt
#' @rdname S3methods-plot.tune
#' @export
plot.tune.spls <- function(x, optimal = TRUE, sd = TRUE, col, ...)
{
    # to satisfy R CMD check that doesn't recognise x, y and group (in aes)
    y = Comp = lwr = upr = NULL

    if (!is.logical(optimal))
    stop("'optimal' must be logical.", call. = FALSE)


    error <- x$error.rate
    if (sd & !is.null(x$error.rate.sd))
    {
        error.rate.sd = x$error.rate.sd
        ylim = range(c(error + error.rate.sd), c(error - error.rate.sd))
    } else {
        error.rate.sd = NULL
        ylim = range(error)
    }

    select.keepX <- x$choice.keepX[colnames(error)]
    comp.tuned = length(select.keepX)
    
    legend = NULL
    measure = x$measure
    
    if (length(select.keepX) < 10)
    {
        #only 10 colors in color.mixo
        if (missing(col))
            col = color.mixo(seq_len(comp.tuned))
    } else {
        #use color.jet
        if (missing(col))
            col = color.jet(comp.tuned)
    }
    if (length(col) != comp.tuned)
        stop("'col' should be a vector of length ", comp.tuned, ".")
    
    if (measure == "overall")
    {
        ylab = "Classification error rate"
    } else if (measure == "BER")
    {
        ylab = "Balanced error rate"
    } else if (measure == "MSE") {
        ylab = "MSE"
    } else if (measure == "MAE") {
        ylab = "MAE"
    } else if (measure == "Bias") {
        ylab = "Bias"
    } else if (measure == "R2") {
        ylab = "R2"
    } else if (measure == "AUC") {
        ylab = "AUC"
    }

    #legend
    names.comp = substr(colnames(error), 5, 10) # remove "comp" from the name
    if (length(x$choice.keepX) == 1) {
        #only first comp tuned
        legend = "1"
    } else if (length(x$choice.keepX) == comp.tuned) {
        # all components have been tuned
        legend = c("1", paste("1 to", names.comp[-1]))
    } else {
        #first components were not tuned
        legend = paste("1 to", names.comp)
    }
    
    
    
    # creating data.frame with all the information
    df = melt(error)
    colnames(df) = c("x", "Comp", "y")
    df$Comp = factor(df$Comp, labels = legend)
    
    p = ggplot(df, aes(x = x, y = y, color = Comp)) +
        labs(x = "Number of selected features", y = ylab) +
        theme_bw() +
        geom_line() + geom_point()
    p = p + scale_x_continuous(trans = 'log10') +
        scale_color_manual(values = col)
    
    # error bar
    if (!is.null(error.rate.sd))
    {
        dferror = melt(error.rate.sd)
        df$lwr = df$y - dferror$value
        df$upr = df$y + dferror$value
        
        #adding the error bar to the plot
        p = p + geom_errorbar(data = df, aes(ymin = lwr, ymax = upr))
    }
    
    if (optimal)
    {
        index = NULL
        for (i in seq_len(comp.tuned))
            index = c(index, which(df$x == select.keepX[i] &
                                       df$Comp == levels(df$Comp)[i]))
        
        # adding the choseen keepX to the graph
        p = p + geom_point(data = df[index, ],
                           size = 7,
                           shape = 18)
        p = p + guides(color = guide_legend(override.aes =
                                                list(size = 0.7, stroke = 1)))
    }
    
    p
}

#' @rdname S3methods-plot.tune
#' @export
plot.tune.splsda <- plot.tune.spls

#' @rdname S3methods-plot.tune
#' @export
plot.tune.block.splsda <- function(x, sd = TRUE, col, ...)
{
    # R check
    error.sd = NULL
    
    error <- x$error.rate
    if (sd & !is.null(x$error.rate.sd))
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
        
        error.plot = data.frame (
            names = names(error.ordered),
            error = error.ordered,
            error.sd = error.sd.ordered,
            color = col[comp]
        )
        
        ## ggplot
        p = ggplot(error.plot, aes(x = reorder(names, -error), y = error)) +
            geom_bar(stat = "identity", fill = error.plot$color)
        if (sd)
            p = p + geom_errorbar(aes(ymin = error - error.sd, ymax = error + error.sd), width =
                                      0.2)
        
        p = p +
            ylab(ylab) +
            xlab("Number of selected features for each block") +
            ggtitle(colnames(error)[comp]) +
            coord_flip()
        
        
        if (comp == 1)
            p1 = p
        if (comp == 2)
            p2 = p
        #+theme(axis.text.x = element_text(angle = 90, hjust = 1))
        
        pp[[comp]] = p#assign(paste0("p", colnames(error)[comp]), p)
        
    }
    
    do.call("grid.arrange", c(pp, nrow=ceiling(comp.tuned/3)))
    
    
}

#' @rdname S3methods-plot.tune
#' @export
plot.tune.rcc <- function(x, col = heat.colors, ...) 
{
    
    opar = par(no.readonly = TRUE)
    
    grid1 = x$grid1
    grid2 = x$grid2
    mat = x$mat
    nlevel = min(255, length(c(unique(mat))))
    if (nlevel / 2 - trunc(nlevel / 2) == 0 & nlevel > 6) nlevel = nlevel + 1
    col = col(nlevel)
    
    def.par = par(no.readonly = TRUE)
    layout(matrix(c(2, 1), ncol = 2, nrow = 1, byrow = FALSE),
           widths = c(1, 0.18))
    
    #-- layout 1 --#
    min.mat = min(mat)
    max.mat = max(mat)
    par(pty = "m", mai = c(1.1, 0.1, 0.95, 0.5))
    z = seq(min.mat, max.mat, length = nlevel)
    breaks = seq(min.mat, max.mat, length = nlevel + 1)
    z = matrix(z, nrow = 1)
    
    image(z, col = col, 
          zlim = c(min.mat, max.mat), oldstyle = TRUE,
          xaxt = "n", yaxt = "n", breaks = breaks)
    box()
    
    par(usr = c(0, 1, min.mat, max.mat))
    binwidth = (max.mat - min.mat) / nlevel
    midpoints = seq(min.mat + binwidth/2, max.mat - binwidth/2, by = binwidth)
    
    if (nlevel <= 6)
    {
        axis(4, at = midpoints, labels = round(midpoints, 3))
    }else{
        binwidth = seq(1, nlevel, by = round(nlevel / 5))
        if (binwidth[length(binwidth)] != nlevel){
            binwidth[length(binwidth) + 1] = nlevel
        }
        
        axis(4, at = midpoints[binwidth], labels = round(midpoints[binwidth],
                                                         3))
        
    }
    
    
    
    #-- layout 2 --#
    par(mai = c(0.9, 0.85, 0.75, 0.2))
    image(grid1, grid2, mat, col = col,
          xlab = expression(lambda[1]), ylab=expression(lambda[2]),
          main = expression(CV(lambda[1], lambda[2])), axes = FALSE,
          zlim = c(min.mat, max.mat), oldstyle = TRUE)
    
    
    if (length(grid1) > 10) {
        grid1 = seq(min(grid1), max(grid1), length = 11)
    }
    if (length(grid2) > 10) {
        grid2 = seq(min(grid2), max(grid2), length = 11)
    }
    
    axis(1, at = grid1, labels = as.character(round(grid1, 4)))
    axis(2, at = grid2, labels = as.character(round(grid2, 4)))
    box()
    
    par(opar) #reset par
    
}


#' @rdname mixOmics-deprecated
#' @seealso \code{plot.tune.rcc}
plot.estim.regul <- function(x, col = heat.colors, ...) {
    .Deprecated('plot.tune.rcc', package = 'mixOmics')
    do.call(plot.tune.rcc, as.list(match.call())[-1])
}
