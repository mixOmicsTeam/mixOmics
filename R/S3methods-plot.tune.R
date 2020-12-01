## ------------------------------------------------------------------------ ##
###                               plot.tune                                ###
## ------------------------------------------------------------------------ ##
#' Plot model performance
#' 
#' Function to plot performance criteria, such as classification error rate or
#' correlation of cross-validated components for different models.
#' 
#' \code{plot.tune.splsda} plots the classification error rate or the balanced
#' error rate from x$error.rate, for each component of the model. A lozenge
#' highlights the optimal number of variables on each component.
#' 
#' \code{plot.tune.block.splsda} plots the classification error rate or the
#' balanced error rate from x$error.rate, for each component of the model. The
#' error rate is ordered by increasing value, the yaxis shows the optimal
#' combination of keepX at the top (e.g. `keepX on block 1'_`keepX on block
#' 2'_`keepX on block 3')
#' 
#' \code{plot.tune.spls} plots either the correlation of cross-validated
#' components or the Residual Sum of Square (RSS) values for these components
#' against those from the full model for both \code{t} (X components) and
#' \code{u} (Y components). The optimal number of features chosen are indicated
#' by squares.
#' 
#' \code{plot.tune.spca} plots  the correlation of cross-validated components from
#' the \code{tune.spca} function with respect to the full model.
#' @inheritParams plotIndiv
#' @param x a \code{tune} object. See details for supported objects.
#' @param measure Character. Measure used for plotting a \code{tune.spls} object.
#' One of c('cor', 'RSS').
#' @param optimal If TRUE, highlights the optimal keepX per component
#' @param sd If \eqn{nrepeat >= 3} was used in the call, error bar
#' shows the standard deviation if sd=TRUE. Note that the values might exceeed
#' the valid performance measures (such as [0, 1] for accuracy)
#' @param col character (or symbol) color to be used, possibly vector. One
#' colour per component.
#' @param \dots Further arguments sent to \code{\link{xyplot}} function.
#' Not currently used for \code{tune.spca}.
#' @return none
#' @author Kim-Anh Lê Cao, Florian Rohart, Francois Bartolo, Al J Abadi
#' @seealso \code{\link{tune.mint.splsda}}, \code{\link{tune.splsda}},
#'   \code{\link{tune.block.splsda}}, \code{\link{tune.spca}} and
#'   http://www.mixOmics.org for more details.
#' @keywords regression multivariate hplot
#' @name plot.tune
#' @example ./examples/plot.tune-examples.R
NULL
## --------------------------- plot.tune.(s)pls --------------------------- ##
#' @method plot tune.spls
#' @rdname plot.tune
#' @export
plot.tune.spls <-
    function(x, measure = NULL, pch = 16, cex = 1.2,...)
    {
        ncomp <- x$call$ncomp

        ## if measure not given, use object's 'measure.tune' for spls
        if (is.null(measure) & is(x, 'tune.spls') )
        {
            measure <- x$call$measure.tune
        } else {
            measure <- match.arg(measure, c('cor', 'RSS'))    
        }
        
        ggplot_measure <- function(x, v = c('u', 't'), title = NULL, measure = 'cor', ncomp) {
            
            pred <- ifelse(measure == 'cor', 'cor.pred', 'RSS.pred')
            
            df_comps <- lapply(seq_len(ncomp), function(comp)
            {
            ut <- lapply(c(u='u', t='t'), function(o){
                mean = x[[pred]][[o]][[comp]]$mean
                sd = x[[pred]][[o]][[comp]]$sd
                list(mean = round(mean, 2), sd = round(sd, 3))
            })
            
            df.list <- lapply(v, function(V) {
                df <- expand.grid(keepX = x$call$test.keepX, keepY = x$call$test.keepY)
                df$V <- V
                df$mean <- as.vector(ut[[V]]$mean)
                df$sd <- as.vector(ut[[V]]$sd)
                df
            })
            
            df <- Reduce(f = rbind, df.list)
            df$comp <- paste0('comp_', comp)
            df
        })
            df <- Reduce(f = rbind, df_comps)
            ## optimal keepX/keepY
            df$optimal <-              df$comp == 'comp_1' & df$keepX == x$choice.keepX[1] & df$keepY == x$choice.keepY[1]
            df$optimal <- df$optimal | df$comp == 'comp_2' & df$keepX == x$choice.keepX[2] & df$keepY == x$choice.keepY[2]
            text.size = as.integer(cex*10)
            p <- ggplot(df, aes(factor(keepX), factor(keepY))) + 
                geom_point(aes_string(size = 'mean', col = 'sd'), shape = pch) + 
                scale_color_gradient(low = 'blue', high = 'red', na.value = color.mixo(1)) + 
                theme(panel.border = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(size = 0.5, linetype = "solid",
                                               colour = "black"),
                      
                      panel.background = element_rect(fill='grey97'),
                      
                      axis.text = element_text( size = text.size ),
                        axis.text.x = element_text( size = text.size, angle = 90, hjust = 1),
                        axis.title = element_text( size = text.size),
                        legend.text = element_text( size = text.size ),
                        legend.title =  element_text( size = text.size),
                        plot.title = element_text(hjust = 0.5),
                        # subtitles
                        strip.text = element_text(size = 1.3*text.size, face = 'bold')
                      
                      )
                
            ## optimal keepX/keepY
                p <- p + geom_point(data = df[df$optimal,], 
                                    aes(factor(keepX), factor(keepY), size = mean*1.3), 
                                    shape = 0, 
                                    col = 'green', 
                                    show.legend = FALSE) +
                
                labs(x = 'keepX', y = 'keepY', size = 'mean', col = 'SD', 
                     title = sprintf("measure = '%s'", measure)) +
                facet_grid(V~comp)
            
            if (measure == 'RSS')
            {
                p <- p + scale_size_continuous(range = c(6,1))
            }
            
            p <- p + guides(colour = guide_legend(order=2, override.aes = list(size=2)),
                            size = guide_legend(order=1))
            
            list(gg.plot = p, df= df)
        }
        
        res <- ggplot_measure(x=x, measure = measure, ncomp = ncomp)
        
        res$gg.plot
    }

## -------------------------- plot.tune.splsda -------------------------- ##
#' @name plot.tune
#' @method plot tune.splsda
#' @importFrom reshape2 melt
#' @export
plot.tune.splsda <-
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
## ------------------------ plot.tune.block.(s)plsda ---------------------- ##
#' @importFrom gridExtra grid.arrange
#' @rdname plot.tune
#' @method plot tune.block.splsda
#' @export
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
        
        # if(FALSE)
        # {
        #     # not ordered graph
        #     
        #     # creating one dataframe with all the comp
        #     error.plot = data.frame(comp = rep(colnames(error), each = nrow(error)), names = do.call("rbind", as.list(rownames(error))), error = do.call("rbind", as.list(error)), error.sd = do.call("rbind", as.list(error.rate.sd)), color = rep(col, each = nrow(error)))
        #     
        #     #    p = ggplot(error.plot, aes(x=reorder(names, -error), y=error)) +
        #     p = ggplot(error.plot, aes(x=names, y=error)) + 
        #         theme_minimal() +
        #         geom_bar(stat="identity", fill = error.plot$color)
        #     if(sd) p = p + geom_errorbar(aes(ymin=error-error.sd, ymax = error+error.sd), width=0.04)
        #     
        #     p= p +
        #         ylab(ylab)+
        #         xlab("Number of selected features for each block")+
        #         coord_flip()+
        #         facet_grid(~comp,scales='free')
        #     p
        # }
        
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
                theme_classic() +
                geom_bar(stat="identity", fill = error.plot$color)
            if(sd) p = p + geom_errorbar(aes(ymin=error-error.sd, ymax = error+error.sd), width=0.04)
            
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

## --------------------------- plot.tune.spca --------------------------- ##
#' @rdname plot.tune
#' @method plot tune.spca
#' @export
plot.tune.spca <-
    function(x, optimal = TRUE, sd = NULL, col=NULL, ...)
    {
        ncomp <- length(x$cor.comp)
        nrepeat <- x$call$nrepeat
        
        if (nrepeat <= 2 & isTRUE(sd))
        {
            cat("nrepeat < 2 so no SD can be calculated.",
                "setting sd to FALSE")
            sd <- FALSE
        } else if (nrepeat > 2 & is.null(sd))
        {
            sd <- TRUE
        }
        
        cors <- mapply(z=x$cor.comp, w=seq_len(ncomp), FUN = function(z, w){
            z$comp = w
            rownames(z) <- NULL
            if (nrepeat > 2)
            {
                z$corQ1 = z$cor.mean - z$cor.sd
                z$corQ3 = z$cor.mean + z$cor.sd
            }
            z
        }, SIMPLIFY = FALSE)
        cors <- Reduce(rbind, cors)
        cors$comp <- factor(cors$comp)
        
        if (is.null(col))
        {
            col <- color.mixo(seq_len(ncomp))
        }
        names(col) <- seq_len(ncomp)
        p <- ggplot(cors, aes_string('keepX', 'cor.mean', col = 'comp')) +
            theme_minimal() +
            geom_line() +
            geom_point() +
            scale_x_continuous(trans='log10', breaks = cors$keepX) +
            ylim(c(min(cors$corQ1, 0), max(cors$corQ3, 1))) +
            labs(x= 'Number of features selected',
                 y = 'Correlation of components',
                 col = 'Comp') +
            scale_color_manual(values = col)
        
        if (nrepeat > 2) {
            p <- p +  geom_point(data=cors[!is.na(cors$opt.keepX),], 
                                 aes_string('keepX', 'cor.mean', col = 'comp'), 
                                 size=6, shape = 18, show.legend = FALSE)
        }
        
        
        ## ----- error bars
        if (isTRUE(sd))
        {
            p <- p + geom_errorbar(aes_string(ymin = 'corQ1', ymax = 'corQ3'), 
                                   # position = position_dodge(0.02),
                                   width = 0.04,
                                   ...)
            ## suppress "position_dodge requires non-overlapping x intervals"
            suppressWarnings(print(p))
            return(invisible(p))
        } else
        {
            return(p)
        }
        
    }
