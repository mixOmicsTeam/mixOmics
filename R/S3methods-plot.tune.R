## ------------------------------------------------------------------------ ##
###                               plot.tune                                ###
## ------------------------------------------------------------------------ ##
#' Plot model performance
#' 
#' Function to plot performance criteria, such as classification error rate or
#' balanced error rate for different models.
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
#' \code{plot.tune.spca} plots  the correlation of cross-validated components from
#' the \code{tune.spca} function with respect to the full model.
#' 
#' @param x an \code{tune.*} object. See details for supported objects.
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
        
        ggplot_measure <- function(x, v = c('u', 't'), title = NULL, measure = 'cor') {
            
            pred <- ifelse(measure == 'cor', 'cor.pred', 'RSS.pred')
            
            ut <- lapply(c(u='u', t='t'), function(o){
                mean = x[[pred]][[o]][[ncomp]]$mean
                sd = x[[pred]][[o]][[ncomp]]$sd
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
            text.size = as.integer(cex*10)
            p <- ggplot(df, aes(factor(keepX), factor(keepY), size = mean, col = sd)) + 
                geom_point(shape = pch) + 
                scale_color_gradient(low = color.mixo(4), high = 'red') + 
                theme(panel.border = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(size = 0.5, linetype = "solid",
                                               colour = "black"),
                      
                      panel.background = element_rect(fill='grey95'),
                      
                      axis.text = element_text( size = text.size ),
                        axis.text.x = element_text( size = text.size ),
                        axis.title = element_text( size = text.size),
                        legend.text = element_text( size = text.size ),
                        legend.title =  element_text( size = text.size),
                        # subtitles
                        strip.text = element_text(size = 1.3*text.size, face = 'bold')
                      
                      ) +
                
                labs(x = 'keepX', y = 'keepY', size = measure) +
                facet_wrap(.~V)
        
            list(gg.plot = p, df= df)
        }
        
        res <- ggplot_measure(x=x, measure = measure)
        
        res$gg.plot
    }

#' @rdname plot.tune
#' @method plot tune.splsda
#' @export
plot.tune.splsda <- plot.tune.spls
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
