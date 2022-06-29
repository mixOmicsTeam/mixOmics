## ------------------------------------------------------------------------ ##
####                               Generic                                ####
## ------------------------------------------------------------------------ ##

#' Plot of Loading vectors
#' 
#' This function provides a horizontal bar plot to visualise loading vectors.
#' For discriminant analysis, it provides visualisation of highest or lowest
#' mean/median value of the variables with color code corresponding to the
#' outcome of interest.
#' 
#' 
#' The contribution of each variable for each component (depending on the
#' object) is represented in a barplot where each bar length corresponds to the
#' loading weight (importance) of the feature. The loading weight can be
#' positive or negative.
#' 
#' For discriminant analysis, the color corresponds to the group in which the
#' feature is most 'abundant'.  Note that this type of graphical output is
#' particularly insightful for count microbial data - in that latter case using
#' the \code{method = 'median'} is advised. Note also that if the parameter
#' \code{contrib} is not provided, plots are white.
#' 
#' For MINT analysis, \code{study="global"} plots the global loadings while
#' partial loadings are plotted when \code{study} is a level of
#' \code{object$study}. Since variable selection in MINT is performed at the
#' global level, only the selected variables are plotted for the partial
#' loadings even if the partial loadings are not sparse. See references.
#' Importantly for multi plots, the legend accounts for one subplot in the
#' layout design.
#' 
#' @aliases plotLoadings.pls plotLoadings.spls
#' @param object object
#' @param contrib a character set to 'max' or 'min' indicating if the color of
#' the bar should correspond to the group with the maximal or minimal
#' expression levels / abundance.
#' @param method a character set to 'mean' or 'median' indicating the criterion
#' to assess the contribution. We recommend using median in the case of count
#' or skewed data.
#' @param study Indicates which study are to be plotted.  A character vector
#' containing some levels of \code{object$study}, "all.partial" to plot all
#' studies or "global" is expected.
#' @param block A single value indicating which block to consider in a
#' \code{sgccda} object.
#' @param comp integer value indicating the component of interest from the
#' object.
#' @param col color used in the barplot, only for object from non Discriminant
#' analysis
#' @param plot Logical indicating of the plot should be output. If set to FALSE
#' the user can extract the contribution matrix, see example. Default value is
#' TRUE.
#' @param show.ties Logical. If TRUE then tie groups appear in the color set by
#' \code{col.ties}, which will appear in the legend. Ties can happen when
#' dealing with count data type. By default set to TRUE.
#' @param col.ties Color corresponding to ties, only used if
#' \code{show.ties=TRUE} and ties are present.
#' @param ndisplay integer indicating how many of the most important variables
#' are to be plotted (ranked by decreasing weights in each PLS-component).
#' Useful to lighten a graph.
#' @param size.name A numerical value giving the amount by which plotting the
#' variable name text should be magnified or reduced relative to the default.
#' @param size.legend A numerical value giving the amount by which plotting the
#' legend text should be magnified or reduced relative to the default.
#' @param name.var A character vector indicating the names of the variables.
#' The names of the vector should match the names of the input data, see
#' example.
#' @param name.var.complete Logical. If \code{name.var} is supplied with some
#' empty names, \code{name.var.complete} allows you to use the initial variable
#' names to complete the graph (from colnames(X)). Defaut to FALSE.
#' @param title A set of characters to indicate the title of the plot. Default
#' value is NULL.
#' @param subtitle subtitle for each plot, only used when several \code{block}
#' or \code{study} are plotted.
#' @param size.title size of the title
#' @param size.subtitle size of the subtitle
#' @param legend Logical indicating if the legend indicating the group outcomes
#' should be added to the plot. Default value is TRUE.
#' @param legend.color A color vector of length the number of group outcomes.
#' See examples.
#' @param legend.title A set of characters to indicate the title of the legend.
#' Default value is NULL.
#' @param layout Vector of two values (rows,cols) that indicates the layout of
#' the plot. If \code{layout} is provided, the remaining empty subplots are
#' still active
#' @param border Argument from \code{\link{barplot}}: indicates whether to draw
#' a border on the barplot.
#' @param xlim Argument from \code{\link{barplot}}: limit of the x-axis. When
#' plotting several \code{block}, a matrix is expected where each row is the
#' \code{xlim} used for each of the blocks.
#' @param \dots not used.
#' @return Invisibly returns a \code{data.frame} containing the contribution of 
#' features on each component. For supervised models the contributions for
#'  each class is also specified. See details.
#' @author Florian Rohart, Kim-Anh Lê Cao, Benoit Gautier, Al J Abadi
#' @seealso \code{\link{pls}}, \code{\link{spls}}, \code{\link{plsda}},
#' \code{\link{splsda}}, \code{\link{mint.pls}}, \code{\link{mint.spls}},
#' \code{\link{mint.plsda}}, \code{\link{mint.splsda}},
#' \code{\link{block.pls}}, \code{\link{block.spls}},
#' \code{\link{block.plsda}}, \code{\link{block.splsda}},
#' \code{\link{mint.block.pls}}, \code{\link{mint.block.spls}},
#' \code{\link{mint.block.plsda}}, \code{\link{mint.block.splsda}}
#' @references Rohart F. et al (2016, submitted). MINT: A multivariate
#' integrative approach to identify a reproducible biomarker signature across
#' multiple experiments and platforms.
#' 
#' Eslami, A., Qannari, E. M., Kohler, A., and Bougeard, S. (2013). Multi-group
#' PLS Regression: Application to Epidemiology. In New Perspectives in Partial
#' Least Squares and Related Methods, pages 243-255. Springer.
#' 
#' Singh A., Shannon C., Gautier B., Rohart F., Vacher M., Tebbutt S.
#' and Lê Cao K.A. (2019), DIABLO: an integrative approach for identifying key 
#' molecular drivers from multi-omics assays, Bioinformatics, 
#' Volume 35, Issue 17, 1 September 2019, Pages 3055–3062.
#' 
#' Lê Cao, K.-A., Martin, P.G.P., Robert-Granie, C. and Besse, P. (2009).
#' Sparse canonical methods for biological data integration: application to a
#' cross-platform study. \emph{BMC Bioinformatics} \bold{10}:34.
#' 
#' Tenenhaus, M. (1998). \emph{La regression PLS: theorie et pratique}. Paris:
#' Editions Technic.
#' 
#' Wold H. (1966). Estimation of principal components and related models by
#' iterative least squares. In: Krishnaiah, P. R. (editors), \emph{Multivariate
#' Analysis}. Academic Press, N.Y., 391-420.
#' @keywords multivariate
#' @export
#' @example ./examples/plotLoadings-examples.R
plotLoadings  <-
    function(object, ...)
        UseMethod("plotLoadings")

## ------------------------------------------------------------------------ ##
####                               Methods                                ####
## ------------------------------------------------------------------------ ##

## --------------------- PLS, sPLS, rCC, rGCCA, sGCCA --------------------- ##

#' @rdname plotLoadings
#' @method plotLoadings mixo_pls
#' @export
plotLoadings.mixo_pls <- 
    function(object,
             block,
             comp = 1,
             col = NULL,
             ndisplay = NULL,
             size.name = 0.7,
             name.var = NULL,
             name.var.complete = FALSE,
             title = NULL,
             subtitle,
             size.title = rel(2),
             size.subtitle = rel(1.5),
             layout = NULL,
             border = NA,
             xlim = NULL,
             ...
             
    )
    {
        
        # -- input checks
        check = check.input.plotLoadings(object = object, block = block, subtitle = subtitle, size.name = size.name, title = title, col = col, name.var = name.var, xlim = xlim)
        
        col = check$col
        size.name = check$size.name
        block = check$block
        xlim = check$xlim
        
        # -- layout
        res = layout.plotLoadings(layout = layout, plot = TRUE, legend = FALSE, block = block)
        reset.mfrow = res$reset.mfrow
        opar = res$opar
        omar = par("mar") #reset mar at the end
        
        if (length(block) == 1 & !is.null(name.var))
            name.var = list(name.var = name.var)
        
        contrib.df <- list()
        
        for (i in 1 : length(block))
        {
            res = get.loadings.ndisplay(object = object, comp = comp, block = block[i], name.var = name.var[[i]], name.var.complete = name.var.complete, ndisplay = ndisplay)
            X = res$X
            names.block = res$names.block
            colnames.X = res$colnames.X
            value.selected.var = res$value.selected.var
            
            df = data.frame(importance = value.selected.var) # contribution of the loading
            
            # barplot with contributions
            if (!is.null(title) & length(block) > 1)
            {
                par(mar = c(4, max(7, max(sapply(colnames.X, nchar),na.rm = TRUE)/3), 6, 2))
            } else {
                par(mar = c(4, max(7, max(sapply(colnames.X, nchar),na.rm = TRUE)/3), 4, 2))
            }
            
            .plotLoadings_barplot(height = df$importance, col = col, names.arg = colnames.X, cex.name = size.name, border = border, xlim = xlim[i, ])
            
            if ( (length(block) == 1 & is.null(title)) | (length(block) > 1 & missing(subtitle)))
            {
                title(paste0('Loadings on comp ', comp, "\nBlock '", names.block,"'"), line=0, cex.main = size.title)
            } else if (length(block) == 1) {
                title(paste(title), line=0, cex.main = size.title)
            } else if (length(block) > 1 & !missing(subtitle)) {
                title(paste(subtitle[i]), line=0, cex.main = size.subtitle)
            }
            
            contrib.df <- c(contrib.df, list(df))
        }
        
        names(contrib.df) <- block
        
        if (length(block) > 1 & !is.null(title))
            title(title, outer=TRUE, line = -2, cex.main = size.title)
        
        if (reset.mfrow)
            par(opar)#par(mfrow = c(1,1))
        
        par(mar = omar) #reset mar
        
        # return the contribution matrix
        return(invisible(contrib.df))
    }

#' @rdname plotLoadings
#' @method plotLoadings mixo_spls
#' @export
plotLoadings.mixo_spls <- plotLoadings.mixo_pls

#' @rdname plotLoadings
#' @method plotLoadings rcc
#' @export
plotLoadings.rcc <- plotLoadings.mixo_pls

#' @rdname plotLoadings
#' @method plotLoadings sgcca
#' @export
plotLoadings.sgcca <- plotLoadings.mixo_pls

#' @rdname plotLoadings
#' @method plotLoadings rgcca
#' @export
plotLoadings.rgcca <- plotLoadings.mixo_pls

## ------------------------ PCA, sPCA, IPCA, sIPCA ------------------------ ##

#' @rdname plotLoadings
#' @method plotLoadings pca
#' @export
plotLoadings.pca <-
    function(object,
             comp = 1,
             col = NULL,
             ndisplay = NULL,
             size.name = 0.7,
             name.var = NULL,
             name.var.complete = FALSE,
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
        par(mar = c(4, max(7, max(sapply(colnames.X, nchar),na.rm = TRUE)/3), 4, 2))
        
        .plotLoadings_barplot(height = df$importance, col = col, names.arg = colnames.X, cex.name = size.name, border = border, xlim = xlim)
        
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

## ----------------------- PLS-DA, SPLS-DA, sGCCDA ------------------------ ##

#' @rdname plotLoadings
#' @method plotLoadings mixo_plsda
#' @export
plotLoadings.mixo_plsda <- 
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
            
            
            contrib.df <- list()
            
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
                if (plot) # condition if all we need is the contribution stats
                {
                    if (!is.null(title) & length(block) > 1)
                    {
                        par(mar = c(4, max(7, max(sapply(colnames.X, nchar),na.rm = TRUE)/3), 6, 2))
                    } else {
                        par(mar = c(4, max(7, max(sapply(colnames.X, nchar), na.rm = TRUE)/3), 4, 2))
                    }
                    
                    .plotLoadings_barplot(height = df$importance, col = df$color, names.arg = colnames.X, cex.name = size.name, border = border, xlim = xlim[i, ])
                    
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
                    
                    contrib.df <- c(contrib.df, list(df))
                } # end if plot
            }
            
            names(contrib.df) <- block
            
            
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
            return(invisible(contrib.df)) # df
        }# end contrib missing
    }


#' @rdname plotLoadings
#' @method plotLoadings mixo_splsda
#' @export
plotLoadings.mixo_splsda <- plotLoadings.mixo_plsda

#' @rdname plotLoadings
#' @method plotLoadings sgccda
#' @export
plotLoadings.sgccda <- plotLoadings.mixo_plsda

## ----------------------------- MINT.(s)PLS ------------------------------ ##

#' @rdname plotLoadings
#' @method plotLoadings mint.pls
#' @export
plotLoadings.mint.pls <-
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
            plotLoadings.mixo_pls(object = object, block = c("X", "Y"), comp = comp, ndisplay = ndisplay,
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
            check = check.input.plotLoadings(object = object, block = c("X", "Y"), study = study, title = title, col = col, size.name = size.name, name.var = name.var)
            
            col = check$col
            size.name = check$size.name
            block = check$block # c("X", "Y")
            
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
            
            if (!missing(subtitle))
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
                .plotLoadings_barplot(height = df$importance, col = df$color, names.arg = colnames.X, cex.name = size.name, border = border, xlim = xlim)
                
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

#' @rdname plotLoadings
#' @method plotLoadings mint.spls
#' @export
plotLoadings.mint.spls <- plotLoadings.mint.pls

## ---------------------------- MINT.(s)PLSDA ----------------------------- ##

#' @rdname plotLoadings
#' @method plotLoadings mint.plsda
#' @export
plotLoadings.mint.plsda <- 
    function(object,
             contrib = NULL,  # choose between 'max" or "min", NULL does not color the barplot
             method = "mean", # choose between 'mean" or "median"
             study = "global",
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
        
        # what I want is to modify the input and call plotLoadings.pls and plotLoadings.splsda where blocks are now studies
        # do not forget to change object$names$block in levels(object$study) and it should work, see you tomorrow
        
        if(any(study == "global"))
        {
            plotLoadings.mixo_plsda(object = object, contrib = contrib, method = method, block = "X", comp = comp, ndisplay = ndisplay,
                                    size.name = size.name,
                                    size.legend = size.legend,
                                    name.var = name.var,
                                    name.var.complete = name.var.complete,
                                    legend = legend,
                                    legend.color = legend.color,
                                    title = if(!is.null(title)){title}else{paste0('Contribution on comp ', comp, "\n All studies")},
                                    subtitle = subtitle,
                                    legend.title = legend.title,
                                    plot = plot,
                                    xlim = xlim,
                                    layout = layout,
                                    size.title = size.title,
                                    size.subtitle = size.subtitle,
                                    border = border,
                                    col.ties = col.ties)
            
        } else {
            
            # -- input checks
            check = check.input.plotLoadings(object = object, block = "X", size.name = size.name, size.legend = size.legend,
                                             title = title, col = NULL, name.var = name.var, contrib = contrib)
            
            size.name = check$size.name
            size.legend = check$size.legend
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
            
            if (!missing(subtitle))
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
            
            X.study = study_split(X, study = object$study)
            Y = object$Y #v6: all $Y are factors for DA methods
            Y.study = study_split(Y, study = object$study)
            
            df.final = list()
            for (i in 1 : length(block))
            {
                
                value.selected.var =  object$loadings.partial$X [[block[i]]][, comp] [name.selected.var]
                
                
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
                
                if(!is.null(contrib))
                {
                    df = get.contrib.df(Y = factor(Y.study[[block[i]]]), X = X.study[[block[i]]], method = method, contrib = contrib, value.selected.var = value.selected.var, colnames.X = colnames.X, name.selected.var = name.selected.var, legend.color = legend.color, col.ties = col.ties)#data.frame(method.group, which.contrib, importance = value.selected.var)
                    # when working with sparse counts in particular and using the median to measure contribution
                    # ties to determine the contribution of a variable may happen, in that case remove them, otherwise they are showns as blank
                    if (show.ties == FALSE)
                    {
                        df = df[!df$color %in% col.ties, ]
                        colnames.X = rownames(df)
                    }
                    
                } else {
                    # if contrib is NULL, then we plot the loadings without colors
                    df = data.frame(importance = value.selected.var, color = "white", stringsAsFactors = FALSE) # contribution of the loading
                    border = TRUE
                }
                #  determine the colors/groups matching max contribution
                
                
                #display barplot with names of variables
                #added condition if all we need is the contribution stats
                if (!is.null(title) & length(block) > 1)
                {
                    par(mar = c(4, max(7, max(sapply(colnames.X, nchar),na.rm = TRUE)/2), 6, 2))
                } else {
                    par(mar = c(4, max(7, max(sapply(colnames.X, nchar),na.rm = TRUE)/2), 4, 2))
                }
                
                .plotLoadings_barplot(height = df$importance, col = df$color, names.arg = colnames.X, cex.name = size.name, border = border, xlim = xlim[i, ])
                
                if ( length(block) == 1 & is.null(title) )
                {
                    title(paste0('Contribution on comp ', comp, "\nStudy '", block[i],"'"), line=0, cex.main = size.title)
                } else if (length(block) == 1) {
                    title(paste(title), line=0, cex.main= size.title)
                } else if ((length(block) > 1 & missing(subtitle))) {
                    title(paste0('Contribution on comp ', comp, "\nStudy '", block[i],"'"), line=0, cex.main = size.subtitle)
                } else if (length(block) > 1 & !missing(subtitle)) {
                    title(paste(subtitle[i]), line=0, cex.main = size.subtitle)
                }
                
                if (legend)
                {
                    par(mar = c(5, 0, 4, 3) + 0.1)
                    plot(1,1, type = "n", axes = FALSE, ann = FALSE)
                    legend(0.8, 1, col = legend.color[1:nlevels(Y)], legend = levels(Y), pch = 19,
                           title = paste(legend.title),
                           cex = size.legend)
                }
                
                df.final[[i]] = df
            }
            names(df.final) = block
            
            # legend
            if (length(block) > 1 & !is.null(title))
                title(title, outer=TRUE, line = -2, cex.main = size.title)
            
            if (reset.mfrow)
                par(opar)#par(mfrow = omfrow)
            
            par(mar = omar) #reset mar
            
            return(invisible(df.final))
            
        }
        
    }

#' @rdname plotLoadings
#' @method plotLoadings mint.splsda
#' @export
plotLoadings.mint.splsda <- plotLoadings.mint.plsda

## ------------------------------------------------------------------------ ##
####                               Helpers                                ####
## ------------------------------------------------------------------------ ##

check.input.plotLoadings <- function(object,
                                     block,
                                     study,
                                     subtitle,
                                     size.name,
                                     size.legend,
                                     title,
                                     col,
                                     contrib,
                                     name.var,
                                     xlim)
{
    
    if (is.null(object$loadings))
        stop("'plotLoadings' should be used on object for which object$loadings is present.")
    
    # block
    # --
    if (missing(block))
    {
        if (!inherits(object, "DA"))
        {
            block = object$names$blocks
        } else  if (inherits(object, c("mixo_plsda", "mixo_splsda"))) {
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
    
    if (inherits(object, c("mixo_plsda", "mixo_splsda")) & (!all(block %in% c(1,"X")) | length(block) > 1 ))
        stop("'block' can only be 'X' or '1' for plsda and splsda object")
    
    if (inherits(object, c("mixo_plsda", "mixo_splsda","pca")))
    {
        object$indY = 2
    } else if (inherits(object, c("mixo_pls", "mixo_spls"))) {
        object$indY = 3 # we don't want to remove anything in that case, and 3 is higher than the number of blocks which is 2
    }
    
    if(!inherits(object, "DA"))
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

layout.plotLoadings <- function(layout, plot, legend, block)
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


get.loadings.ndisplay <- function(object,
                                  comp,
                                  block,
                                  name.var,
                                  name.var.complete,
                                  ndisplay)
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
    if (inherits(object, c("mixo_pls","mixo_spls", "rcc")))# cause pls methods just have 1 ncomp, block approaches have different ncomp per block
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



get.contrib.df <- function(Y,
                           X,
                           method,
                           contrib,
                           value.selected.var,
                           colnames.X,
                           name.selected.var,
                           legend.color,
                           col.ties)
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
        which.contrib[k, 1:nlevels(Y)] = (method.group[[k]]) == get(contrib)((method.group[[k]]), na.rm=TRUE) # contrib is either min or max
        
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
