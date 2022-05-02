#----------------------------------------------------------------------------------------------------------#
#-- Includes plotIndiv for PLS, sPLS, PLS-DA, SPLS-DA,  --#
#----------------------------------------------------------------------------------------------------------#
#' @rdname plotIndiv
#' @method plotIndiv mixo_pls
#' @export
plotIndiv.mixo_pls <- 
    function(object,
             comp  = NULL,
             rep.space  = NULL,
             ind.names  = TRUE,
             group,
             col.per.group,
             style = "ggplot2",
             ellipse  = FALSE,
             ellipse.level  = 0.95,
             centroid = FALSE,
             star = FALSE,
             title = NULL,
             subtitle,
             legend = FALSE,
             X.label  = NULL,
             Y.label  = NULL,
             Z.label  = NULL,
             abline  = FALSE,
             xlim  = NULL,
             ylim  = NULL,
             col,
             cex,
             pch,
             pch.levels,
             alpha = 0.2,
             axes.box  = "box",
             layout = NULL,
             size.title = rel(2),
             size.subtitle = rel(1.5),
             size.xlabel = rel(1),
             size.ylabel = rel(1),
             size.axis = rel(0.8),
             size.legend = rel(1),
             size.legend.title = rel(1.1),
             legend.title = "Legend",
             legend.title.pch = "Legend",
             legend.position = "right",
             point.lwd = 1,
             background = NULL,
             ...
             
    )
    {
        plot_parameters = list(size.title = size.title, size.subtitle = size.subtitle, size.xlabel = size.xlabel, size.ylabel = size.ylabel,
                               size.axis = size.axis, size.legend = size.legend, size.legend.title = size.legend.title, legend.title = legend.title,
                               legend.title.pch = legend.title.pch, legend.position = legend.position, point.lwd = point.lwd)
        
        if (inherits(object, c("mint.block.pls", "mint.block.spls", "mint.block.plsda", "mint.block.splsda")))
            stop("No plotIndiv for the following functions at this stage: mint.block.pls, mint.block.spls, mint.block.plsda, mint.block.splsda.")
        
        #-- choose rep.space
        if (is.null(rep.space) && inherits(object, "DA"))#"splsda", "plsda", "mlsplsda")))
        {
            rep.space = "X-variate"
        } else if (is.null(rep.space)) {
            rep.space = "multi"
        }
        rep.space  = match.arg(rep.space, c("XY-variate", "X-variate", "Y-variate", "multi"))
        #c("XY-variate", "X-variate", "Y-variate", "multi")[pmatch(rep.space, c("XY-variate", "X-variate", "Y-variate", "multi"))]
        
        if (rep.space  == "multi")
        {
            blocks = c("X", "Y")
            object$variates  = object$variates[names(object$variates) %in% blocks]
        }
        
        if (rep.space  == "X-variate")
        {
            object$variates  = object$variates["X"]
            blocks  = "X"
        }
        
        if (rep.space  == "Y-variate")
        {
            object$variates  = object$variates["Y"]
            blocks  = "Y"
        }
        
        if (rep.space  == "XY-variate")
        {
            object$variates$XYvariates  = (object$variates$X + object$variates$Y)/2
            object$variates  = object$variates["XYvariates"]
            blocks  = "XY combined"
        }
        
        if (length(blocks)!= length(unique(blocks)))
            stop("Duplicate in 'blocks' not allowed")
        
        if (!missing(subtitle))
        {
            if (length(subtitle)!= length(blocks) | length(subtitle)!= length(unique(subtitle)))
                stop("'subtitle' indicates the subtitle of the plot for each 'blocks'; it needs to be the same length as 'blocks' and duplicate are not allowed.")
        }
        
        if(!is.null(background) &&  !is(background, "background.predict"))
            stop("'background' must have been obtained with the 'background.predict' function")
        
        #-- check inputs
        check  = check.input.plotIndiv(object = object, comp  = comp , blocks  = blocks, ind.names  = ind.names, 
                                       style = style, ellipse  = ellipse, ellipse.level  = ellipse.level, centroid = centroid, 
                                       star = star, legend = legend, X.label  = X.label, Y.label  = Y.label, Z.label  = Z.label, abline  = abline, 
                                       xlim  = xlim, ylim  = ylim, alpha = alpha, axes.box  = axes.box, plot_parameters = plot_parameters)
        #-- retrieve outputs from the checks
        axes.box = check$axes.box
        comp = check$comp
        xlim = check$xlim
        ylim = check$ylim
        ind.names = check$ind.names
        display.names = check$display.names
        
        #-- get the variates
        variate = internal_getVariatesAndLabels(object, comp, blocks = blocks, rep.space = rep.space, style = style, X.label = X.label,
                                                Y.label = Y.label, Z.label = Z.label)
        #-- retrieve outputs
        x = variate$x
        y = variate$y
        z = variate$z
        X.label = variate$X.label
        Y.label = variate$Y.label
        Z.label = variate$Z.label
        
        n = nrow(object$X)
        
        # create data frame df that contains (almost) all the ploting information
        out = shape.input.plotIndiv(object = object, n = n, blocks  = blocks, x = x, y = y, z = z, ind.names  = ind.names, group = group,
                                    col.per.group = col.per.group, style = style, study = "global", ellipse  = ellipse, ellipse.level  = ellipse.level,
                                    centroid = centroid, star = star, title = title, xlim  = xlim, ylim  = ylim, 
                                    col = col, cex = cex, pch = pch, pch.levels = pch.levels, display.names = display.names, plot_parameters = plot_parameters)
        #-- retrieve outputs
        df = out$df
        df.ellipse = out$df.ellipse
        col.per.group = out$col.per.group
        title = out$title
        display.names = out$display.names
        xlim = out$xlim
        ylim = out$ylim
        #missing.col = out$missing.col
        ellipse = out$ellipse
        centroid = out$centroid
        star = out$star
        plot_parameters = out$plot_parameters
        
        # change the levels of df$Block to "subtitle"
        if (!missing(subtitle) & nlevels(df$Block)>1)#& !is.null(title)) # commented so that subtitle can be change without changing the title
        {
            df$Block = factor(df$Block, labels = subtitle)
            if (ellipse)
                df.ellipse$Block = factor(df.ellipse$Block, labels = subtitle)
        }
        
        # match background color to col.per.group, the color of the groups
        if(!is.null(background))
        {
            ind.match = match(names(background), levels(df$group))
            names(background) = adjustcolor(col.per.group[ind.match],alpha.f=0.1)
        }
        
        #save(list = ls(), file = "temp.Rdata")
        
        #call plot module (ggplot2, lattice, graphics, 3d)
        res = internal_graphicModule(df = df, centroid = centroid, col.per.group = col.per.group, title = title,
                                     X.label = X.label, Y.label = Y.label, Z.label = Z.label, xlim = xlim, ylim = ylim, class.object = class(object),
                                     display.names = display.names, legend = legend, abline = abline, star = star,
                                     ellipse = ellipse, df.ellipse = df.ellipse, style = style, layout = layout, #missing.col = missing.col,
                                     axes.box = axes.box, plot_parameters = plot_parameters, alpha = alpha, background = background)
        
        
        return(invisible(list(df = df, df.ellipse = df.ellipse, graph = res)))
        
    }

#' @method plotIndiv mixo_spls
#' @export
plotIndiv.mixo_spls <- plotIndiv.mixo_pls

#' @method plotIndiv rcc
#' @export
plotIndiv.rcc <- plotIndiv.mixo_pls