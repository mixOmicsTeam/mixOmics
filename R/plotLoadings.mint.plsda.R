## ---------------------------- MINT.(s)PLSDA ----------------------------- ##

#' @rdname plotLoadings
#' @method plotLoadings mint.plsda
#' @export
plotLoadings.mint.plsda <- 
    function(object,
             comp = 1,
             style = "graphics",
             ndisplay = NULL,
             xlim = NULL,
             layout = NULL,
             border = NA,
             name.var = NULL,
             size.name = 0.7,
             title = NULL,
             subtitle,
             size.title = 2,
             size.subtitle = 1.7,
             size.axis = 0.7,
             X.label = NULL,
             Y.label = NULL,
             size.labs = 1,
             contrib = NULL,  # choose between 'max" or "min", NULL does not color the barplot
             method = "mean", # choose between 'mean" or "median"
             show.ties = TRUE,
             col.ties = "white",
             legend = TRUE,
             legend.color = NULL,
             legend.title = 'Outcome',
             size.legend = 0.8,
             study = "global",
             ...
    ) {

        ## Check input args
      # Input checks
      if (!is.numeric(comp) || length(comp) != 1 || comp <= 0)
        stop("'comp' must be a positive integer.")
      if (!style %in% c('graphics', 'ggplot2'))
        stop("'style' must be either 'graphics' or 'ggplot2'.")
      if (!is.null(ndisplay) && (!is.numeric(ndisplay) || length(ndisplay) != 1 || ndisplay <= 0))
        stop("'ndisplay' must be a positive integer.")
      
      if (!is.null(title) && !is.character(title))
        stop("'title' must be NULL or a character string.")
      if (!is.null(xlim) && (!is.numeric(xlim) || length(xlim) != 2))
        stop("'xlim' must be a numeric vector of length 2.")
      if (!is.null(X.label) && !is.character(X.label))
        stop("'X.label' must be NULL or a character string.")
      if (!is.null(Y.label) && !is.character(Y.label))
        stop("'Y.label' must be NULL or a character string.")
      
      if (!is.numeric(size.name) || size.name <= 0)
        stop("'size.name' must be a positive numeric value.")
      if (!is.numeric(size.title) || size.title <= 0)
        stop("'size.title' must be a positive numeric value.")
      if (!is.numeric(size.subtitle) || size.subtitle <= 0)
        stop("'size.subtitle' must be a positive numeric value.")
      if (!is.numeric(size.labs) || size.labs <= 0)
        stop("'size.labs' must be a positive numeric value.")
      if (!is.numeric(size.axis) || size.axis <= 0)
        stop("'size.axis' must be a positive numeric value.")
      
      # check for inappropriate args
      extra_args <- list(...)
      if ("name.var.complete" %in% names(extra_args)) {
        warning("'name.var.complete' argument is deprecated")
      }

      if ("plot" %in% names(extra_args)) {
        warning("'plot' argument is deprecated")
      }
        
        if(any(study == "global"))
        {
            plotLoadings.mixo_plsda(object = object, 
                                    style = style,
                                    contrib = contrib, 
                                    method = method, 
                                    block = "X", 
                                    comp = comp, 
                                    ndisplay = ndisplay,
                                    size.name = size.name,
                                    size.legend = size.legend,
                                    name.var = name.var,
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
            object$loadings = object$loadings.partial[[block]]
            object$names$block = levels(object$study)
            
            X.study = study_split(X, study = object$study)
            Y = object$Y #v6: all $Y are factors for DA methods
            Y.study = study_split(Y, study = object$study)
            
            df.final = list()
            for (i in 1 : length(block))
            {
                
                value.selected.var =  object$loadings.partial[[block]][[block[i]]][, comp] [name.selected.var]
                
                
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
                
                .plotLoadings_barplot(height = df$importance, 
                                    col = df$color, 
                                    names.arg = colnames.X, 
                                    cex.name = size.name, 
                                    border = border, 
                                    xlim = xlim[i, ],
                                    xlab = X.label, ylab = Y.label, cex.lab = size.labs, cex.axis = size.axis)
                
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
