## ----------------------------- MINT.(s)PLS ------------------------------ ##

#' @rdname plotLoadings
#' @method plotLoadings mint.pls
#' @export
plotLoadings.mint.pls <-
    function(object,
             style = "graphics",
             comp = 1,
             ndisplay = NULL,
             xlim = NULL,
             layout = NULL,
             col = NULL,
             border = NA,
             name.var = NULL,
             size.name = 0.7,
             title = NULL,
             subtitle,
             size.title = rel(1.8),
             size.subtitle = rel(1.4),
             size.axis = 0.7,
             X.label = NULL,
             Y.label = NULL,
             size.labs = 1,
             block,
             study = "global",
             ...
    ) {
      
      ## Check input args
      # Input checks
      if (!is.numeric(comp) || length(comp) != 1 || comp <= 0)
        stop("'comp' must be a positive integer.")
      if (!style %in% c('graphics', 'ggplot2'))
        stop("'style' must be either 'graphics' or 'ggplot2'.")
      
      if (!is.null(col) && (length(col) != 1 || !col %in% colors()))
        stop("'col' must be a single valid color.")
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
      name.var.complete <- FALSE
        
        if(any(study == "global"))
        {
            # if study == "global" then we plot the results on the concatenated data, thus direct call to plotLoadings.plsda
            plotLoadings.mixo_pls(object = object, 
                                  comp = comp, 
                                  ndisplay = ndisplay,
                                  size.name = size.name,
                                  name.var = name.var,
                                  title = title,
                                  subtitle = subtitle,
                                  layout = layout,
                                  size.title = size.title,
                                  size.subtitle = size.subtitle,
                                  border = border,
                                  xlim = xlim,
                                  col = col,
                                  X.label = X.label,
                                  Y.label = Y.label,
                                  size.labs = size.labs,
                                  size.axis = size.axis,
                                  block = block)
            
        } else {
            # if study != "global" then we plot the results on each study
            
            # -- input checks
            check = check.input.plotLoadings(object = object, block = object$names$blocks, title = title, col = col, size.name = size.name, name.var = name.var)
            
            col = check$col
            size.name = check$size.name
            # block = check$block # uses actual block names from object
            
            #study needs to be either: from levels(object$study), "all.partial" or "global"
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
            # study = study.init
            
            if (!missing(subtitle))
            {
                if (length(subtitle)!=length(study))
                    stop("'subtitle' indicates the subtitle of the plot for each study and it needs to be the same length as 'study' (", length(study),"), which includes: ", paste(study, collapse = ", "))
            }
            
            # Handle block selection for all.partial case
            if (any(study == "all.partial") | any(study %in% levels(object$study))) {
                if (missing(block)) {
                    block = object$names$blocks[1]  # default to first block if not specified
                } else if (length(block) > 1) {
                    stop("When study = 'all.partial' or a specific study is specified, only one block can be plotted at a time. Please specify one of: ", 
                         paste(object$names$blocks, collapse = ", "))
                }
            }
            
            # # swap block for study
            # block = study
            
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
            res = layout.plotLoadings(layout = layout, plot = TRUE, legend = FALSE, block = study.init)
            reset.mfrow = res$reset.mfrow
            opar = res$opar
            omar = par("mar") #reset mar at the end
            
            # get the selected variables on the concatenated data
            res = get.loadings.ndisplay(object = object, comp = comp, block = block, name.var = name.var, name.var.complete = name.var.complete, ndisplay = ndisplay)
            X = res$X
            colnames.X = res$colnames.X
            name.selected.var = res$name.selected.var
            value.selected.var = res$value.selected.var
            
            
            # swap loadings partial for loadings
            object$loadings.global = object$loadings
            object$loadings = object$loadings.partial[[block]]
            object$names$block = levels(object$study)
            
            df.final = list()
            for (i in 1 : length(study.init))
            {
                value.selected.var = object$loadings.partial[[block]][[study.init[i]]][, comp] [name.selected.var]
                
                df = data.frame(importance = value.selected.var, color = col, stringsAsFactors = FALSE) # contribution of the loading
                
                #display barplot with names of variables
                #added condition if all we need is the contribution stats
                if (!is.null(title) & length(block) > 1)
                {
                    par(mar = c(4, max(7, max(sapply(colnames.X, nchar),na.rm = TRUE)/2), 6, 2))
                } else {
                    par(mar = c(4, max(7, max(sapply(colnames.X, nchar),na.rm = TRUE)/2), 4, 2))
                }
                .plotLoadings_barplot(height = df$importance, col = df$color, names.arg = colnames.X, 
                        cex.name = size.name, border = border, xlim = xlim,
                        xlab = X.label, ylab = Y.label, cex.lab = size.labs, cex.axis = size.axis)
                
                if ( length(study.init) == 1 & is.null(title) )
                {
                    title(paste0('Loadings on comp ', comp), line=1, cex.main = size.title)
                } else if (length(study.init) == 1) {
                    title(paste(title), line=0, cex.main = size.title)
                } else if ((length(study.init) > 1 & missing(subtitle))) {
                    title(paste0('Loadings on comp ', comp, "\nStudy '", study.init[i],"'"), line=0, cex.main = size.subtitle)
                } else if (length(study.init) > 1 & !missing(subtitle)) {
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