## --------------------------- plot.perf.(s)pls --------------------------- ##
#' Plot for model performance
#' 
#' Function to plot performance criteria, such as MSEP, RMSEP, \eqn{R^2},
#' \eqn{Q^2} for s/PLS methods, and classification performance for supervised
#' methods, as a function of the number of components.
#' 
#' \code{plot.perf} creates one plot for each response variable in the model,
#' laid out in a multi panel display.  It uses \code{\link{xyplot}} for
#' performing the actual plotting.
#' 
#' More details about the prediction distances in \code{?predict} and the
#' supplemental material of the mixOmics article (Rohart et al. 2017).
#' 
#' @param x an \code{perf} object.
#' @param criterion character string. What type of validation criterion to plot
#' for \code{pls} or \code{spls}. One of \code{"MSEP"}, \code{"RMSEP"},
#' \code{"R2"} or \code{"Q2"}. See \code{\link{perf}}.
#' @param dist prediction method applied in \code{perf} for \code{plsda} or
#' \code{splsda}. See \code{\link{perf}}.
#' @param measure Two misclassification measure are available: overall
#' misclassification error \code{overall} or the Balanced Error Rate \code{BER}
#' @param col character (or symbol) color to be used, possibly vector. One
#' color per distance \code{dist}.
#' @param weighted plot either the performance of the Majority vote or the
#' Weighted vote.
#' @param study Indicates which study-specific outputs to plot. A character
#' vector containing some levels of \code{object$study}, "all.partial" to plot
#' all studies or "global" is expected. Default to "global".
#' @param overlay parameter to overlay graphs; if 'all', only one graph is
#' shown with all outputs; if 'measure', a graph is shown per distance; if
#' 'dist', a graph is shown per measure.
#' @param legend.position position of the legend, one of "vertical" (only one
#' column) or "horizontal" (two columns).
#' @param xlab,ylab titles for \eqn{x} and \eqn{y} axes.  Typically character
#' strings, but can be expressions (e.g., \code{expression(R^2)}).
#' @param LimQ2 numeric value. Signification limit for the components in the
#' model. Default is \code{LimQ2 = 0.0975}.
#' @param LimQ2.col character string specifying the color for the \code{LimQ2}
#' line to be plotted. If \code{"none"} the line will not be plotted.
#' @param cTicks integer vector. Axis tickmark locations for the used number of
#' components. Default is \code{1:ncomp} (see \code{\link{perf}}).
#' @param layout numeric vector of length two giving the number of rows and
#' columns in a multi panel display. If not specified, \code{plot.perf} tries
#' to be intelligent.
#' @param sd If 'nrepeat' was used in the call to 'perf', error bar shows the
#' standard deviation if sd=TRUE
#' @param \dots Further arguments sent to \code{\link{xyplot}} function.
#' @return none
#' @author Ignacio González, Florian Rohart, Francois Bartolo, Kim-Anh Lê Cao, Al J Abadi
#' @seealso \code{\link{pls}}, \code{\link{spls}}, \code{\link{plsda}},
#' \code{\link{splsda}}, \code{\link{perf}}.
#' @references
#' 
#' Rohart F, Gautier B, Singh A, Lê Cao K-A. mixOmics: an R package for 'omics
#' feature selection and multiple data integration. PLoS Comput Biol 13(11):
#' e1005752
#' @keywords regression multivariate hplot
#' @name plot.perf
#' @example ./examples/plot.perf-examples.R
NULL
## --------------------------- plot.perf.(s)pls --------------------------- ##
#' @method plot perf.pls.mthd
#' @rdname plot.perf
#' @export
plot.perf.pls.mthd <-
    function (x,
              criterion = "MSEP",
              #c("MSEP", "RMSEP", "R2", "Q2"),
              xlab = "number of components",
              ylab = NULL,
              LimQ2 = 0.0975,
              LimQ2.col = "darkgrey",
              cTicks = NULL,
              layout = NULL,
              ...
    )
    {
        
        
        if (!any(criterion %in% c("MSEP", "RMSEP", "R2", "Q2")) || length(criterion) > 1)
            stop("Choose one validation criterion among MSEP, RMSEP, R2 or Q2.")
        
        y = switch(criterion, MSEP = x$MSEP, RMSEP = sqrt(x$MSEP), R2 = x$R2, Q2 = x$Q2)
        
        Q2.total = NULL
        if ((criterion == "Q2") & is.list(y)) {
            Q2.total = y$Q2.total
            y = y$variables
        }
        
        if (is.null(ylab))
            ylab = switch(criterion, MSEP = "MSEP", RMSEP = "RMSEP",
                          R2 = expression(R^~2), Q2 = expression(Q^~2))
        
        nResp = nrow(y)  # Number of response variables
        nComp = ncol(y)  # Number of components
        
        #def.par = par(no.readonly = TRUE)
        
        if (nResp > 1) {
            if (is.null(layout)) {
                nRows = min(c(3, nResp))
                nCols = min(c(3, ceiling(nResp / nRows)))
                layout = c(nRows, nCols)
            }
            else {
                if (length(layout) != 2 || !is.numeric(layout) || any(is.na(layout)))
                    stop("'layout' must be a numeric vector of length 2.")
                nRows = layout[1]
                nCols = layout[2]
            }
            
            if (nRows * nCols < nResp) devAskNewPage(TRUE)
            ynames = rownames(y)
        } else {
            ynames = "Y"
        }
        
        val = comps = vector("numeric")
        varName = vector("character")
        
        for (i in 1:nResp) {
            val = c(val, y[i, ])
            comps = c(comps, 1:nComp)
            varName = c(varName, rep(ynames[i], nComp))
        }
        
        df = data.frame(val = val, comps = comps, varName = varName)
        if (is.null(cTicks)) cTicks = 1:ncol(y)
        yList = list(relation = "free")
        
        
        if (criterion == "Q2")
        {
            plt = xyplot(val ~ comps | varName, data = df, xlab = xlab, ylab = ylab,
                         scales = list(y = yList, x = list(at = cTicks)),
                         as.table = TRUE, layout = layout,
                         panel = function(x, y) {
                             if (LimQ2.col != "none") panel.abline(h = LimQ2, col = LimQ2.col)
                             panel.xyplot(x, y, ...)})
            plot(plt)
            
            if (!is.null(Q2.total)) {
                devAskNewPage(TRUE)
                Q2.df = data.frame(Q2 = Q2.total, comps = 1:nComp, varName = rep("Total", nComp))
                xyplot(Q2 ~ comps | varName, data = Q2.df, xlab = xlab, ylab = ylab,
                       scales = list(y = yList, x = list(at = cTicks)), as.table = TRUE,
                       panel = function(x, y) {
                           if (LimQ2.col != "none") panel.abline(h = LimQ2, col = LimQ2.col)
                           panel.xyplot(x, y, ...)})
            }
        } else {
            plt = xyplot(val ~ comps | varName, data = df, xlab = xlab, ylab = ylab,
                         scales = list(y = yList, x = list(at = cTicks)),
                         as.table = TRUE, layout = layout, ...)
            plot(plt)
            
        }
        
        if (nResp > 1) {
            if (nRows * nCols < nResp) devAskNewPage(FALSE)
        }
        
        #par(def.par)
        
        
    }

#' @method plot perf.spls.mthd
#' @rdname plot.perf
#' @export
plot.perf.spls.mthd <- plot.perf.pls.mthd

## -------------------------- plot.perf.(s)plsda -------------------------- ##
#' @method plot perf.plsda.mthd
#' @rdname plot.perf
#' @importFrom methods hasArg
#' @export
plot.perf.plsda.mthd <- 
    function (x,
              dist = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"),
              measure = c("all", "overall", "BER"),
              col,
              xlab = NULL,
              ylab = NULL,
              overlay = c("all", "measure", "dist"),
              legend.position = c("vertical", "horizontal"),
              sd = TRUE,
              ...
    )
    {
        # maybe later, so far we set type = "l"
        type = "l"
        
        if (hasArg(pred.method))
            stop("'pred.method' argument has been replaced by 'dist' to match the 'tune' and 'perf' functions")
        pred.method = NULL # to pass R CMD check
        
        
        if (any(measure == "all"))
            measure = names(x$error.rate)
        
        if (is.null(measure) || !any(measure %in% names(x$error.rate)))
            stop("'measure' should be among the ones used in your call to 'perf': ", paste(names(x$error.rate),collapse = ", "),".")
        
        if (any(dist == "all"))
            dist = colnames(x$error.rate[[1]])
        
        
        if (is.null(dist) || !any(dist %in% colnames(x$error.rate[[1]])))
            stop("'dist' should be among the ones used in your call to 'perf': ", paste(colnames(x$error.rate[[1]]),collapse = ", "),".")
        
        if(missing(col)) #one col per distance
        {
            col = color.mixo(1:length(dist))
        } else {
            if(length(col) != length(dist))
                stop("'col' should be a vector of length ", length(dist),".")
        }
        
        if (is.null(ylab))
            ylab = 'Classification error rate'
        
        if (is.null(xlab))
            xlab = 'Component'
        
        if(length(overlay) >1 )
            overlay = overlay[1]
        
        if(length(legend.position) >1 )
            legend.position = legend.position[1]
        
        # error.rate is a list [[measure]]
        # error.rate[[measure]] is a matrix of dist columns and ncomp rows
        # same for error.rate.sd, if any
        error.rate = x$error.rate
        if(sd)
        {
            error.rate.sd = x$error.rate.sd
        } else {
            error.rate.sd = NULL
        }
        def.par = par(no.readonly = TRUE)
        
        internal_graphic.perf(error.rate = error.rate, error.rate.sd = error.rate.sd,
                              overlay = overlay, type = type, measure = measure, dist = dist, legend.position = legend.position,
                              xlab = xlab, ylab = ylab, sd = sd, color = col, ...)
        
        par(def.par)
        # error.bar(out,as.vector(mat.error.plsda),as.vector(cbind(x$error.rate.sd$overall,x$error.rate.sd$BER)))
        
        return(invisible())
        
    }


#' @method plot perf.splsda.mthd
#' @rdname plot.perf
#' @export
plot.perf.splsda.mthd <- plot.perf.plsda.mthd

## ----------------------- plot.perf.mint.(s)plsda ------------------------ ##
#' @rdname plot.perf
#' @method plot perf.mint.plsda.mthd
#' @importFrom methods hasArg
#' @export
plot.perf.mint.plsda.mthd <-
    function (x,
              dist = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"),
              measure = c("all", "overall", "BER"),
              col,
              xlab = NULL,
              ylab = NULL,
              study = "global",
              overlay = c("all", "measure", "dist"),
              legend.position = c("vertical", "horizontal"),
              ...
    )
    {
        # maybe later, so far we set type = "l"
        type = "l"
        
        if (hasArg(pred.method))
            stop("'pred.method' argument has been replaced by 'dist' to match the 'tune' and 'perf' functions")
        pred.method = NULL # to pass R CMD check
        
        
        if (any(measure == "all"))
            measure = c("BER","overall")
        
        if (is.null(measure) || !any(measure %in% c("BER","overall")))
            stop("'measure' should be among the ones used in your call to 'perf': ", paste(c("BER","overall"),collapse = ", "),".")
        
        
        if (any(dist == "all"))
            dist = colnames(x$global.error[[1]])
        
        if(length(overlay) >1 )
            overlay = overlay[1]
        
        if(length(legend.position) >1 )
            legend.position = legend.position[1]
        
        if(missing(col)) #one col per distance
        {
            col = color.mixo(1:length(dist))
        } else {
            if(length(col) != length(dist))
                stop("'col' should be a vector of length ", length(dist),".")
        }
        
        
        if(any(study == "global"))
        {
            if (is.null(dist) || !any(dist %in% colnames(x$global.error[[1]])))
                stop("'dist' should be among the ones used in your call to 'perf': ", paste(colnames(x$global.error[[1]]),collapse = ", "),".")
            
            
            if (is.null(ylab))
                ylab = 'Classification error rate'
            
            if (is.null(xlab))
                xlab = 'Component'
            
            # error.rate is a list [[measure]]
            # error.rate[[measure]] is a matrix of dist columns and ncomp rows
            # same for error.rate.sd, if any
            
            error.rate = x$global.error
            
            def.par = par(no.readonly = TRUE)
            
            internal_graphic.perf(error.rate = error.rate, error.rate.sd = NULL,
                                  overlay = overlay, type = type, measure = measure, dist = dist, legend.position = legend.position,
                                  xlab = xlab, ylab = ylab, color = col, ...)
            
            par(def.par)
            
        } else {
            
            def.par = par(no.readonly = TRUE)
            
            
            if (any(study == "all.partial"))
                study = 1:length(x$study.specific.error)
            
            
            if (any(dist == "all"))
                dist = colnames(x$study.specific.error[[1]][[1]])
            
            if((length(study) >1) & (overlay != "all"))
                stop("When more than one study is plotted, overlay must be 'all'")
            
            
            if (is.null(dist) || !any(dist %in% colnames(x$study.specific.error[[1]][[1]])))
                stop("'dist' should be among the ones used in your call to 'perf': ", paste(colnames(x$study.specific.error[[1]][[1]]),collapse = ", "),".")
            
            if (is.null(ylab))
                ylab = 'Classification error rate'
            
            if (is.null(xlab))
                xlab = 'Component'
            
            
            if(overlay=="all")
            {
                par(mfrow=c(1,length(study)))
                
            } else if(overlay=="measure") {
                par(mfrow=c(length(study),length(dist)))
            } else if(overlay=="dist") {
                par(mfrow=c(length(study),length(measure)))
            }
            
            
            for(stu in study)
            {
                error.rate = x$study.specific.error[[stu]]
                
                internal_graphic.perf(error.rate = error.rate, error.rate.sd = NULL,
                                      overlay = overlay, type = type, measure = measure, dist = dist, legend.position = legend.position,
                                      xlab = xlab, ylab = ylab, color = col, ...)
                
                if (overlay == "all")
                    title(stu, line = 1)
            }
            
            if((length(study)==1) & (length(measure) > 1) & overlay != "all")
                title(stu, outer=TRUE, line = -1)#,...)
            
            
            par(def.par)
            
        }
        return(invisible())
        
    }

#' @method plot perf.mint.splsda.mthd
#' @rdname plot.perf
#' @export
plot.perf.mint.splsda.mthd <- plot.perf.mint.plsda.mthd

## --------------------------- plot.perf.sgccda --------------------------- ##
#' @method plot perf.sgccda.mthd 
#' @importFrom methods hasArg
#' @rdname plot.perf
#' @export
plot.perf.sgccda.mthd <-
    function (x,
              dist = c("all","max.dist","centroids.dist","mahalanobis.dist"),
              measure = c("all","overall","BER"),
              col,
              weighted = TRUE,
              xlab = NULL,
              ylab = NULL,
              overlay= c("all", "measure", "dist"),
              legend.position=c("vertical","horizontal"),
              sd = TRUE,
              ...)
    {
        # maybe later, so far we set type = "l"
        type = "l"
        
        if (hasArg(pred.method))
            stop("'pred.method' argument has been replaced by 'dist' to match the 'tune' and 'perf' functions")
        pred.method = NULL # to pass R CMD check
        
        
        measure.input = measure
        measure = NULL
        if(any(measure.input == "all"))
            measure.input = c("BER", "overall")
        
        if(any(measure.input == "BER"))
            measure = c(measure, "Overall.BER")
        
        if (any(measure.input == "overall"))
            measure = c(measure, "Overall.ER")
        
        if(!all(measure.input %in% c("all", "overall", "BER")))
            stop("'measure' must be 'all', 'overall' or 'BER'")
        
        if (any(dist == "all"))
            dist = colnames(x$error.rate[[1]])
        
        if(length(overlay) >1 )
            overlay = overlay[1]
        
        if(length(legend.position) >1 )
            legend.position = legend.position[1]
        
        
        if (is.null(dist) || !any(dist %in% colnames(x$error.rate[[1]])))
            stop("'dist' should be among the ones used in your call to 'perf': ", paste(colnames(x$error.rate[[1]]),collapse = ", "),".")
        
        if(missing(col)) #one col per distance
        {
            col = color.mixo(1:length(dist))
        } else {
            if(length(col) != length(dist))
                stop("'col' should be a vector of length ", length(dist),".")
        }
        
        if (is.null(ylab))
            ylab = 'Classification error rate'
        
        if (is.null(xlab))
            xlab = 'Component'
        
        if(weighted == TRUE)
        {
            perfo = "WeightedVote.error.rate"
            perfo.sd = "WeightedVote.error.rate.sd"
        } else {
            perfo = "MajorityVote.error.rate"
            perfo.sd = "MajorityVote.error.rate.sd"
        }
        
        if(sd == TRUE)
        {
            if(is.null(x[[perfo.sd]]))
                sd = FALSE
        }
        
        # error.rate is a list [[measure]]
        # error.rate[[measure]] is a matrix of dist columns and ncomp rows
        # same for error.rate.sd, if any
        error.rate = error.rate.sd = list()
        for(mea in measure)
        {
            error.temp = error.temp.sd = NULL
            for(di in dist)
            {
                temp = t(x[[perfo]][[di]][mea, , drop=FALSE])
                colnames(temp) = di
                error.temp = cbind(error.temp, temp)
                if(sd)
                {
                    temp.sd = t(x[[perfo.sd]][[di]][mea, , drop=FALSE])
                    colnames(temp.sd) = di
                    error.temp.sd = cbind(error.temp.sd, temp.sd)
                }
                
            }
            error.rate[[mea]] = error.temp
            if(sd)
            {
                error.rate.sd[[mea]] = error.temp.sd
            } else {
                error.rate.sd = NULL
            }
        }
        
        
        
        def.par = par(no.readonly = TRUE)
        
        internal_graphic.perf(error.rate = error.rate, error.rate.sd = error.rate.sd,
                              overlay = overlay, type = type, measure = measure, dist = dist, legend.position = legend.position,
                              xlab = xlab, ylab = ylab, color = col, ...)
        
        par(def.par)
        return(invisible())
        
    }


