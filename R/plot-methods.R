## ------------------------------- plot.pca ------------------------------- ##
#' Show PCA output
#' 
#' @param x A \code{pca} object
#' @param ncomp Integer, the number of components
#' @param type Character, default "barplot" or any other type available in plot, as "l","b","p",..
#' @param explained.var Logical, whether to plot the explained variance
#' @param ... Not used
#' @author Kim-Anh Lê Cao, Florian Rohart, Leigh Coonan, Al J Abadi
#' @method plot pca
#' @export
plot.pca  <- function(x,
                      ncomp = min(10, length(x$sdev)),
                      type = "barplot",
                      # either barplot or any other type available in plot, as "l","b","p",..
                      explained.var = TRUE,
                      ...)
{
    #-- checking general input parameters --------------------------------------#
    #---------------------------------------------------------------------------#
    
    #-- ncomp
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp < 1 || !is.finite(ncomp))
        stop("invalid value for 'ncomp'.", call. = FALSE)
    
    ncomp = round(ncomp)
    
    if (ncomp > length(x$sdev))
        stop("'ncomp' must be lower or equal than ", length(x$sdev), ".",
             call. = FALSE)
    
    #-- end checking --#
    #------------------#
    
    #-- scree plot -------------------------------------------------------------#
    #---------------------------------------------------------------------------#
    
    variances = (x$sdev^2)[1:ncomp] # relative variance
    ylab = "Variance"
    if(explained.var==TRUE)
    {
        variances=variances/x$var.tot #explained variances
        ylab = "Explained Variance"
    }
    if (type == "barplot")
    {
        barplot(variances, names.arg = seq(1, ncomp), xlab = "Principal Components", ylab = ylab,...)
    } else {
        plot(variances, type = type, axes = FALSE,
             xlab = "Principal Components",
             ylab = ylab,... )
        axis(1, at = 1:ncomp)
        axis(2)
    }
    
}

## ------------------------------ plot.spca ------------------------------- ##
#' Show sPCA output
#' 
#' @param x A \code{spca} object
#' @param ncomp Integer, the number of components
#' @param type Character, default "barplot" or any other type available in plot, as "l","b","p",..
#' @param ... Not used
#' @author Al J Abadi, Kim-Anh Lê Cao, Florian Rohart, Leigh Coonan
#' @method plot spca
#' @export
plot.spca  <- function(x,
                      ncomp = length(x$explained_variance),
                      type = "barplot",
                      # either barplot or any other type available in plot, as "l","b","p",..
                      ...)
{
    #-- checking general input parameters --------------------------------------#
    #---------------------------------------------------------------------------#
    
    #-- ncomp
    ncomp_model <- length(x$explained_variance)
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp < 1 || !is.finite(ncomp))
        stop("invalid value for 'ncomp'.", call. = FALSE)
    
    ncomp = round(ncomp)
    
    if (ncomp > ncomp_model)
        stop("'ncomp' must be lower or equal than ", length(ncomp_model), ".",
             call. = FALSE)
    
    #-- end checking --#
    #------------------#
    
    #-- scree plot -------------------------------------------------------------#
    #---------------------------------------------------------------------------#
    
    expl_vars = (x$explained_variance)[seq_len(ncomp)] # relative variance
    ylab = "Explained Variance"
    if (type == "barplot")
    {
        barplot(expl_vars, names.arg = seq_len(ncomp), xlab = "Principal Components", ylab = ylab,...)
    } else {
        plot(expl_vars, type = type, axes = FALSE,
             xlab = "Principal Components",
             ylab = ylab,... )
        axis(1, at = seq_len(ncomp))
        axis(2)
    }
    
}

## ------------------------------- plot.rcc ------------------------------- ##
#' Canonical Correlations Plot
#' 
#' This function provides scree plot of the canonical correlations.
#' 
#' 
#' @param x object of class inheriting from \code{"rcc"}.
#' @param scree.type character string, (partially) matching one of
#' \code{"pointplot"} or \code{"barplot"}, determining the kind of scree plots
#' to be produced.
#' @param ... arguments to be passed to other methods. For the
#' \code{"pointplot"} type see \code{\link{points}}, for \code{"barplot"} type
#' see \code{\link{barplot}}.
#' @return none
#' @author Sébastien Déjean, Ignacio González, Al J Abadi
#' @seealso \code{\link{points}}, \code{\link{barplot}}, \code{\link{par}}.
#' @keywords multivariate hplot
#' @method plot rcc
#' @export
#' @examples
#' 
#' data(nutrimouse)
#' X <- nutrimouse$lipid
#' Y <- nutrimouse$gene
#' nutri.res <- rcc(X, Y, lambda1 = 0.064, lambda2 = 0.008)
#' 
#' ## 'pointplot' type scree
#' plot(nutri.res) #(default)
#' 
#' \dontrun{
#' plot(nutri.res, pch = 19, cex = 1.2,
#' col = c(rep("red", 3), rep("darkblue", 18)))
#' 
#' ## 'barplot' type scree
#' plot(nutri.res, scree.type = "barplot")
#' 
#' plot(nutri.res, scree.type = "barplot", density = 20, col = "black")
#' }
#' 
plot.rcc <-
    function(x, scree.type = c("pointplot", "barplot"), ...) 
    {
        
        scree.type = match.arg(scree.type)
        if (scree.type == "pointplot") {
            plot(x$cor, xlab = "Dimension", ylim = c(0, 1),
                 ylab = "Canonical correlation", ...)
        }
        else {
            barplot(x$cor, xlab = "Dimension", ylim = c(0, 1),
                    ylab = "Canonical correlation", ...)
        }
    }

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
#' @method plot perf.pls.mthd
#' @export
#' @example ./examples/plot.perf-examples.R
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

#' @rdname plot.perf
#' @method plot perf.spls.mthd
#' @export
plot.perf.spls.mthd <- plot.perf.pls.mthd

## -------------------------- plot.perf.(s)plsda -------------------------- ##
#' @rdname plot.perf
#' @method plot perf.plsda.mthd
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


#' @rdname plot.perf
#' @method plot perf.splsda.mthd
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

#' @rdname plot.perf
#' @method plot perf.mint.splsda.mthd
#' @export
plot.perf.mint.splsda.mthd <- plot.perf.mint.plsda.mthd

## --------------------------- plot.perf.sgccda --------------------------- ##
#' @rdname plot.perf
#' @method plot perf.sgccda.mthd 
#' @importFrom methods hasArg
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

## --------------------------- plot.tune.(s)pls --------------------------- ##
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
#' combination of keepX at the top (e.g. `keepX on block 1'_`keepX on block
#' 2'_`keepX on block 3')
#' 
#' @param x an \code{tune.splsda} object.
#' @param optimal If TRUE, highlights the optimal keepX per component
#' @param sd If 'nrepeat' was used in the call to 'tune.splsda', error bar
#' shows the standard deviation if sd=TRUE
#' @param col character (or symbol) color to be used, possibly vector. One
#' color per component.
#' @param \dots Further arguments sent to \code{\link{xyplot}} function.
#' @return none
#' @author Kim-Anh Lê Cao, Florian Rohart, Francois Bartolo, AL J Abadi
#' @seealso \code{\link{tune.mint.splsda}}, \code{\link{tune.splsda}}
#' \code{\link{tune.block.splsda}} and http://www.mixOmics.org for more
#' details.
#' @keywords regression multivariate hplot
#' @name plot.tune
#' @method plot tune.spls
#' @importFrom reshape2 melt
#' @export
#' @example ./examples/plot.tune-examples.R
plot.tune.spls <-
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

## --------------------------- plot.tune.spca --------------------------- ##
#' Plot for tune.spca model performance
#' 
#' Function to plot the correlation of cross-validated components from
#' the \code{tune.spca} function with respect to the full model.
#' 
#' 
#' @param x an \code{tune.splsda} object.
#' @param optimal If TRUE, highlights the optimal keepX per component
#' @param sd If \code{nrepeat>2} was used in the call to \code{tune.spca}, error bar
#' shows the standard deviation if \code{sd=TRUE}.
#' @param col character (or symbol) colour to be used, possibly vector. One
#' colour per component.
#' @param \dots Further arguments sent to \code{\link{xyplot}} function.
#' @return None
#' @author Kim-Anh Lê Cao, Al J Abadi
#' @seealso \code{\link{spca}}, \code{\link{tune.pca}}
#' @keywords regression multivariate hplot
#' @name plot.tune
#' @method plot tune.spca
#' @export
plot.tune.spca <-
    function(x, optimal = TRUE, sd = TRUE, col=NULL, ...)
    {
        ncomp <- length(x$cor.comp)
        nrepeat <- x$call$nrepeat
        
        if (nrepeat <= 2 & isTRUE(sd))
        {
            cat("nrepeat < 2 so no SD can be calculated.",
                "setting sd to FALSE")
            sd <- FALSE
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
        if (any(cors$corQ1 < -1 | cors$corQ3 > 1))
        { #TODO truncated distribution assumption needed?
            message("The interquartile ranges for correlations ",
                    "contain values outside valid range [-1, 1]. Consider larger 'nrepeat'.")
        }
        
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
            scale_x_log10() +
            ylim(c(min(cors$corQ1, 0), max(cors$corQ3, 0))) +
            labs(x= 'Number of features selected',
                 y = 'Correlation of components',
                 col = 'Comp') +
            scale_color_manual(values = col)
        
        p <- p +  geom_point(data=cors[!is.na(cors$opt.keepX),], 
                             aes_string('keepX', 'cor.mean', col = 'comp'), 
                             size=6, shape = 18, show.legend = FALSE)
        
        ## ----- error bars
        if (isTRUE(sd))
        {
            p <- p + geom_errorbar(aes_string(ymin = 'corQ1', ymax = 'corQ3'), 
                                   position = position_dodge(0.02),
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
