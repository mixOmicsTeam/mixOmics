#############################################################################################################
# Authors:
#   Sebastien Dejean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
#   Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
#   Kim-Anh Le Cao, French National Institute for Agricultural Research and
#   Queensland Facility for Advanced Bioinformatics, University of Queensland, Australia
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 2011
# last modified: 19-04-2016
#
# Copyright (C) 2011
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
#############################################################################################################

#############################################################################################################
# Authors:
#   Sebastien Dejean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
#   Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
#   Kim-Anh Le Cao, French National Institute for Agricultural Research and
#   Queensland Facility for Advanced Bioinformatics, University of Queensland, Australia
# Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 2011
# last modified: 19-04-2016
#
# Copyright (C) 2011
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
#############################################################################################################



# PLS object
# ----------------------

plot.perf.pls.mthd = plot.perf.spls.mthd =
function (x,
criterion = "MSEP",#c("MSEP", "RMSEP", "R2", "Q2"),
xlab = "number of components",
ylab = NULL,
LimQ2 = 0.0975,
LimQ2.col = "darkgrey",
cTicks = NULL,
layout = NULL,
...)
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

# PLSDA object
# ----------------------

plot.perf.plsda.mthd = plot.perf.splsda.mthd =
function (x,
dist = c("all","max.dist","centroids.dist","mahalanobis.dist"),
measure = c("all","overall","BER"),
col,
xlab = NULL,
ylab = NULL,
overlay=c("all", "measure", "dist"),
legend.position=c("vertical", "horizontal"),
sd = TRUE,
...)
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

# mint.PLSDA object
# ----------------------

plot.perf.mint.plsda.mthd = plot.perf.mint.splsda.mthd =
function (x,
dist = c("all","max.dist","centroids.dist","mahalanobis.dist"),
measure = c("all","overall","BER"),
col,
xlab = NULL,
ylab = NULL,
study = "global",
overlay= c("all", "measure", "dist"),
legend.position=c("vertical", "horizontal"),
...)
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

# SGCCDA object
# ----------------------

plot.perf.sgccda.mthd =
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

