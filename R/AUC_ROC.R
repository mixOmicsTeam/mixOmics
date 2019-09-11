###############################################################################
#Authors:
#    Francois Bartolo,
#    Benoit Gautier,
#    Florian Rohart,
#    Kim-Anh Le Cao
#
# created: 23-08-2016
# last modified: 23-08-2016
#
# Copyright (C) 2016
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

auroc = function(object, ...)
UseMethod("auroc")


# PLSDA object
# ----------------------
auroc.mixo_plsda = auroc.mixo_splsda = function(
object,
newdata = object$input.X,
outcome.test = as.factor(object$Y),
multilevel = NULL,
plot = TRUE,
roc.comp = 1,
title=paste("ROC Curve Comp",roc.comp),
...)
{
    if(dim(newdata)[[1]] != length(outcome.test))
    stop("Factor outcome.test must be a factor with ",dim(newdata)[[1]],
    " elements.",call. = FALSE)
    
    data = list()
    statauc = graph = list()
    data$outcome=factor(outcome.test)
    
    # note here: the dist does not matter as we used the predicted scores only
    res.predict = predict.mixo_spls(object, newdata = newdata,
    dist = "max.dist", multilevel = multilevel)$predict
    
    for (i in seq_len(object$ncomp))
    {
        data$data=res.predict[,,i]
        temp = statauc(data, plot = ifelse(i%in%roc.comp,plot,FALSE), title=title,...)
        statauc[[paste0("Comp", i, sep = "")]] = temp[[1]]
        graph[[paste0("Comp", i, sep = "")]] = temp$graph
    }
    print(statauc)
    return(invisible(c(statauc,graph=graph)))
}


# MINT object
# ----------------------
auroc.mint.plsda = auroc.mint.splsda = function(
object,
newdata = object$X,
outcome.test = as.factor(object$Y),
study.test = object$study,
multilevel = NULL,
plot = TRUE,
roc.comp = 1,
roc.study = "global",
title=NULL,
...)
{
    if(length(roc.study) != 1)
    stop("`roc.study' must be a single entry,
    either `global' or one of levels(object$study)")
    
    if(roc.study == "global"){
        if(dim(newdata)[[1]] != length(outcome.test))
        stop("Factor outcome.test must be a factor with ",dim(newdata)[[1]],
        " elements.",call. = FALSE)
        
        if(dim(newdata)[[1]]!=length(study.test))
        stop("Factor study.test must be a factor with ",dim(newdata)[[1]],
        " elements.",call. = FALSE)
        study.test=factor(study.test)
        title.temp = NULL
        
    } else {
        
        # check study
        if (!roc.study%in%c(levels(object$study)))
        stop("'roc.study' must be one of 'levels(object$study)'")
        
        ind.study = object$study == roc.study
        newdata = object$X[ind.study, ]
        outcome.test = as.factor(object$Y[ind.study])
        study.test = factor(object$study[ind.study])
        title.temp = paste0(", Study ", roc.study)
        
    }
    
    data=list()
    statauc = graph = list()
    data$outcome=factor(outcome.test)
    
    # note here: the dist does not matter as we used the predicted scores only
    res.predict = predict.mixo_spls(object, newdata = newdata, dist = "max.dist",
    multilevel = multilevel, study.test = study.test)$predict
    
    for (i in seq_len(object$ncomp))
    {
        data$data=res.predict[,,i]
        
        if (is.null(title)) {
            title=paste0("ROC Curve Comp ",i, title.temp)
        }
        
        temp = statauc(data, plot = ifelse(i%in%roc.comp,plot,FALSE), title=title,...)
        statauc[[paste0("Comp", i, sep = "")]] = temp[[1]]
        graph[[paste0("Comp", i, sep = "")]] = temp$graph
    }
    print(statauc)
    return(invisible(c(statauc,graph=graph)))
    
}


# block.splsda object
# ----------------------
auroc.sgccda = function(
object,
newdata = object$X,
outcome.test = as.factor(object$Y),
multilevel = NULL,
plot = TRUE,
roc.block = 1L,
roc.comp = 1L,
title=NULL,
...)
{
    
    data=list()
    auc.mean = graph=list()
    data$outcome=factor(outcome.test)
    
    # note here: the dist does not matter as we used the predicted scores only
    res.predict  =  predict.block.spls(object, newdata = newdata,
        dist = "max.dist", multilevel = multilevel)$predict
    block.all = names(res.predict)
    if (is(roc.block, "numeric")) {
        roc.block <- as.integer(roc.block)
        lb <- length(names(res.predict))
        if (roc.block > lb)
            stop(sprintf("roc.block cannot be greater than %s", lb ))
        block.temp = names(res.predict[roc.block])
    } else if (is(roc.block, "character")) {
        block.temp = roc.block 
    } else {
        stop("'roc.block' should be an integer or character")
    }
    title.temp = title
    for(j in seq_len(length(res.predict)))
    {
        for (i in seq_len(object$ncomp[j]))
        {
            data$data=res.predict[[j]][,,i]
            if (is.null(title.temp)) {
                title=paste("ROC Curve\nBlock: ", names(res.predict)[j],
                            ", comp: ",i, sep="")
            }
            
            plot.temp =
                ifelse(i%in%roc.comp && names(res.predict)[j]%in%block.temp,
                plot, FALSE)
            temp = statauc(data, plot = plot.temp, title = title, ...)
            auc.mean[[names(res.predict)[j]]][[paste0("comp",i,sep = "")]] =
                temp[[1]]
            graph[[names(res.predict)[j]]][[paste0("comp",i,sep = "")]] =
                temp$graph
            
        }
        out = c(auc.mean,graph=graph)
    }
    print(auc.mean)
    return(invisible(out))
}

# mint.block.splsda object
# ----------------------
auroc.mint.block.splsda=auroc.mint.block.plsda = function(
object,
newdata = object$X,
study.test = object$study,
outcome.test = as.factor(object$Y),
multilevel = NULL,
plot = TRUE,
roc.block = 1,
roc.comp = 1,
title=NULL,
...)
{
    
    data=list()
    auc.mean = graph=list()
    data$outcome=factor(outcome.test)
    study.test=factor(study.test)
    
    # note here: the dist does not matter as we used the predicted scores only
    res.predict  =  predict.mixo_spls(object, newdata = newdata,
    study.test=study.test,dist = "max.dist", multilevel = multilevel)$predict
    block.all = names(res.predict)
    block.temp = names(res.predict[roc.block])
    
    for(j in seq_len(length(res.predict)))
    {
        for (i in seq_len(object$ncomp[j]))
        {
            data$data=res.predict[[j]][,,i]
            if (is.null(title)) {
                title=paste("ROC Curve\nBlock: ", names(res.predict)[j],
                            ", comp: ",i, sep="")
            }

            plot.temp =
                ifelse(i%in%roc.comp && names(res.predict)[j]%in%block.temp,
                plot, FALSE)
            temp = statauc(data, plot = plot.temp, title = title, ...)
            auc.mean[[names(res.predict)[j]]][[paste0("comp",i,sep = "")]] =
                temp[[1]]
            graph[[names(res.predict)[j]]][[paste0("comp",i,sep = "")]] =
                temp$graph
            
        }
        out = c(auc.mean,graph=graph)
    }
    print(auc.mean)
    return(invisible(out))
}

