#' Area Under the Curve (AUC) and Receiver Operating Characteristic (ROC)
#' curves for supervised classification
#' 
#' Calculates the AUC and plots ROC for supervised models from s/plsda,
#' mint.s/plsda and block.plsda, block.splsda or wrapper.sgccda functions.
#' 
#' For more than two classes in the categorical outcome Y, the AUC is
#' calculated as one class vs. the other and the ROC curves one class vs. the
#' others are output.
#' 
#' The ROC and AUC are calculated based on the predicted scores obtained from
#' the \code{predict} function applied to the multivariate methods
#' (\code{predict(object)$predict}). Our multivariate supervised methods
#' already use a prediction threshold based on distances (see \code{predict})
#' that optimally determine class membership of the samples tested. As such AUC
#' and ROC are not needed to estimate the performance of the model (see
#' \code{perf}, \code{tune} that report classification error rates). We provide
#' those outputs as complementary performance measures.
#' 
#' The pvalue is from a Wilcoxon test between the predicted scores between one
#' class vs the others.
#' 
#' External independent data set (\code{newdata}) and outcome
#' (\code{outcome.test}) can be input to calculate AUROC. The external data set
#' must have the same variables as the training data set (\code{object$X}).
#' 
#' If \code{newdata} is not provided, AUROC is calculated from the training
#' data set, and may result in overfitting (too optimistic results).
#' 
#' Note that for mint.plsda and mint.splsda objects, if \code{roc.study} is
#' different from "global", then \code{newdata}), \code{outcome.test} and
#' \code{sstudy.test} are not used.
#' 
#' @aliases auroc auroc.mixo_plsda auroc.mixo_splsda auroc.mint.plsda
#' auroc.mint.splsda auroc.sgccda
#' 
#' @param object Object of class inherited from one of the following supervised
#' analysis function: "plsda", "splsda", "mint.plsda", "mint.splsda",
#' "block.splsda" or "wrapper.sgccda"
#' @param newdata numeric matrix of predictors, by default set to the training
#' data set (see details).
#' @param outcome.test Either a factor or a class vector for the discrete
#' outcome, by default set to the outcome vector from the training set (see
#' details).
#' @param study.test For MINT objects, grouping factor indicating which samples
#' of \code{newdata} are from the same study. Overlap with \code{object$study}
#' are allowed.
#' @param multilevel Sample information when a newdata matrix is input and when
#' multilevel decomposition for repeated measurements is required. A numeric
#' matrix or data frame indicating the repeated measures on each individual,
#' i.e. the individuals ID. See examples in \code{splsda}.
#' @param plot Whether the ROC curves should be plotted, by default set to TRUE
#' (see details).
#' @param roc.comp Specify the component (integer) for which the ROC will be
#' plotted from the multivariate model, default to 1.
#' @param roc.block Specify the block number (integer) or the name of the block
#' (set of characters) for which the ROC will be plotted for a block.plsda or
#' block.splsda object, default to 1.
#' @param roc.study Specify the study for which the ROC will be plotted for a
#' mint.plsda or mint.splsda object, default to "global".
#' @param title Character, specifies the title of the plot.
#' @param print Logical, specifies whether the output should be printed.
#' @param ... external optional arguments for plotting - \code{line.col} for
#' custom colors and \code{legend.title} for custom legend title
#' @return Depending on the type of object used, a list that contains: The AUC
#' and Wilcoxon test pvalue for each 'one vs other' classes comparison
#' performed, either per component (splsda, plsda, mint.plsda, mint.splsda), or
#' per block and per component (wrapper.sgccda, block.plsda, blocksplsda).
#' @author Benoit Gautier, Francois Bartolo, Florian Rohart, Al J Abadi
#' @seealso \code{\link{tune}}, \code{\link{perf}}, and http://www.mixOmics.org
#' for more details.
#' @keywords regression multivariate
#' @example ./examples/auroc-examples.R
#' @export
auroc <- function(object, ...)
    UseMethod("auroc")


# PLSDA object
# ----------------------
#' @rdname auroc
#' @method auroc mixo_plsda
#' @export
auroc.mixo_plsda <- 
    function(
        object,
        newdata = object$input.X,
        outcome.test = as.factor(object$Y),
        multilevel = NULL,
        plot = TRUE,
        roc.comp = 1,
        title=paste("ROC Curve Comp",roc.comp),
        print=TRUE,
        ...)
    {
        if(dim(newdata)[[1]] != length(outcome.test))
            stop("Factor outcome.test must be a factor with ",dim(newdata)[[1]],
                 " elements.",call. = FALSE)
        
        data = list()
        statauc.res = graph = list()
        data$outcome=factor(outcome.test)
        
        # note here: the dist does not matter as we used the predicted scores only
        res.predict = predict.mixo_spls(object, newdata = newdata,
                                        dist = "max.dist", multilevel = multilevel)$predict
        
        for (i in seq_len(object$ncomp))
        {
            data$data=res.predict[,,i]
            temp = statauc(data, plot = ifelse(i%in%roc.comp,plot,FALSE), title=title,...)
            statauc.res[[paste0("Comp", i, sep = "")]] = temp[[1]]
            graph[[paste0("Comp", i, sep = "")]] = temp$graph
        }
        if (isTRUE(print))
            print(statauc.res)
        
        return(invisible(c(statauc.res,graph=graph)))
    }

#' @rdname auroc
#' @export
auroc.mixo_splsda <- auroc.mixo_plsda

# MINT object
# ----------------------
#' @rdname auroc
#' @method auroc mint.plsda
#' @export
auroc.mint.plsda <- 
    function(
        object,
        newdata = object$X,
        outcome.test = as.factor(object$Y),
        study.test = object$study,
        multilevel = NULL,
        plot = TRUE,
        roc.comp = 1,
        roc.study = "global",
        title=NULL,
        print=TRUE,
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
        statauc.res = graph = list()
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
            statauc.res[[paste0("Comp", i, sep = "")]] = temp[[1]]
            graph[[paste0("Comp", i, sep = "")]] = temp$graph
        }
        if (isTRUE(print))
            print(statauc.res)
        
        return(invisible(c(statauc.res,graph=graph)))
        
    }

#' @rdname auroc
#' @method auroc mint.splsda
#' @export
auroc.mint.splsda <- auroc.mint.plsda

# block.splsda object
# ----------------------
#' @rdname auroc
#' @importFrom methods is
#' @method auroc sgccda
#' @export
auroc.sgccda <- function(
    object,
    newdata = object$X,
    outcome.test = as.factor(object$Y),
    multilevel = NULL,
    plot = TRUE,
    roc.block = 1L,
    roc.comp = 1L,
    title=NULL,
    print=TRUE,
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
    if (isTRUE(print))
        print(auc.mean)
    
    return(invisible(out))
}

# mint.block.splsda object
# ----------------------
#' @rdname auroc
#' @method auroc mint.block.plsda
#' @export
auroc.mint.block.plsda <- function(
    object,
    newdata = object$X,
    study.test = object$study,
    outcome.test = as.factor(object$Y),
    multilevel = NULL,
    plot = TRUE,
    roc.block = 1,
    roc.comp = 1,
    title=NULL,
    print=TRUE,
    ...)
{
    
    data=list()
    auc.mean = graph=list()
    data$outcome=factor(outcome.test)
    study.test=factor(study.test)
    
    # note here: the dist does not matter as we used the predicted scores only
    res.predict  =  predict.mixo_spls(object, newdata = newdata,
                                      study.test=study.test,dist = "max.dist", multilevel = multilevel)$predict
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
    if (isTRUE(print))
        print(auc.mean)
    
    return(invisible(out))
}

#' @rdname auroc
#' @method auroc mint.block.splsda
#' @export
auroc.mint.block.splsda <- auroc.mint.block.plsda
