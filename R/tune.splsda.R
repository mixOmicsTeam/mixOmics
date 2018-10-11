#############################################################################################################
# Authors:
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Benoit Gautier, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Francois Bartolo, Institut National des Sciences Appliquees et Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
#
# created: 2013
# last modified: 05-10-2017
#
# Copyright (C) 2013
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


# ========================================================================================================
# tune.splsda: chose the optimal number of parameters per component on a splsda method
# ========================================================================================================

# X: numeric matrix of predictors
# Y: a factor or a class vector for the discrete outcome
# ncomp: the number of components to include in the model. Default to 1.
# test.keepX: grid of keepX among which to chose the optimal one
# already.tested.X: a vector giving keepX on the components that were already tuned
# validation: Mfold or loo cross validation
# folds: if validation=Mfold, how many folds?
# dist: distance to classify samples. see predict
# measure: one of c("overall","BER"). Accuracy measure used in the cross validation processs
# scale: boleean. If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
# auc: calculate AUC
# progressBar: show progress,
# tol: Convergence stopping value.
# max.iter: integer, the maximum number of iterations.
# near.zero.var: boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations
# nrepeat: number of replication of the Mfold process
# logratio: one of "none", "CLR"
# multilevel: repeated measurement. `multilevel' is passed to multilevel(design = ) in withinVariation. Y is ommited and shouldbe included in `multilevel'
# light.output: if FALSE, output the prediction and classification of each sample during each folds, on each comp, for each repeat
# cpus: number of cpus to use. default to no parallel

tune.splsda = function (X, Y,
ncomp = 1,
test.keepX = c(5, 10, 15),
already.tested.X,
validation = "Mfold",
folds = 10,
dist = "max.dist",
measure = "BER", # one of c("overall","BER")
scale = TRUE,
auc = FALSE,
progressBar = TRUE,
tol = 1e-06,
max.iter = 100,
near.zero.var = FALSE,
nrepeat = 1,
logratio = c('none','CLR'),
multilevel = NULL,
light.output = TRUE,
cpus
)
{    #-- checking general input parameters --------------------------------------#
    #---------------------------------------------------------------------------#
    
    #------------------#
    #-- check entries --#
    if(!any(class(X) == "matrix"))
    X = as.matrix(X)
    
    if (length(dim(X)) != 2 || !is.numeric(X))
    stop("'X' must be a numeric matrix.")
    
    
    # Testing the input Y
    if (is.null(multilevel))
    {
        if (is.null(Y))
        stop("'Y' has to be something else than NULL.")
        
        if (is.null(dim(Y)))
        {
            Y = factor(Y)
        }  else {
            stop("'Y' should be a factor or a class vector.")
        }
        
        if (nlevels(Y) == 1)
        stop("'Y' should be a factor with more than one level")
        
    } else {
        # we expect a vector or a 2-columns matrix in 'Y' and the repeated measurements in 'multilevel'
        multilevel = data.frame(multilevel)
        
        if ((nrow(X) != nrow(multilevel)))
        stop("unequal number of rows in 'X' and 'multilevel'.")
        
        if (ncol(multilevel) != 1)
        stop("'multilevel' should have a single column for the repeated measurements, other factors should be included in 'Y'.")
        
        if (!is.null(ncol(Y)) && !ncol(Y) %in% c(0,1,2))# multilevel 1 or 2 factors
        stop("'Y' should either be a factor, a single column data.frame containing a factor, or a 2-columns data.frame containing 2 factors.")
        
        multilevel = data.frame(multilevel, Y)
        multilevel[, 1] = as.numeric(factor(multilevel[, 1])) # we want numbers for the repeated measurements
        
    }
    
    
    #-- progressBar
    if (!is.logical(progressBar))
    stop("'progressBar' must be a logical constant (TRUE or FALSE).", call. = FALSE)
    
    
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0)
    stop("invalid number of variates, 'ncomp'.")
    
    
    #-- validation
    choices = c("Mfold", "loo")
    validation = choices[pmatch(validation, choices)]
    if (is.na(validation))
    stop("'validation' must be either 'Mfold' or 'loo'")
    
    if (validation == 'loo')
    {
        if (nrepeat != 1)
        warning("Leave-One-Out validation does not need to be repeated: 'nrepeat' is set to '1'.")
        nrepeat = 1
    }
    
    #-- logratio
    if (length(logratio) > 1)
    logratio = logratio[1]
    
    if (!logratio %in% c("none","CLR"))
    stop("'logratio must be one of 'none' or 'CLR'")
    
    #if ((!is.null(already.tested.X)) && (length(already.tested.X) != (ncomp - 1)) )
    #stop("The number of already tested parameters should be NULL or ", ncomp - 1, " since you set ncomp = ", ncomp)
    
    if (missing(already.tested.X))
    {
        already.tested.X = NULL
    } else {
        if(is.null(already.tested.X) | length(already.tested.X)==0)
        stop("''already.tested.X' must be a vector of keepX values")

        if(is.list(already.tested.X))
        stop("''already.tested.X' must be a vector of keepX values")

        message(paste("Number of variables selected on the first", length(already.tested.X), "component(s):", paste(already.tested.X,collapse = " ")))
    }
    
    if(length(already.tested.X) >= ncomp)
    stop("'ncomp' needs to be higher than the number of components already tuned, which is length(already.tested.X)=",length(already.tested.X) , call. = FALSE)
    
    #-- measure
    choices = c("BER", "overall","AUC")
    measure = choices[pmatch(measure, choices)]
    if (is.na(measure))
    stop("'measure' must be either 'BER', 'overall' or 'AUC' ")

    if (any(is.na(validation)) || length(validation) > 1)
    stop("'validation' should be one of 'Mfold' or 'loo'.", call. = FALSE)
    
    #-- test.keepX
    if (is.null(test.keepX) | length(test.keepX) == 1 | !is.numeric(test.keepX))
    stop("'test.keepX' must be a numeric vector with more than two entries", call. = FALSE)
    
    # remove some test.keepX if needed
    if (any(test.keepX > ncol(X))){
        test.keepX = test.keepX[-which(test.keepX>ncol(X))]
        if (length(test.keepX) < 2)
        stop("Some entries of 'test.keepX' were higher than the number of
        variables in 'X' and were removed, the resulting 'test.keepX' has now
        too few entries (<2)", call. = FALSE)
    }
    
    
    if(!missing(cpus))
    {
        if(!is.numeric(cpus) | length(cpus)!=1)
        stop("'cpus' must be a numerical value")
        
        parallel = TRUE
        cl = makeCluster(cpus, type = "SOCK")
        #clusterExport(cl, c("splsda","selectVar"))
        clusterEvalQ(cl, library(mixOmics))

        if(progressBar == TRUE)
        message(paste("As code is running in parallel, the progressBar will only show 100% upon completion of each nrepeat/ component.",sep=""))

    } else {
        parallel = FALSE
        cl = NULL
    }
    
    # add colnames and rownames if missing
    X.names = dimnames(X)[[2]]
    if (is.null(X.names))
    {
        X.names = paste("X", 1:ncol(X), sep = "")
        dimnames(X)[[2]] = X.names
    }
    
    ind.names = dimnames(X)[[1]]
    if (is.null(ind.names))
    {
        ind.names = 1:nrow(X)
        rownames(X)  = ind.names
    }
    
    if (length(unique(rownames(X))) != nrow(X))
    stop("samples should have a unique identifier/rowname")
    if (length(unique(X.names)) != ncol(X))
    stop("Unique indentifier is needed for the columns of X")
    
    rm(X.names);rm(ind.names)

    #-- end checking --#
    #------------------#
    
   
    #---------------------------------------------------------------------------#
    #-- logration + multilevel approach ----------------------------------------#
    # we can do logratio and multilevel on the whole data as these transformation are done per sample
    X = logratio.transfo(X = X, logratio = logratio)
    
    if (!is.null(multilevel)) # logratio is applied per sample, multilevel as well, so both can be done on the whole data
    {

        Xw = withinVariation(X, design = multilevel)
        X = Xw

        #-- Need to set Y variable for 1 or 2 factors
        Y = multilevel[, -1, drop=FALSE]
        if (ncol(Y) >= 1)
        Y = apply(Y, 1, paste, collapse = ".")  #  paste is to combine in the case we have 2 levels
        
        Y = as.factor(Y)
    }
    #-- logration + multilevel approach ----------------------------------------#
    #---------------------------------------------------------------------------#
    
    
    #---------------------------------------------------------------------------#
    #-- NA calculation      ----------------------------------------------------#
    
    misdata = c(X=anyNA(X), Y=FALSE) # Detection of missing data. we assume no missing values in the factor Y
    
    if (any(misdata))
    {
        is.na.A = is.na(X)
        
        #ind.NA = which(apply(is.na.A, 1, sum) > 0) # calculated only once
        #ind.NA.col = which(apply(is.na.A, 2, sum) >0) # indice of the col that have missing values. used in the deflation
    } else {
        is.na.A = NULL
        #ind.NA = ind.NA.col = NULL
    }
    #-- NA calculation      ----------------------------------------------------#
    #---------------------------------------------------------------------------#



    test.keepX = sort(unique(test.keepX)) #sort test.keepX so as to be sure to chose the smallest in case of several minimum
    names(test.keepX) = test.keepX
    # if some components have already been tuned (eg comp1 and comp2), we're only tuning the following ones (comp3 comp4 .. ncomp)
    if ((!is.null(already.tested.X)))
    {
        comp.real = (length(already.tested.X) + 1):ncomp
    } else {
        comp.real = 1:ncomp
    }
    
    
    choices = c("all", "max.dist", "centroids.dist", "mahalanobis.dist")
    dist = match.arg(dist, choices, several.ok = TRUE)
    
    
    mat.error.rate = list()
    error.per.class = list()
    AUC = list()
    
    mat.sd.error = matrix(0,nrow = length(test.keepX), ncol = ncomp-length(already.tested.X),
    dimnames = list(c(test.keepX), c(paste('comp', comp.real, sep=''))))
    mat.mean.error = matrix(nrow = length(test.keepX), ncol = ncomp-length(already.tested.X),
    dimnames = list(c(test.keepX), c(paste('comp', comp.real, sep=''))))        
   
    # first: near zero var on the whole data set
    if(near.zero.var == TRUE)
    {
        nzv = nearZeroVar(X)
        if (length(nzv$Position > 0))
        {
            warning("Zero- or near-zero variance predictors.\nReset predictors matrix to not near-zero variance predictors.\nSee $nzv for problematic predictors.")
            X = X[, -nzv$Position, drop=TRUE]
            
            if(ncol(X)==0)
            stop("No more predictors after Near Zero Var has been applied!")
            
        }
    }

    if (is.null(multilevel) | (!is.null(multilevel) && ncol(multilevel) == 2))
    {
        if(light.output == FALSE)
        prediction.all = class.all = list()
        
        if(auc)
        {
            auc.mean.sd=list()
            if(light.output == FALSE)
            auc.all=list()
        }
        
        class.object=c("mixo_splsda","DA")
        if(!missing(cpus))
            clusterExport(cl, c("X","Y","is.na.A","misdata","scale","near.zero.var","class.object","test.keepX"),envir=environment())

        error.per.class.keepX.opt = list()
        error.per.class.keepX.opt.mean = matrix(0, nrow = nlevels(Y), ncol = length(comp.real),
        dimnames = list(c(levels(Y)), c(paste('comp', comp.real, sep=''))))
        # successively tune the components until ncomp: comp1, then comp2, ...
        for(comp in 1:length(comp.real))
        {

            if (progressBar == TRUE)
            cat("\ncomp",comp.real[comp], "\n")
            
            result = MCVfold.spls (X, Y, multilevel = multilevel, validation = validation, folds = folds, nrepeat = nrepeat, ncomp = 1 + length(already.tested.X),
            choice.keepX = already.tested.X,
            test.keepX = test.keepX, test.keepY = nlevels(Y), measure = measure, dist = dist, scale=scale,
            near.zero.var = near.zero.var, progressBar = progressBar, tol = tol, max.iter = max.iter, auc = auc,
            cl = cl, parallel = parallel,
            misdata = misdata, is.na.A = is.na.A, class.object=class.object)
            
            # in the following, there is [[1]] because 'tune' is working with only 1 distance and 'MCVfold.splsda' can work with multiple distances
            mat.error.rate[[comp]] = result[[measure]]$mat.error.rate[[1]]
            mat.mean.error[, comp]=result[[measure]]$error.rate.mean[[1]]
            if (!is.null(result[[measure]]$error.rate.sd[[1]]))
            mat.sd.error[, comp]=result[[measure]]$error.rate.sd[[1]]
            
            # confusion matrix for keepX.opt
            if (all(measure!="AUC")){
                error.per.class.keepX.opt[[comp]]=result[[measure]]$confusion[[1]]
                error.per.class.keepX.opt.mean[, comp]=apply(result[[measure]]$confusion[[1]], 1, mean)
            }

            # best keepX
            already.tested.X = c(already.tested.X, result[[measure]]$keepX.opt[[1]])
             
            if(light.output == FALSE)
            {
                #prediction of each samples for each fold and each repeat, on each comp
                class.all[[comp]] = result$class.comp[[1]]
                prediction.all[[comp]] = result$prediction.comp
            }
            
            if(auc)
            {
                auc.mean.sd[[comp]] = result$auc
                if(light.output == FALSE)
                auc.all[[comp]] = result$auc.all
            }
            
        } # end comp
        if (parallel == TRUE)
        stopCluster(cl)
        names(mat.error.rate) = c(paste('comp', comp.real, sep=''))
        if (all(measure!="AUC"))
            names(error.per.class.keepX.opt) = c(paste('comp', comp.real, sep=''))
        
        names(already.tested.X) = c(paste('comp', 1:ncomp, sep=''))
        
        if (progressBar == TRUE)
        cat('\n')

        # calculating the number of optimal component based on t.tests and the error.rate.all, if more than 3 error.rates(repeat>3)
        if(nrepeat > 2 & length(comp.real) >1 & all(measure!="AUC"))
        {
            keepX = already.tested.X
            error.keepX = NULL
            for(comp in 1:length(comp.real))
            {
                ind.row = match(keepX[[comp.real[comp]]],test.keepX)
                error.keepX = cbind(error.keepX, mat.error.rate[[comp]][ind.row,])
            }
            colnames(error.keepX) = c(paste('comp', comp.real, sep=''))
            
            opt = t.test.process(error.keepX)
            
            ncomp_opt = comp.real[opt]
        } else {
            ncomp_opt = error.keepX = NULL
        }
        
        result = list(
        error.rate = mat.mean.error,
        error.rate.sd = mat.sd.error,
        error.rate.all = mat.error.rate,
        choice.keepX = already.tested.X,
        choice.ncomp = list(ncomp = ncomp_opt, values = error.keepX),
        error.rate.class = error.per.class.keepX.opt.mean,
        error.rate.class.all = error.per.class.keepX.opt)

        if(light.output == FALSE)
        {
            names(class.all) = names(prediction.all) = c(paste('comp', comp.real, sep=''))
            result$predict = prediction.all
            result$class = class.all
        }
        if(auc)
        {
            
            #we add the AUC outputs only if it was not the measure
            #otherwise there is twice the same output
            if(all(measure != "AUC")){
                names(auc.mean.sd) = c(paste('comp', comp.real, sep=''))
                result$auc = auc.mean.sd
            }
            if(light.output == FALSE)
            {
                names(auc.all) = c(paste('comp', comp.real, sep=''))
                result$auc.all =auc.all
            }
        }
        result$measure = measure
        result$call = match.call()

        class(result) = "tune.splsda"

        return(result)
    } else {
        # if multilevel with 2 factors, we can not do as before because withinvariation depends on the factors, we maximase a correlation
        message("For a two-factor analysis, the tuning criterion is based on the maximisation of the correlation between the components on the whole data set")

        cor.value = vector(length = length(test.keepX))
        names(cor.value) = test.keepX

        for (i in 1:length(test.keepX))
        {
            spls.train = mixOmics::splsda(X, Y, ncomp = ncomp, keepX = c(already.tested.X, test.keepX[i]), logratio = logratio, near.zero.var = FALSE, mode = "regression")
            
            # Note: this is performed on the full data set
            # (could be done with resampling (bootstrap) (option 1) and/or prediction (option 2))
            cor.value[i] = cor(spls.train$variates$X[, ncomp], spls.train$variates$Y[, ncomp])
            #
        }
        return(list(cor.value = cor.value))

    }
}
