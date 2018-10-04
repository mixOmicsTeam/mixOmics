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

# X: a list of data sets (called 'blocks') matching on the same samples. Data in the list should be arranged in samples x variables, with samples order matching in all data sets. \code{NA}s are not allowed.
# Y: a factor or a class vector for the discrete outcome.
# indY: to supply if Y is missing, indicate the position of the outcome in the list X.
# ncomp: numeric vector of length the number of blocks in \code{X}. The number of components to include in the model for each block (does not necessarily need to take the same value for each block). By default set to 2 per block.
# test.keepX: list of length, length(X). each test.keepX[[i]] is a grid of keepX among which to chose the optimal one
# already.tested.X: list of length, length(X). Each already.tested.X[[i]] is a vector giving the keepX on the components that were already tuned
# validation: Mfold or loo cross validation
# folds: if validation=Mfold, how many folds?
# dist: distance to classify samples. see predict
# measure: one of c("overall","BER"). Accuracy measure used in the cross validation processs
# weighted: optimise the weighted or not-weighted prediction
# progressBar: show progress,
# tol: Convergence stopping value.
# max.iter: integer, the maximum number of iterations.
# near.zero.var: boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations
# nrepeat: number of replication of the Mfold process
# design: the input design.
# scheme: the input scheme, one of "horst", "factorial" or ""centroid". Default to "centroid"
# scale: boleean. If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
# init: intialisation of the algorithm, one of "svd" or "svd.single". Default to "svd"
# light.output: if FALSE, output the classification of each sample during each folds, on each comp, for each repeat
# cpus: number of cpus to use. default to no parallel
# name.save: if saving a file after each component


tune.block.splsda = function (
X,
Y,
indY,
ncomp = 2,
test.keepX,
already.tested.X,
validation = "Mfold",
folds = 10,
dist = "max.dist",
measure = "BER", # one of c("overall","BER")
weighted = TRUE, # optimise the weighted or not-weighted prediction
progressBar = TRUE,
tol = 1e-06,
max.iter = 100,
near.zero.var = FALSE,
nrepeat = 1,
design,
scheme= "horst",
scale = TRUE,
init = "svd",
light.output = TRUE, # if FALSE, output the prediction and classification of each sample during each folds, on each comp, for each repeat
cpus,
name.save = NULL)
{
    #-- checking general input parameters --------------------------------------#
    #---------------------------------------------------------------------------#
    
    #------------------#
    #-- check entries --#
    
    
    # check inpuy 'Y' and transformation in a dummy matrix
    if(!missing(Y))
    {
        if (is.null(dim(Y)))
        {
            Y = factor(Y)
        } else {
            stop("'Y' should be a factor or a class vector.")
        }
        
        if (nlevels(Y) == 1)
        stop("'Y' should be a factor with more than one level")
        
    } else if(!missing(indY)) {
        Y = X[[indY]]
        if (is.null(dim(Y)))
        {
            Y = factor(Y)
        } else {
            stop("'Y' should be a factor or a class vector.")
        }
        
        if (nlevels(temp) == 1)
        stop("'X[[indY]]' should be a factor with more than one level")
        
        X = X[-indY] #remove Y from X to pass the arguments simpler to block.splsda
        
    } else if(missing(indY)) {
        stop("Either 'Y' or 'indY' is needed")
        
    }
    
    # there's a X and a Y, we force the data to be matrices
    cl = lapply(X,class)
    ind.no.matrix = which(sapply(cl, function(x) !any(x == "matrix")))
    if(length(ind.no.matrix)>0){
        X[ind.no.matrix] = lapply(X[ind.no.matrix], function(x) as.matrix(x))
    }

    
    #-- dist
    dist = match.arg(dist, choices = c("max.dist", "centroids.dist", "mahalanobis.dist"), several.ok = FALSE)

    #-- progressBar
    if (!is.logical(progressBar))
    stop("'progressBar' must be a logical constant (TRUE or FALSE).", call. = FALSE)
    
    #-- ncomp
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
    
    
    #-- measure
    measure.input = measure
    if(! measure %in% c("overall", "BER"))
    stop("'measure' must be 'overall' or 'BER'")
    
    
    #-- already.tested.X

    if (missing(already.tested.X))
    {
        already.tested.X = NULL
    } else {

        if(is.null(already.tested.X))
        stop("'already.tested.X' must be a vector of keepX values ")
        
        # we require the same number of already tuned components on each block
        if(length(unique(sapply(already.tested.X, length))) > 1)
        stop("The same number of components must be already tuned for each block, in 'already.tested.X'")
        
        if(any(sapply(already.tested.X, function(x) is.list(x))) == TRUE)
        stop(" Each entry of 'already.tested.X' must be a vector of keepX values")

        if(length(already.tested.X[[1]]) >= ncomp)
        stop("'ncomp' needs to be higher than the number of components already tuned, which is length(already.tested.X)=",length(already.tested.X) , call. = FALSE)
    }
    
    if (any(is.na(validation)) || length(validation) > 1)
    stop("'validation' should be one of 'Mfold' or 'loo'.", call. = FALSE)
    
    #-- test.keepX
    if(missing(test.keepX))
    {
        test.keepX = lapply(1:length(X),function(x){c(5,10,15)[which(c(5,10,15)<ncol(X[[x]]))]})
        names(test.keepX) = names(X)
        
    } else {
        if(length(test.keepX) != length(X))
        stop(paste("test.keepX should be a list of length ", length(X),", corresponding to the blocks: ", paste(names(X),collapse=", "), sep=""))
        
        aa = sapply(test.keepX, length)
        if (any(is.null(aa) | aa == 1 | !is.numeric(aa)))
        stop("Each entry of 'test.keepX' must be a numeric vector with more than two values", call. = FALSE)
        
    }
    
    l = sapply(test.keepX,length)
    n = names(test.keepX)
    temp = data.frame(l, n)
    
    
    message(paste("You have provided a sequence of keepX of length: ", paste(apply(temp, 1, function(x) paste(x,collapse=" for block ")), collapse= " and "), ".\nThis results in ",prod(sapply(test.keepX,length)), " models being fitted for each component and each nrepeat, this may take some time to run, be patient!",sep=""))
    
    if(missing(cpus))
    {
        parallel = FALSE
        message(paste("You can look into the 'cpus' argument to speed up computation time.",sep=""))

    } else {
        parallel = TRUE
        if(progressBar == TRUE)
        message(paste("As code is running in parallel, the progressBar will only show 100% upon completion of each nrepeat/ component.",sep=""))

    }
    
    #-- end checking --#
    #------------------#
    
    
    #---------------------------------------------------------------------------#
    #-- NA calculation      ----------------------------------------------------#
    
    misdata = c(sapply(X,anyNA), Y=FALSE) # Detection of missing data. we assume no missing values in the factor Y
    
    is.na.A = vector("list", length = length(X))
    for(q in 1:length(X))
    {
        if(misdata[q])
        {
            is.na.A[[q]] = is.na(X[[q]])
            #ind.NA[[q]] = which(apply(is.na.A[[q]], 1, sum) > 0) # calculated only once
            #ind.NA.col[[q]] = which(apply(is.na.A[[q]], 2, sum) >0) # indice of the col that have missing values. used in the deflation
       }
    }

    #-- NA calculation      ----------------------------------------------------#
    #---------------------------------------------------------------------------#
    
    
    # if some components have already been tuned (eg comp1 and comp2), we're only tuning the following ones (comp3 comp4 .. ncomp)
    if ((!is.null(already.tested.X)) & length(already.tested.X) > 0)
    {
        comp.real = (length(already.tested.X[[1]]) + 1):ncomp
        #check and match already.tested.X to X
        if(length(already.tested.X[[1]]) >0)
        {
            if(length(unique(names(already.tested.X)))!=length(already.tested.X) | sum(is.na(match(names(already.tested.X),names(X)))) > 0)
            stop("Each entry of 'already.tested.X' must have a unique name corresponding to a block of 'X'")
            
        }
        
    } else {
        comp.real = 1:ncomp
    }
    
    # near zero var on the whole data sets. It will be performed inside each fold as well
    if(near.zero.var == TRUE)
    {
        nzv.A = lapply(X, nearZeroVar)
        for(q in 1:length(X))
        {
            if (length(nzv.A[[q]]$Position) > 0)
            {
                names.remove.X = colnames(X[[q]])[nzv.A[[q]]$Position]
                X[[q]] = X[[q]][, -nzv.A[[q]]$Position, drop=FALSE]
                warning("Zero- or near-zero variance predictors.\n Reset predictors matrix to not near-zero variance predictors.\n See $nzv for problematic predictors.")
                if (ncol(X[[q]]) == 0)
                stop(paste0("No more variables in",X[[q]]))
                
                #need to check that the keepA[[q]] is now not higher than ncol(A[[q]])
                if (any(test.keepX[[q]] > ncol(X[[q]])))
                test.keepX[[q]][which(test.keepX[[q]] > ncol(X[[q]]))] = ncol(X[[q]])
            }
            
        }
    }
    
    
    if (parallel == TRUE)
    {
        cl <- makeCluster(cpus, type = "SOCK")
        clusterEvalQ(cl, library(mixOmics))
    } else{cl=NULL}
    
    
    N.test.keepX = nrow(expand.grid(test.keepX))
    
    mat.error.rate = list()
    error.per.class = list()
    
    mat.sd.error = matrix(0,nrow = N.test.keepX, ncol = ncomp-length(already.tested.X[[1]]))#,
    #    dimnames = list(c(test.keepX), c(paste('comp', comp.real, sep=''))))
    mat.mean.error = matrix(nrow = N.test.keepX, ncol = ncomp-length(already.tested.X[[1]]))#,
    #dimnames = list(c(test.keepX), c(paste('comp', comp.real, sep=''))))
    
    
    mat.error.rate = list()
    error.per.class.keepX.opt=list()
    error.per.class.keepX.opt.mean = matrix(0, nrow = nlevels(Y), ncol = length(comp.real),
    dimnames = list(c(levels(Y)), c(paste('comp', comp.real, sep=''))))

    error.opt.per.comp = matrix(nrow = nrepeat, ncol = length(comp.real), dimnames=list(paste("nrep",1:nrepeat,sep="."), paste("comp", comp.real, sep='')))
    
    if(light.output == FALSE)
    prediction.all = class.all = list()

    #save(list=ls(),file="temp3.Rdata")
    # successively tune the components until ncomp: comp1, then comp2, ...
    for(comp in 1:length(comp.real))
    {
        if (progressBar == TRUE)
        cat("\ncomp",comp.real[comp], "\n")
        
        result = MCVfold.block.splsda (X, Y, validation = validation, folds = folds, nrepeat = nrepeat, ncomp = 1 + length(already.tested.X[[1]]),
        choice.keepX = already.tested.X, scheme = scheme, design=design, init=init, tol=tol,
        test.keepX = test.keepX, measure = measure, dist = dist, scale=scale, weighted=weighted,
        near.zero.var = near.zero.var, progressBar = progressBar, max.iter = max.iter, cl = cl,
        misdata = misdata, is.na.A = is.na.A, parallel = parallel)
        
        #returns error.rate for all test.keepX
    
        
        # in the following, there is [[1]] because 'tune' is working with only 1 distance and 'MCVfold.block.splsda' can work with multiple distances
        mat.error.rate[[comp]] = result[[measure]]$mat.error.rate[[1]]
        mat.mean.error[, comp]=result[[measure]]$error.rate.mean[[1]]
        if (!is.null(result[[measure]]$error.rate.sd[[1]]))
        mat.sd.error[, comp]=result[[measure]]$error.rate.sd[[1]]
        
        # confusion matrix for keepX.opt
        error.per.class.keepX.opt[[comp]]=result[[measure]]$confusion[[1]]
        error.per.class.keepX.opt.mean[, comp]=apply(result[[measure]]$confusion[[1]], 1, mean)
        
        # error rate for best keepX
        error.opt.per.comp[,comp] = mat.error.rate[[comp]][result[[measure]]$ind.keepX.opt[[1]],]
        
        # best keepX
        already.tested.X = result[[measure]]$choice.keepX
        
        if(light.output == FALSE)
        {
            #prediction of each samples for each fold and each repeat, on each comp
            class.all[[comp]] = result$class.comp[[1]]
            #prediction.all[[comp]] = result$prediction.comp
        }
        
        # prepping the results and save a file, if necessary
        if(!is.null(name.save))
        {
            rownames(mat.mean.error) = rownames(result[[measure]]$mat.error.rate[[1]])
            colnames(mat.mean.error) = paste("comp", comp.real, sep='')
            names(mat.error.rate) = c(paste("comp", comp.real[1:comp], sep=''))
            names(error.per.class.keepX.opt) = c(paste("comp", comp.real[1:comp], sep=''))
            if(nrepeat > 1)
            {
                rownames(mat.sd.error) = rownames(result[[measure]]$mat.error.rate[[1]])
                colnames(mat.sd.error) = paste("comp", comp.real, sep='')
            }
            
            
            result = list(
            error.rate = mat.mean.error,
            error.rate.sd = mat.sd.error,
            error.rate.all = mat.error.rate,
            choice.keepX = already.tested.X,
            error.rate.class = error.per.class.keepX.opt)
            
            result$measure = measure.input
            result$call = match.call()
            
            class(result) = "tune.block.splsda"
            
            save(result, file = paste0(name.save,".comp",comp.real[1],"to",comp.real[comp],".Rdata"))
        }

    }
    rownames(mat.mean.error) = rownames(result[[measure]]$mat.error.rate[[1]])
    colnames(mat.mean.error) = paste("comp", comp.real, sep='')
    names(mat.error.rate) = c(paste("comp", comp.real, sep=''))
    names(error.per.class.keepX.opt) = c(paste("comp", comp.real, sep=''))
    if(nrepeat > 1)
    {
        rownames(mat.sd.error) = rownames(result[[measure]]$mat.error.rate[[1]])
        colnames(mat.sd.error) = paste("comp", comp.real, sep='')
    }
    
    #close the cluster after ncomp
    if (parallel == TRUE)
    stopCluster(cl)
    
    cat("\n")
    
    
    # calculating the number of optimal component based on t.tests and the error.rate.all, if more than 3 error.rates(repeat>3)
    if(nrepeat > 2 & length(comp.real) >1)
    {
       error.keepX = error.opt.per.comp
       opt = t.test.process(error.opt.per.comp)
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
    error.rate.class = error.per.class.keepX.opt)
    
    result$measure = measure.input
    result$call = match.call()
    
    class(result) = "tune.block.splsda"
    
    return(result)
    
}

