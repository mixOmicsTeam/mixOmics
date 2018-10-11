#############################################################################################################
# Authors:
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Francois Bartolo, Institut National des Sciences Appliquees et Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 24-08-2016
# last modified: 05-10-2017
#
# Copyright (C) 2015
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
# multilevel: repeated measurement only
# validation: Mfold or loo cross validation
# folds: if validation = Mfold, how many folds?
# nrepeat: number of replication of the Mfold process
# ncomp: the number of components to include in the model. Default to 1.
# choice.keepX: a vector giving keepX on the components that were already tuned
# test.keepX: grid of keepX among which to chose the optimal one
# measure: one of c("overall","BER"). Accuracy measure used in the cross validation processs
# dist: distance to classify samples. see predict
# auc:
# tol: Convergence stopping value.
# max.iter: integer, the maximum number of iterations.
# near.zero.var: boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations
# progressBar: show progress,
# class.object
# cl: if parallel, the clusters
# scale: boleean. If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
# misdata: optional. any missing values in the data? list, misdata[[q]] for each data set
# is.na.A: optional. where are the missing values? list, is.na.A[[q]] for each data set (if misdata[[q]] == TRUE)
# ind.NA: optional. which rows have missing values? list, ind.NA[[q]] for each data set.
# ind.NA.col: optional. which col have missing values? list, ind.NA.col[[q]] for each data set.
# parallel: logical.

stratified.subsampling = function(Y, folds = 10)
{
    stop = 0
    for(i in 1:nlevels(Y))
    {
        ai=sample(which(Y==levels(Y)[i]),replace=FALSE) # random sampling of the samples from level i
        aai=suppressWarnings(split(ai,factor(1:min(folds,length(ai)))))                       # split of the samples in k-folds
        if(length(ai)<folds)                                                # if one level doesn't have at least k samples, the list is completed with "integer(0)"
        {
            for(j in (length(ai)+1):folds)
            aai[[j]]=integer(0)
            stop = stop +1
        }
        assign(paste("aa",i,sep="_"),sample(aai,replace=FALSE))         # the `sample(aai)' is to avoid the first group to have a lot more data than the rest
    }
    
    # combination of the different split aa_i into SAMPLE
    SAMPLE=list()
    for(j in 1:folds)
    {
        SAMPLE[[j]]=integer(0)
        for(i in 1:nlevels(Y))
        {
            SAMPLE[[j]]=c(SAMPLE[[j]],get(paste("aa",i,sep="_"))[[j]])
        }
    }# SAMPLE is a list of k splits
    
    ind0 = sapply(SAMPLE, length)
    if(any(ind0 == 0))
    {
        SAMPLE = SAMPLE [-which(ind0 == 0)]
        message("Because of a too high number of 'folds' required, ",length(which(ind0 == 0))," folds were randomly assigned no data: the number of 'folds' is reduced to ", length(SAMPLE))
    }
    
    return(list(SAMPLE = SAMPLE, stop = stop))
}


MCVfold.spls = function(
X,
Y,
multilevel = NULL, # repeated measurement only
validation,
folds,
nrepeat = 1,
ncomp,
choice.keepX = NULL, # keepX chosen on the first components
test.keepX, # a vector of value(keepX) to test on the last component. There needs to be names(test.keepX)
test.keepY, # a vector of value(keepX) to test on the last component. There needs to be names(test.keepX)
measure = c("overall"), # one of c("overall","BER")
dist = "max.dist",
auc = FALSE,
tol = 1e-06,
max.iter = 100,
near.zero.var = FALSE,
progressBar = TRUE,
class.object,
cl,
scale,
misdata,
is.na.A,
#ind.NA,
#ind.NA.col,
parallel
)
{    #-- checking general input parameters --------------------------------------#
    #---------------------------------------------------------------------------#
    #-- set up a progress bar --#
    if (progressBar ==  TRUE)
    {
        pb = txtProgressBar(style = 3)
        nBar = 1
    } else {
        pb = FALSE
    }
    if(any(measure == "AUC")) auc=TRUE

    design = matrix(c(0,1,1,0), ncol = 2, nrow = 2, byrow = TRUE)
    
    if(ncomp>1 &  any(class.object == "DA")){keepY = rep(nlevels(Y), ncomp-1)} else {keepY=rep(ncol(Y),ncomp-1)}
    
    rownames.X = rownames(X)
    M = length(folds)
    features = features.j = NULL
    auc.all = prediction.comp = class.comp = list()
    

    if(any(class.object == "DA")){
        for(ijk in dist)
        class.comp[[ijk]] = array(0, c(nrow(X), nrepeat, length(test.keepX)))# prediction of all samples for each test.keepX and  nrep at comp fixed
    } else {
        prediction.keepX = vector("list", length=length(test.keepX))
        for(i in 1:length(test.keepX))
        prediction.keepX[[i]] = array(0,c(nrow(X), ncol(Y), nrepeat), dimnames = list(rownames(X), colnames(Y), paste0("nrep.", 1:nrepeat)))
    }
    
    folds.input = folds # save fold number to be used for each nrepeat
    for(nrep in 1:nrepeat)
    {
        n = nrow(X)
        if(any(class.object == "DA")){
            prediction.comp[[nrep]] = array(0, c(n, nlevels(Y), length(test.keepX)), dimnames = list(rownames.X, levels(Y), names(test.keepX)))
            colnames(prediction.comp[[nrep]]) = levels(Y)
            
            if(nlevels(Y)>2)
            {
                auc.all[[nrep]] = array(0, c(nlevels(Y),2, length(test.keepX)), dimnames = list(paste(levels(Y), "vs Other(s)"), c("AUC","p-value"), names(test.keepX)))
            }else{
                auc.all[[nrep]] = array(0, c(1,2, length(test.keepX)), dimnames = list(paste(levels(Y)[1], levels(Y)[2], sep = " vs "), c("AUC","p-value"), names(test.keepX)))
            }
        } else {
            prediction.comp[[nrep]] = array(0, c(nrow(X), ncol(Y), length(test.keepX)), dimnames = list(rownames(X), colnames(Y), test.keepX))
            colnames(prediction.comp[[nrep]]) = colnames(Y)
        }
        rownames(prediction.comp[[nrep]]) = rownames.X
        
        
        repeated.measure = 1:n
        if (!is.null(multilevel))
        {
            repeated.measure = multilevel[,1]
            n = length(unique(repeated.measure)) # unique observation: we put every observation of the same "sample" in the either the training or test set
        }
        
        
        #-- define the folds --#
        if (validation ==  "Mfold")
        {
            
            if (nrep > 1) # reinitialise the folds
            folds = folds.input
            
            if (is.null(folds) || !is.numeric(folds) || folds < 2 || folds > n)
            {
                stop("Invalid number of folds.")
            } else {
                M = round(folds)
                if (is.null(multilevel))
                {
                    if(any(class.object == "DA")){
                        temp = stratified.subsampling(Y, folds = M)
                        folds = temp$SAMPLE
                        if(temp$stop > 0 & nrep == 1) # to show only once
                        warning("At least one class is not represented in one fold, which may unbalance the error rate.\n  Consider a number of folds lower than the minimum in table(Y): ", min(table(Y)))
                        rm(temp)
                    } else {
                        folds = suppressWarnings(split(sample(1:n), rep(1:n, length = M)))
                    }
                    
                } else {
                    folds = split(sample(1:n), rep(1:M, length = n)) # needs to have all repeated samples in the same fold
                }
            }
        } else if (validation ==  "loo") {
            folds = split(1:n, rep(1:n, length = n))
            M = n
        }
        
        M = length(folds)
        
        error.sw = matrix(0,nrow = M, ncol = length(test.keepX))
        rownames(error.sw) = paste0("fold",1:M)
        colnames(error.sw) = test.keepX
        # for the last keepX (i) tested, prediction combined for all M folds so as to extract the error rate per class
        # prediction.all = vector(length = nrow(X))
        # in case the test set only includes one sample, it is better to advise the user to
        # perform loocv
        stop.user = FALSE

        # function instead of a loop so we can use lapply and parLapply. Can't manage to put it outside without adding all the arguments
        #result.all=list()
        #save(list=ls(),file="temp22.Rdata")
        
        fonction.j.folds = function(j)#for (j in 1:M)
        {
            if (progressBar ==  TRUE)
            setTxtProgressBar(pb, (M*(nrep-1)+j-1)/(M*nrepeat))
            
            #print(j)
            #---------------------------------------#
            #-- set up leave out samples. ----------#

            omit = which(repeated.measure %in% folds[[j]] == TRUE)
            
            # get training and test set
            X.train = X[-omit, ]
            X.test = X[omit, , drop = FALSE]#matrix(X[omit, ], nrow = length(omit)) #removed to keep the colnames in X.test
            
            if(any(class.object == "DA")){
                Y.train = Y[-omit]
                Y.train.mat = unmap(Y.train)
                q = ncol(Y.train.mat)
                colnames(Y.train.mat) = levels(Y.train)
                Y.test = Y[omit]
            } else {
                Y.train.mat = Y[-omit, , drop = FALSE]
                q = ncol(Y.train.mat)
                Y.test = Y[omit, , drop = FALSE]
            }
            
            #-- set up leave out samples. ----------#
            #---------------------------------------#


            #---------------------------------------#
            #-- near.zero.var ----------------------#
            remove = NULL
            # first remove variables with no variance
            var.train = colVars(X.train, na.rm=TRUE)#apply(X.train, 2, var)
            ind.var = which(var.train == 0)
            if (length(ind.var) > 0)
            {
                remove = c(remove, colnames(X.train)[ind.var])

                X.train = X.train[, -c(ind.var),drop = FALSE]
                X.test = X.test[, -c(ind.var),drop = FALSE]
                
                # reduce choice.keepX and test.keepX if needed
                if (any(choice.keepX > ncol(X.train)))
                choice.keepX[which(choice.keepX>ncol(X.train))] = ncol(X.train)

                # reduce test.keepX if needed
                if (any(test.keepX > ncol(X.train)))
                test.keepX[which(test.keepX>ncol(X.train))] = ncol(X.train)
                
            }
            
            if(near.zero.var == TRUE)
            {
                remove.zero = nearZeroVar(X.train)$Position
                
                if (length(remove.zero) > 0)
                {
                    remove = c(remove, colnames(X.train)[remove.zero])

                    X.train = X.train[, -c(remove.zero),drop = FALSE]
                    X.test = X.test[, -c(remove.zero),drop = FALSE]
                    
                    # reduce choice.keepX and test.keepX if needed
                    if (any(choice.keepX > ncol(X.train)))
                    choice.keepX[which(choice.keepX>ncol(X.train))] = ncol(X.train)
                    
                    # reduce test.keepX if needed
                    if (any(test.keepX > ncol(X.train)))
                    test.keepX[which(test.keepX>ncol(X.train))] = ncol(X.train)
                    
                }
                #print(remove.zero)
            }
            
            #-- near.zero.var ----------------------#
            #---------------------------------------#
            
            
            #------------------------------------------#
            #-- split the NA in training and testing --#
            if(any(misdata))
            {
                if(any(class.object == "DA")){
                    if(length(remove)>0){
                        ind.remove = which(colnames(X) %in% remove)
                        is.na.A.train = is.na.A[-omit, -ind.remove, drop=FALSE]
                        is.na.A.test = list(X=is.na.A[omit, -ind.remove, drop=FALSE])
                    } else {
                        is.na.A.train = is.na.A[-omit,, drop=FALSE]
                        is.na.A.test = list(X=is.na.A[omit,, drop=FALSE])
                    }
                    
                    temp = which(is.na.A.train, arr.ind=TRUE)
                    ind.NA.train = unique(temp[,1])
                    ind.NA.col.train = unique(temp[,2])
                    
                    is.na.A.train = list(X=is.na.A.train, Y=NULL)
                    ind.NA.train = list(X=ind.NA.train, Y=NULL)
                    ind.NA.col.train = list(X=ind.NA.col.train, Y=NULL)
                } else{
                    if(length(remove)>0){
                        ind.remove = which(colnames(X) %in% remove)
                        is.na.A.train = list(X=is.na.A[[1]][omit, -ind.remove, drop=FALSE], Y=is.na.A[[1]][omit,, drop=FALSE])
                        #lapply(is.na.A, function(x){x[-omit,, drop=FALSE]})
                        is.na.A.test = list(X=is.na.A[[1]][omit, -ind.remove, drop=FALSE]) #only for X
                    }else {
                        is.na.A.train = lapply(is.na.A, function(x){x[-omit,, drop=FALSE]})
                        is.na.A.test = list(X=is.na.A[[1]][omit,, drop=FALSE]) #only for X
                        
                    }
                    temp = lapply(is.na.A.train,function(x){which(x,arr.ind=TRUE)})
                    
                    ind.NA.train = lapply(temp,function(x){unique(x[,1])})
                    ind.NA.col.train = lapply(temp,function(x){unique(x[,2])})
                }
                #ind.NA.train = which(apply(is.na.A.train, 1, sum) > 0) # calculated only once
                #ind.NA.test = which(apply(is.na.A.test, 1, sum) > 0) # calculated only once
                
                #ind.NA.col.train = which(apply(is.na.A.train, 2, sum) > 0) # calculated only once
                #ind.NA.col.test = which(apply(is.na.A.test, 2, sum) > 0) # calculated only once
                
            } else {
                is.na.A.train = is.na.A.test =NULL
                ind.NA.train = NULL
                ind.NA.col.train = NULL
            }
            #-- split the NA in training and testing --#
            #------------------------------------------#


            class.comp.j = list()
            if(any(class.object == "DA")){
                prediction.comp.j = array(0, c(length(omit), nlevels(Y), length(test.keepX)), dimnames = list(rownames(X.test), levels(Y), names(test.keepX)))
                
                for(ijk in dist)
                class.comp.j[[ijk]] = matrix(0, nrow = length(omit), ncol = length(test.keepX))# prediction of all samples for each test.keepX and  nrep at comp fixed
            } else{
                prediction.comp.j = array(0, c(length(omit), ncol(Y), length(test.keepX)), dimnames = list(rownames(X.test), colnames(Y), test.keepX))
            }
            
            
            # shape input for `internal_mint.block' (keepA, test.keepA, etc)
            result = suppressWarnings(internal_wrapper.mint(X=X.train, Y=Y.train.mat, study=factor(rep(1,nrow(X.train))), ncomp=ncomp,
            keepX=choice.keepX, keepY=rep(ncol(Y.train.mat), ncomp-1), test.keepX=test.keepX, test.keepY=test.keepY,
            mode="regression", scale=scale, near.zero.var=near.zero.var,
            max.iter=max.iter, logratio="none", DA=TRUE, multilevel=NULL,
            misdata = misdata, is.na.A = is.na.A.train, ind.NA = ind.NA.train,
            ind.NA.col = ind.NA.col.train, all.outputs=FALSE))
            
            # `result' returns loadings and variates for all test.keepX on the ncomp component
            
            # need to find the best keepX/keepY among all the tested models
            
            #---------------------------------------#
            #-- scaling X.test ---------------------#
            
            # we prep the test set for the successive prediction: scale and is.na.newdata
            if (!is.null(attr(result$A[[1]], "scaled:center")))
            X.test = sweep(X.test, 2, STATS = attr(result$A[[1]], "scaled:center"))
            if (scale)
            X.test = sweep(X.test, 2, FUN = "/", STATS = attr(result$A[[1]], "scaled:scale"))
            
            means.Y = matrix(attr(result$A[[2]], "scaled:center"),nrow=nrow(X.test),ncol=q,byrow=TRUE);
            if (scale)
            {sigma.Y = matrix(attr(result$A[[2]], "scaled:scale"),nrow=nrow(X.test),ncol=q,byrow=TRUE)}else{sigma.Y=matrix(1,nrow=nrow(X.test),ncol=q)}
            
            #-- scaling X.test ---------------------#
            #---------------------------------------#


            #-----------------------------------------#
            #-- prediction on X.test for all models --#

            # record prediction results for each test.keepX
            keepA = result$keepA
            test.keepA = keepA[[ncomp]]
            
            #save variates and loadings for all test.keepA
            result.temp = list(variates = result$variates, loadings = result$loadings)

            # creates temporary splsda object to use the predict function
            result$X = result$A$X
            
            # add the "splsda" or "spls" class
            if(any(class.object == "DA")){
                class(result) = c("mixo_splsda","mixo_spls","DA")
                result$ind.mat = result$A$Y
                result$Y = factor(Y.train)
            } else{
                class(result) = c("mixo_spls")
                result$Y = result$A$Y
            }
            result$A = NULL




            
            
            for(i in 1:nrow(test.keepA))
            {
                #print(i)
                
                # only pick the loadings and variates relevant to that test.keepX
                
                names.to.pick = NULL
                if(ncomp>1)
                names.to.pick = unlist(lapply(1:(ncomp-1), function(x){
                    paste(paste0("comp",x),apply(keepA[[x]],1,function(x) paste(x,collapse="_")), sep=":")
                    
                }))
                
                names.to.pick.ncomp = paste(paste0("comp",ncomp),paste(as.numeric(keepA[[ncomp]][i,]),collapse="_"), sep=":")
                names.to.pick = c(names.to.pick, names.to.pick.ncomp)

                #change variates and loadings for each test.keepA
                result$variates = lapply(result.temp$variates, function(x){if(ncol(x)!=ncomp) {x[,colnames(x)%in%names.to.pick, drop=FALSE]}else{x}})
                result$loadings = lapply(result.temp$loadings, function(x){if(ncol(x)!=ncomp) {x[,colnames(x)%in%names.to.pick, drop=FALSE]}else{x}})
                
                # added: record selected features
                if (any(class.object == "mixo_splsda") & length(test.keepX) ==  1) # only done if splsda and if only one test.keepX as not used if more so far
                # note: if plsda, 'features' includes everything: to optimise computational time, we don't evaluate for plsda object
                features.j = selectVar(result, comp = ncomp)$name

                # do the prediction, we are passing to the function some invisible parameters:
                # the scaled newdata and the missing values
                #save(list=ls(),file="temp.Rdata")

                test.predict.sw <- predict.mixo_spls(result, newdata.scale = X.test, dist = dist, misdata.all=misdata[1], is.na.X = is.na.A.train, is.na.newdata = is.na.A.test)
                prediction.comp.j[, , i] =  test.predict.sw$predict[, , ncomp]
                
                if(any(class.object == "DA")){
                    for(ijk in dist)
                    class.comp.j[[ijk]][, i] =  test.predict.sw$class[[ijk]][, ncomp] #levels(Y)[test.predict.sw$class[[ijk]][, ncomp]]
                }
            } # end i
            
            #-- prediction on X.test for all models --#
            #-----------------------------------------#


            return(list(class.comp.j = class.comp.j, prediction.comp.j = prediction.comp.j, features = features.j, omit = omit))
            #result.all[[j]] = list(class.comp.j = class.comp.j,  prediction.comp.j = prediction.comp.j, features = features.j, omit = omit)

        } # end fonction.j.folds
        
        if (parallel == TRUE)
        {
            clusterExport(cl, c("folds","choice.keepX","ncomp"),envir=environment())
            #clusterExport(cl, ls(), envir=environment())
            #print(clusterEvalQ(cl,ls()))
            
            result.all = parLapply(cl, 1: M, fonction.j.folds)
        } else {
            result.all = lapply(1: M, fonction.j.folds)
            
        }
        #---------------------------#
        #--- combine the results ---#

        for(j in 1:M)
        {
            omit = result.all[[j]]$omit
            prediction.comp.j = result.all[[j]]$prediction.comp.j

            prediction.comp[[nrep]][omit, , ] = prediction.comp.j
            
            if(any(class.object == "DA")){
                class.comp.j = result.all[[j]]$class.comp.j
                for(ijk in dist)
                class.comp[[ijk]][omit,nrep, ] = class.comp.j[[ijk]]
                
                if (length(test.keepX) ==  1) # only done if splsda and if only one test.keepX as not used if more so far
                features = c(features, result.all[[j]]$features)
            }
            

        }
        if(!any(class.object == "DA")){
            # create a array, for each keepX: n*q*nrep
            for(i in 1:length(test.keepX))
                prediction.keepX[[i]][,,nrep] = prediction.comp[[nrep]][,,i, drop=FALSE]
        }


        #--- combine the results ---#
        #---------------------------#
        

        if (progressBar ==  TRUE)
        setTxtProgressBar(pb, (M*nrep)/(M*nrepeat))
        
        #---------------------------#
        #----- AUC on the test -----#
        if(auc)
        {
            data=list()
            for (i in 1:length(test.keepX))
            {
                data$outcome = Y
                data$data = prediction.comp[[nrep]][, , i]
                auc.all[[nrep]][, , i] = as.matrix(statauc(data)[[1]]) # [[1]] because [[2]] is the graph
            }
        }
        #----- AUC on the test -----#
        #---------------------------#

    } #end nrep 1:nrepeat
    if(auc) names(auc.all) = paste0("nrep.", 1:nrepeat)
    
    names(prediction.comp) = paste0("nrep.", 1:nrepeat)
    # class.comp[[ijk]] is a matrix containing all prediction for test.keepX, all nrepeat and all distance, at comp fixed
    
    #-------------------------------------------------------------#
    #----- average AUC over the nrepeat, for each test.keepX -----#
    if(auc)
    {
        
        if(nlevels(Y)>2)
        {
            auc.mean.sd =  array(0, c(nlevels(Y),2, length(test.keepX)), dimnames = list(rownames(auc.all[[1]]), c("AUC.mean","AUC.sd"), names(test.keepX)))
        }else{
            auc.mean.sd =  array(0, c(1,2, length(test.keepX)), dimnames = list(rownames(auc.all[[1]]), c("AUC.mean","AUC.sd"), names(test.keepX)))
        }
        
        for(i in 1:length(test.keepX))
        {
            temp = NULL
            for(nrep in 1:nrepeat)
            {
                temp = cbind(temp, auc.all[[nrep]][, 1, i])
            }
            auc.mean.sd[, 1, i] = apply(temp,1,mean)
            auc.mean.sd[, 2, i] = apply(temp,1,sd)
        }
    } else {
        auc.mean.sd = auc.all = NULL
        
    }
    #----- average AUC over the nrepeat, for each test.keepX -----#
    #-------------------------------------------------------------#
    result = list()
    error.mean = error.sd = error.per.class.keepX.opt.comp = keepX.opt = test.keepX.out = test.keepY.out = mat.error.final = choice.keepX.out = choice.keepY.out = list()

    if (any(measure == "overall"))
    {
        for(ijk in dist)
        {
            rownames(class.comp[[ijk]]) = rownames.X
            colnames(class.comp[[ijk]]) = paste0("nrep.", 1:nrepeat)
            dimnames(class.comp[[ijk]])[[3]] = paste0("test.keepX.",names(test.keepX))
            
            #finding the best keepX depending on the error measure: overall or BER
            # classification error for each nrep and each test.keepX: summing over all samples
            error = apply(class.comp[[ijk]],c(3,2),function(x)
            {
                sum(as.character(Y) != x)
            })
            rownames(error) = names(test.keepX)
            colnames(error) = paste0("nrep.",1:nrepeat)
            
            # we want to average the error per keepX over nrepeat and choose the minimum error
            error.mean[[ijk]] = apply(error,1,mean)/length(Y)
            if (!nrepeat ==  1)
            error.sd[[ijk]] = apply(error,1,sd)/length(Y)
            
            mat.error.final[[ijk]] = error/length(Y)  # percentage of misclassification error for each test.keepX (rows) and each nrepeat (columns)
            
            keepX.opt[[ijk]] = which(error.mean[[ijk]] ==  min(error.mean[[ijk]]))[1] # chose the lowest keepX if several minimum
            
            # confusion matrix for keepX.opt
            error.per.class.keepX.opt.comp[[ijk]] = apply(class.comp[[ijk]][, , keepX.opt[[ijk]], drop = FALSE], 2, function(x)
            {
                conf = get.confusion_matrix(truth = factor(Y), predicted = x)
                out = (apply(conf, 1, sum) - diag(conf)) / summary(Y)
            })
            
            rownames(error.per.class.keepX.opt.comp[[ijk]]) = levels(Y)
            colnames(error.per.class.keepX.opt.comp[[ijk]]) = paste0("nrep.", 1:nrepeat)
            
            
            test.keepX.out[[ijk]] = test.keepX[keepX.opt[[ijk]]]
            choice.keepX.out[[ijk]] = c(choice.keepX, test.keepX.out)
            
            result$"overall"$error.rate.mean = error.mean
            if (!nrepeat ==  1)
            result$"overall"$error.rate.sd = error.sd
            
            result$"overall"$confusion = error.per.class.keepX.opt.comp
            result$"overall"$mat.error.rate = mat.error.final
            result$"overall"$keepX.opt = test.keepX.out
        }
    }
    
    if (any(measure ==  "BER"))
    {
        for(ijk in dist)
        {
            rownames(class.comp[[ijk]]) = rownames.X
            colnames(class.comp[[ijk]]) = paste0("nrep.", 1:nrepeat)
            dimnames(class.comp[[ijk]])[[3]] = paste0("test.keepX.",names(test.keepX))
            
            error = apply(class.comp[[ijk]],c(3,2),function(x)
            {
                conf = get.confusion_matrix(truth = factor(Y),predicted = x)
                get.BER(conf)
            })
            rownames(error) = names(test.keepX)
            colnames(error) = paste0("nrep.",1:nrepeat)
            
            # average BER over the nrepeat
            error.mean[[ijk]] = apply(error,1,mean)
            if (!nrepeat ==  1)
            error.sd[[ijk]] = apply(error,1,sd)
            
            mat.error.final[[ijk]] = error  # BER for each test.keepX (rows) and each nrepeat (columns)
            
            keepX.opt[[ijk]] = which(error.mean[[ijk]] ==  min(error.mean[[ijk]]))[1]
            
            # confusion matrix for keepX.opt
            error.per.class.keepX.opt.comp[[ijk]] = apply(class.comp[[ijk]][, , keepX.opt[[ijk]], drop = FALSE], 2, function(x)
            {
                conf = get.confusion_matrix(truth = factor(Y), predicted = x)
                out = (apply(conf, 1, sum) - diag(conf)) / summary(Y)
            })
            
            rownames(error.per.class.keepX.opt.comp[[ijk]]) = levels(Y)
            colnames(error.per.class.keepX.opt.comp[[ijk]]) = paste0("nrep.", 1:nrepeat)
            
            test.keepX.out[[ijk]] = test.keepX[keepX.opt[[ijk]]]
            choice.keepX.out[[ijk]] = c(choice.keepX, test.keepX.out)

            result$"BER"$error.rate.mean = error.mean
            if (!nrepeat ==  1)
            result$"BER"$error.rate.sd = error.sd
            
            result$"BER"$confusion = error.per.class.keepX.opt.comp
            result$"BER"$mat.error.rate = mat.error.final
            result$"BER"$keepX.opt = test.keepX.out
            
        }
        
        
    }
    
    if (any(measure ==  "AUC"))
    {
        # no loop in ijk because AUC does not use any distances
        # still need to put [[1]] to be similar to other distances
        
        #save(list=ls(),file="temp.Rdata")
        # average the AUCs over the class (AUC per keepX)
        error.mean[[1]] = apply(auc.mean.sd,c(3),function(x){mean(x[,1])})

        #choose highest AUC
        keepX.opt[[1]] = which(error.mean[[1]] ==  max(error.mean[[1]]))[1]
        
        # match to the best keepX
        test.keepX.out[[1]]= test.keepX[keepX.opt[[1]]]
        choice.keepX.out[[1]] = c(choice.keepX, test.keepX.out)

        result$"AUC"$error.rate.mean = error.mean
        if (!nrepeat ==  1)
        result$"AUC"$error.rate.sd = NULL
        
        result$"AUC"$confusion = NULL
        result$"AUC"$mat.error.rate = auc.all
        result$"AUC"$keepX.opt = test.keepX.out
    }
    
    
    
    if (any(measure == "MSE"))
    {
        ijk=1 # in case more measure later on
        names(prediction.keepX) = paste0("test.keepX.",test.keepX)
        
        # MSE error for each nrep and each test.keepX: summing over all samples
        varY=apply(Y,2,var)
        error = lapply(prediction.keepX, function(z){apply(z,c(3),function(x)
            {
                temp = (Y-x)^2
                temp2=t(t(temp)/varY)
                temp3=apply(temp2,2, sum)
                temp3
            })})
        mat.error.final[[ijk]] = lapply(error, function(x){x/nrow(Y)})  # percentage of misclassification error for each test.keepX (rows) and each nrepeat (columns)
        
        # we want to average the error per keepX over nrepeat and choose the minimum error
        error.mean[[ijk]] = sapply(error, mean)/nrow(Y)
        if (!nrepeat ==  1)
        error.sd[[ijk]] = sapply(lapply(error, function(x){apply(matrix(x,nrow=ncol(Y)), 1, sd)/nrow(Y)}), mean)
        
        keepX.opt[[ijk]] = which(error.mean[[ijk]] ==  min(error.mean[[ijk]]))[1] # chose the lowest keepX if several minimum
        
        
        test.keepX.out[[ijk]] = test.keepX[keepX.opt[[ijk]]]

        choice.keepX.out[[ijk]] = c(choice.keepX, test.keepX.out)

        result$"MSE"$error.rate.mean = error.mean
        if (!nrepeat ==  1)
        result$"MSE"$error.rate.sd = error.sd
        
        result$"MSE"$mat.error.rate = mat.error.final
        result$"MSE"$keepX.opt = test.keepX.out
    }
    
    if (any(measure == "MAE")) # MAE (Mean Absolute Error: MSE without the square)
    {
        ijk=1 # in case more measure later on
        names(prediction.keepX) = paste0("test.keepX.",test.keepX)

        # MAE error for each nrep and each test.keepX: summing over all samples
        varY=apply(Y,2,var)
        error = lapply(prediction.keepX, function(z){apply(z,c(3),function(x)
            {
                temp = abs(Y-x)
                temp2=t(t(temp)/varY)
                temp3=apply(temp2,2, sum)
                temp3
            })})
        
        mat.error.final[[ijk]] = lapply(error, function(x){x/nrow(Y)})  # percentage of misclassification error for each test.keepX (rows) and each nrepeat (columns)
        
        # we want to average the error per keepX over nrepeat and choose the minimum error
        error.mean[[ijk]] = sapply(error, mean)/nrow(Y)
        if (!nrepeat ==  1)
        error.sd[[ijk]] = sapply(lapply(error, function(x){apply(matrix(x,nrow=ncol(Y)), 1, sd)/nrow(Y)}), mean)
        
        keepX.opt[[ijk]] = which(error.mean[[ijk]] ==  min(error.mean[[ijk]]))[1] # chose the lowest keepX if several minimum
        
        test.keepX.out[[ijk]] = test.keepX[keepX.opt[[ijk]]]
        
        choice.keepX.out[[ijk]] = c(choice.keepX, test.keepX.out)

        result$"MAE"$error.rate.mean = error.mean
        if (!nrepeat ==  1)
        result$"MAE"$error.rate.sd = error.sd
        
        result$"MAE"$mat.error.rate = mat.error.final
        result$"MAE"$keepX.opt = test.keepX.out
    }
    
    if (any(measure == "Bias")) # Bias (average of the differences)
    {
        ijk=1 # in case more measure later on
        names(prediction.keepX) = paste0("test.keepX.",test.keepX)

        # Bias error for each nrep and each test.keepX: summing over all samples
        varY=apply(Y,2,var)
        error = lapply(prediction.keepX, function(z){apply(z,c(3),function(x)
            {
                temp = (Y-x)
                temp2=t(t(temp)/varY)
                temp3=abs(apply(temp2,2, sum)) # absolute value of the bias
                temp3
            })})
        
        mat.error.final[[ijk]] = lapply(error, function(x){x/nrow(Y)})  # percentage of misclassification error for each test.keepX (rows) and each nrepeat (columns)
        
        # we want to average the error per keepX over nrepeat and choose the minimum error
        error.mean[[ijk]] = sapply(error, mean)/nrow(Y)
        if (!nrepeat ==  1)
        error.sd[[ijk]] = sapply(lapply(error, function(x){apply(matrix(x,nrow=ncol(Y)), 1, sd)/nrow(Y)}), mean)
        
        
        keepX.opt[[ijk]] = which(error.mean[[ijk]] ==  min(error.mean[[ijk]]))[1] # chose the lowest keepX if several minimum
        
        
        test.keepX.out[[ijk]] = test.keepX[keepX.opt[[ijk]]]
        
        choice.keepX.out[[ijk]] = c(choice.keepX, test.keepX.out)

        result$"Bias"$error.rate.mean = error.mean
        if (!nrepeat ==  1)
        result$"Bias"$error.rate.sd = error.sd
        
        result$"Bias"$mat.error.rate = mat.error.final
        result$"Bias"$keepX.opt = test.keepX.out
    }
    
    if (any(measure == "R2")) # R2 (square of the correlation of the truth and the predicted values, averaged over the columns of Y)
    {
        ijk=1 # in case more measure later on
        names(prediction.keepX) = paste0("test.keepX.",test.keepX)

        # R2 error for each test.keepX (list),each Y (rows) and each nrepeat (columns)
        error = lapply(prediction.keepX, function(z){apply(z,c(3),function(x)
            {
                temp = diag(cor(Y,x))^2
                temp
            })}) # list for each test.keepX of nrow=ncol(Y) and ncol=nrepeat
        
        mat.error.final[[ijk]] = error
        
        # we want to average the error over the nrepeat and the columns of Y
        error.mean[[ijk]] = sapply(error, mean)
        if (!nrepeat ==  1)
        error.sd[[ijk]] = sapply(lapply(error, function(x){apply(matrix(x,nrow=ncol(Y)), 1, sd)}), mean)
        
        keepX.opt[[ijk]] = which(error.mean[[ijk]] ==  max(error.mean[[ijk]]))[1] # chose the lowest keepX if several minimum
        
        
        test.keepX.out[[ijk]] = test.keepX[keepX.opt[[ijk]]]
        
        choice.keepX.out[[ijk]] = c(choice.keepX, test.keepX.out)

        result$"R2"$error.rate.mean = error.mean
        if (!nrepeat ==  1)
        result$"R2"$error.rate.sd = error.sd
        
        result$"R2"$mat.error.rate = mat.error.final
        result$"R2"$keepX.opt = test.keepX.out
    }


    result$prediction.comp = prediction.comp
    if(any(class.object == "DA")){
        if(auc){
            result$auc = auc.mean.sd
            result$auc.all = auc.all
        }
        result$class.comp = class.comp
        if (length(test.keepX) ==  1)
        result$features$stable = sort(table(as.factor(features))/M/nrepeat, decreasing = TRUE)
    }
    return(result)
}
