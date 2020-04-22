tune.spca <- function(X, ncomp, nrepeat, kfold, grid.keepX, center = TRUE, scale = TRUE) {
    ## optimal keepX for all components
    keepX.opt <- NULL
    ## a data.frame to store correlations for each keepX at each repeat
    cor.df <- data.frame(matrix(ncol = nrepeat, nrow = length(grid.keepX), 
                                dimnames = list(
                                    paste0('keepX_', grid.keepX),
                                    paste0('repeat_', seq_len(nrepeat)))))
    ## a list of cor.df for each component
    cor.df.list <- .name_list(char = paste0('comp', seq_len(ncomp)))
    cor.df.list <- lapply(cor.df.list, function(x) cor.df)
    
    ## ------ component loop
    for(ncomp in seq_len(ncomp)) {
        for (keepX_i in seq_along(grid.keepX)) {
            keepX.value <- grid.keepX[keepX_i]
            cat('KeepX = ', keepX.value, '\n')  # to remove in the final function
            
            ## ------ repeated cv
            for (j in seq_len(nrepeat)) {
                folds = split(sample(seq_len(nrow(X))),seq_len(kfold))
               ## ------ mean cor for CV
                cor.pred = sapply(folds, function(test.fold.inds){
                    ## split data to train/test
                    X.train = X[-test.fold.inds,]
                    X.test = X[test.fold.inds,]
                    X.test = scale(X.test, center = center, scale = scale)
                    
                    # ---- run sPCA 
                    suppressWarnings({
                        ## train
                        spca.train = mixOmics::spca(X.train, ncomp = ncomp, keepX = c(keepX.opt, keepX.value), center = center, scale = scale)
                        # full data
                        spca.full = mixOmics::spca(X, ncomp = ncomp, keepX = c(keepX.opt, keepX.value), center = center, scale = scale)
                    })
                    # ---- deflation on X with only the fold left out
                    for(k in seq_len(ncomp)){ 
                        # loop to calculate deflated matrix and predicted comp
                        # calculate the predicted comp on the fold left out
                        # calculate reg coeff, then deflate
                        if(k == 1){
                            t.comp.sub = X.test %*% spca.train$loadings$X[,k]
                        } else{
                            # calculate deflation beyond comp 1
                            # recalculate the loading vector (here c.sub) on the test set (perhaps we could do this instead on the training set by extracting from spca.train$loadings$X[,k]?)
                            c.sub = crossprod(X.test, t.comp.sub) / drop(crossprod(t.comp.sub)) 
                            X.test = X.test - t.comp.sub %*% t(c.sub) 
                            # update predicted comp based on deflated matrix
                            t.comp.sub = X.test %*% spca.train$loadings$X[,k]
                        }
                    }
                    
                    # calculate predicted component and compare with component from sPCA on full data on the left out set
                    # cor with the component on the full data, abs value
                    cor.comp = abs(cor(t.comp.sub, spca.full$variates$X[,ncomp][test.fold.inds], use = 'pairwise.complete.obs'))
                    # need to flip the sign in one of the comp to calculate the RSS
                    #RSS.comp = sum((c(sign(cor.comp))*(t.comp.sub) - spca.full$variates$X[,comp][i])^2)
                    # also look at correlation between loading vectors
                    # cor.loading = abs(cor(spca.full$loadings$X[,comp], spca.train$loadings$X[,comp]))
                    # get the cor(pred component, comp) per fold for a given comp
                    return(cor.comp)
                })  # end CV, get the cor(pred component, comp) per fold for a given comp
                
                # output is correlation between reconstructed components cor.comp on the left out samples and if we had those samples in PCA,
                # average correlations across folds for the repeat
                cor.df.list[[ncomp]][keepX_i,j] <- mean(cor.pred, na.rm = TRUE)
            } # end repeats loop
        } # end keepX loop
        ## use a one-sided t.test using repeat correlations to assess if addition of keepX improved the correlation
        ## and get the index of optimum keepX
        keepX.opt.comp.ind <-  t.test.process(t(cor.df.list[[ncomp]]), alpha = 0.05, alternative = 'less')
        keepX.opt.comp <- grid.keepX[keepX.opt.comp.ind]
        ## update keepX.optimum for next comp
        keepX.opt <- c(keepX.opt, keepX.opt.comp)
    }
    choice.keepX <- keepX.opt
    names(choice.keepX) <- paste0('comp', seq_len(ncomp))
    
    cor.comp = lapply(cor.df.list, function(df){
        return(data.frame(cor.mean = apply(df, 1, mean), cor.sd = apply(df, 1, sd)))
    })
    result <- list(choice.keepX = choice.keepX,
                   cor.comp = cor.comp)
    result
}
