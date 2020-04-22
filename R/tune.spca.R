tune.spca <- function(X, ncomp, nrepeat, kfold, grid.keepX, center = TRUE, scale = TRUE) {
    keepX.opt <- NULL
    cor.df <- data.frame(matrix(ncol = nrepeat, nrow = length(grid.keepX), 
                                dimnames = list(
                                    paste0('keepX_', grid.keepX),
                                    paste0('repeat_', seq_len(nrepeat)))))
    cor.df.list <- .name_list(char = paste0('comp', seq_len(ncomp)))
    cor.df.list <- lapply(cor.df.list, function(x) cor.df)
    
    for(ncomp in seq_len(ncomp)) {
        # 1 - a foreach list for each keepX value tested
        for (keepX_i in seq_along(grid.keepX)) {
            keepX.value <- grid.keepX[keepX_i]
            cat('KeepX = ', keepX.value, '\n')  # to remove in the final function
            
            # 2 - a foreach list for repeated CV
            for (j in seq_len(nrepeat)) {
                folds = split(sample(seq_len(nrow(X))),seq_len(kfold))
                
                # 3 -  a foreach list for k-fold CV
                #  if small n, then need to define the fold based on boostrapping (with replacement) as we need to center / scale the data for prediction
                cor.pred = sapply(folds, function(test.fold.inds){
                    # determine matrix without the fold and with the fold
                    X.train = X[-test.fold.inds,]  # could rename as train
                    X.test = X[test.fold.inds,]  # used for prediction, could rename as test
                    X.test = scale(X.test, center = center, scale = scale) # used for deflation
                    
                    # ---- run sPCA ------------ #
                    # spca on the data minus the subsample
                    suppressWarnings({
                        spca.res.sub = mixOmics::spca(X.train, ncomp = ncomp, keepX = c(keepX.opt, keepX.value), center = center, scale = scale)
                        # spca on all data 
                        spca.res.full = mixOmics::spca(X, ncomp = ncomp, keepX = c(keepX.opt, keepX.value), center = center, scale = scale)
                    })
                    # ---- deflation on X with only the fold left out------------ #
                    for(k in seq_len(ncomp)){ # loop to calculate deflated matrix and predicted comp
                        # calculate the predicted comp on the fold left out
                        # calculate reg coeff, then deflate
                        if(k == 1){
                            t.comp.sub = X.test %*% spca.res.sub$loadings$X[,k]
                        }else{ # calculate deflation beyond comp 1
                            # recalculate the loading vector (here c.sub) on the test set (perhaps we could do this instead on the training set by extracting from spca.res.sub$loadings$X[,k]?)
                            c.sub = crossprod(X.test, t.comp.sub) / drop(crossprod(t.comp.sub)) 
                            X.test = X.test - t.comp.sub %*% t(c.sub) 
                            # update predicted comp based on deflated matrix
                            t.comp.sub = X.test %*% spca.res.sub$loadings$X[,k]
                        }
                    }
                    
                    # ---- calculate predicted component and compare with component from sPCA on full data on the left out set------------ #
                    # cor with the component on the full data, abs value
                    cor.comp = abs(cor(t.comp.sub, spca.res.full$variates$X[,ncomp][test.fold.inds], use = 'pairwise.complete.obs'))
                    # need to flip the sign in one of the comp to calculate the RSS
                    #RSS.comp = sum((c(sign(cor.comp))*(t.comp.sub) - spca.res.full$variates$X[,comp][i])^2)
                    # also look at correlation between loading vectors
                    #cor.loading = abs(cor(spca.res.full$loadings$X[,comp], spca.res.sub$loadings$X[,comp]))
                    
                    # return
                    return(cor.comp)  #, cor.loading, RSS.comp))
                })  # end foreach 3 (on folds), get the cor(pred component, comp) per fold for a given comp
                
                # output is correlation between reconstructed components cor.comp on the left out samples and if we had those samples in PCA, as well as correlation between loading vectors cor.loading just to check both trends
                # rownames(cor.pred) = c('cor.comp')  #, 'cor.loading', 'RSS.comp')
                # average correlations across folds
                cor.df.list[[ncomp]][keepX_i,j] <- mean(cor.pred, na.rm = TRUE)
            } # end foreach 2 (on repeats), get the cor(pred component, comp) averaged across folds for a given comp
        } # end foreach 1 (on keepX), get the cor(pred component, comp) averaged across repeats for a given comp and each keepX
        keepX.opt.comp.ind <-  t.test.process(t(cor.df.list[[ncomp]]), alpha = 0.05, alternative = 'less')
        keepX.opt.comp <- grid.keepX[keepX.opt.comp.ind]
        keepX.opt <- c(keepX.opt, keepX.opt.comp)
    }
    choice.keepX <- lapply(cor.df.list, function(x) {
        keepX.ind <- t.test.process(t(x), alpha = 0.05, alternative = 'less')
        return(grid.keepX[keepX.ind])
    })
    
    cor.comp = lapply(cor.df.list, function(df){
        return(data.frame(cor.mean = apply(df, 1, mean), cor.sd = apply(df, 1, sd)))
    })
    result <- list(choice.keepX = choice.keepX,
                   cor.comp = cor.comp)
    result
}
