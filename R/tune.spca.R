tune.spca <- function(X, ncomp, nrepeat, kfold, grid.keepX, center = TRUE, scale = TRUE) {
    cor.pred.per.dim <- list()
    keepX.dim <- NULL
    result <- list()
    for(ncomp in seq_len(ncomp)) {
        # 1 - a foreach list for each keepX value tested
        cor.pred.repeat.keepX.dim = foreach(keepX.value = as.list(grid.keepX),.combine=cbind) %do% {
            cat('KeepX = ', keepX.value, '\n')  # to remove in the final function
            
            # 2 - a foreach list for repeated CV
            cor.pred.repeat = foreach(j = as.list(c(seq_len(nrepeat))),.combine=cbind) %do% {
                folds = split(sample(seq_len(nrow(X))),seq_len(kfold))
                
                # 3 -  a foreach list for k-fold CV
                #  if small n, then need to define the fold based on boostrapping (with replacement) as we need to center / scale the data for prediction
                cor.pred = foreach(i = folds,.combine=cbind) %do% {
                    # determine matrix without the fold and with the fold
                    X.train = X[-i,]  # could rename as train
                    X.test = X[i,]  # used for prediction, could rename as test
                    X.test = scale(X.test, center = center, scale = scale) # used for deflation
                    
                    # ---- run sPCA ------------ #
                    # spca on the data minus the subsample
                    suppressWarnings({
                        spca.res.sub = mixOmics::spca(X.train, ncomp = ncomp, keepX = c(keepX.dim, keepX.value), center = center, scale = scale)
                        # spca on all data 
                        spca.res.full = mixOmics::spca(X, ncomp = ncomp, keepX = c(keepX.dim, keepX.value), center = center, scale = scale)
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
                    cor.comp = abs(cor(t.comp.sub, spca.res.full$variates$X[,ncomp][i], use = 'pairwise.complete.obs'))
                    # need to flip the sign in one of the comp to calculate the RSS
                    #RSS.comp = sum((c(sign(cor.comp))*(t.comp.sub) - spca.res.full$variates$X[,comp][i])^2)
                    # also look at correlation between loading vectors
                    #cor.loading = abs(cor(spca.res.full$loadings$X[,comp], spca.res.sub$loadings$X[,comp]))
                    
                    # return
                    return(cor.comp)  #, cor.loading, RSS.comp))
                }  # end foreach 3 (on folds), get the cor(pred component, comp) per fold for a given comp
                
                # output is correlation between reconstructed components cor.comp on the left out samples and if we had those samples in PCA, as well as correlation between loading vectors cor.loading just to check both trends
                rownames(cor.pred) = c('cor.comp')  #, 'cor.loading', 'RSS.comp')
                # average correlations across folds
                return(abs(apply(cor.pred, 1, mean)))
            } # end foreach 2 (on repeats), get the cor(pred component, comp) averaged across folds for a given comp
            
            # average correlation across repeats 
            return(abs(apply(cor.pred.repeat, 1, mean)))
        } # end foreach 1 (on keepX), get the cor(pred component, comp) averaged across repeats for a given comp and each keepX
        
        # # correlation for each keepX value
        cor.pred.per.dim[[ncomp]]  = cor.pred.repeat.keepX.dim
        
        # working out the max based on cor.comp and append to keepX.dim
        keepX.dim = c(keepX.dim, grid.keepX[which.max(cor.pred.repeat.keepX.dim['cor.comp',])])
        result[[paste0("comp", ncomp)]] <- keepX.dim  # I think this output is irrelevant, as it is combining across comp per column (consider the last column)
    }
    result
}
